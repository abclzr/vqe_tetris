from synthesis_sd import *
from utils.hardware import *
from synthesis_FT import assign_time_parameter
from functools import partial
from utils.scheduler import Scheduler
import random
import pdb

from qiskit import QuantumCircuit

def pauli_single_gates(qc, pauli_map, ps, left=True):
    if left == True:
        for i in range(len(ps)):
            if ps[i] == 'X':
                qc.u(np.pi/2, 0, np.pi, pauli_map[i])
            elif ps[i] == 'Y':
                qc.u(np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
    else:
        for i in range(len(ps)):
            if ps[i] == 'X':
                qc.u(np.pi/2, 0, np.pi, pauli_map[i])
            elif ps[i] == 'Y':
                qc.u(-np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])

def synthesis_initial(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    assign_time_parameter(pauli_layers, 1)
    lnq = len(pauli_layers[0][0][0]) # logical qubits
    if graph == None:
        G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
        graph = pGraph(G, C)
    if pauli_map == None:
        pauli_map = dummy_qubit_mapping(graph, lnq)
    else:
        add_pauli_map(graph, pauli_map)
    pnq = len(graph) # physical qubits
    if qc == None:
        qc = QuantumCircuit(pnq)
    return pauli_map, graph, qc

def synthesis(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    pauli_map, graph, qc = synthesis_initial(pauli_layers, pauli_map, graph, qc, arch)
    scheduler = Scheduler(pauli_map, graph, qc)
    n_qubits = len(pauli_layers[0][0][0].ps)
    for blocks in pauli_layers:
        for block in blocks:
            level = [-1 for i in range(n_qubits)]
            prior = ['' for i in range(n_qubits)]
            # level 0: always I
            # level 1: always X, Y or Z
            # level 2: not always the same pauli
            for pauli_string in block:
                for wire, pauli_op in enumerate(pauli_string.ps):
                    if prior[wire] == '':
                        if pauli_op == 'I':
                            level[wire] = 0
                            prior[wire] = 'I'
                        elif pauli_op == 'X' or pauli_op == 'Y' or pauli_op == 'Z':
                            level[wire] = 1
                            prior[wire] = pauli_op
                        else:
                            raise Exception('None I, X, Y or Z character in ' + pauli_string.ps)
                        continue
                    
                    if level[wire] == 2:
                        continue
                    
                    if pauli_op == 'I':
                        if level[wire] == 1:
                            level[wire] = 2
                    elif pauli_op == 'X' or pauli_op == 'Y' or pauli_op == 'Z':
                        if level[wire] == 0:
                            level[wire] = 2
                        elif level[wire] == 1 and prior[wire] != pauli_op:
                            level[wire] = 2
                    else:
                        raise Exception('None I, X, Y or Z character in ' + pauli_string.ps)
                        
            # assign level 1 qubits as flower_head
            # assign level 2 qubits as stalk
            flower_head = []
            stalk = []
            for i, l in enumerate(level):
                if l == 1:
                    flower_head.append(i)
                elif l == 2:
                    stalk.append(i)
            
            # do MST algorithm 3 times:
            # 1st time: connect flower_head
            # 2nd time: connect stalk
            # 3rd time: connect flower_head and stalk by only one edge
            
            scheduler.MST_init(n_qubits)
            
            mst_edges1 = scheduler.MST(flower_head, edges=[(flower_head[i], flower_head[j])\
                                                for i in range(len(flower_head))\
                                                    for j in range(i + 1, len(flower_head))])
            mst_edges2 = scheduler.MST(stalk, edges=[(stalk[i], stalk[j])\
                                                for i in range(len(stalk))\
                                                    for j in range(i + 1, len(stalk))])
            mst_edges3 = scheduler.MST(flower_head + stalk, edges=[(f, s) for f in flower_head for s in stalk])
            
            # decide root:
            # pick a non-I node in stalk
            find_root = True
            # scheduler.Tree_init(mst_edges1 + mst_edges2 + mst_edges3, root)
            
            # execution of each pauli string
            for pauli_string in block:
                if find_root == True:
                    root = None
                    for i in stalk:
                        if pauli_string.ps[i] != 'I':
                            root = i
                            break
                    scheduler.Tree_init(mst_edges2, root)
                
                # the left side of a pauli string circuit
                scheduler.enable_cancel = True
                for i in flower_head + stalk:
                    pauli = pauli_string.ps[i]
                    if pauli == 'I' or pauli == 'Z':
                        pass
                    elif pauli == 'X':
                        scheduler.add_instruction('Logical_left_X', i)
                    elif pauli == 'Y':
                        scheduler.add_instruction('Logical_left_Y', i)
                    else:
                        raise Exception('Illegal pauli operator: ' + pauli)
                
                for node in scheduler.tree.node_list:
                    if node.parent != -1:
                        scheduler.add_instruction('Logical_CNOT', (node.idx, node.parent))
                    else:
                        scheduler.add_instruction('Logical_RZ', node.idx)
                
                scheduler.clear_uncompiled_logical_instructions()
                # the right side of a pauli string circuit
                scheduler.enable_cancel = False
                for node in reversed(scheduler.tree.node_list):
                    if node.parent != -1:
                        scheduler.add_instruction('Logical_CNOT', (node.idx, node.parent))
                
                for i in reversed(flower_head + stalk):
                    pauli = pauli_string.ps[i]
                    if pauli == 'I' or pauli == 'Z':
                        pass
                    elif pauli == 'X':
                        scheduler.add_instruction('Logical_right_X', i)
                    elif pauli == 'Y':
                        scheduler.add_instruction('Logical_right_Y', i)
                    else:
                        raise Exception('Illegal pauli operator: ' + pauli)

    scheduler.clear_uncompiled_logical_instructions()
    # pdb.set_trace()
    return scheduler.qc

def dummy_synthesis(pauli_layers, pauli_map=None, graph=None, qc=None, arch='manhattan'):
    lnq = len(pauli_layers[0][0][0]) # logical qubits
    if graph == None:
        G, C = load_graph(arch, dist_comp=True) # G is adj, C is dist
        graph = pGraph(G, C)
    if pauli_map == None:
        pauli_map = dummy_qubit_mapping(graph, lnq)
    else:
        add_pauli_map(graph, pauli_map)
    pnq = len(graph) # physical qubits
    if qc == None:
        qc = QuantumCircuit(pnq)
    for i1 in pauli_layers: # i1 is layer of blocks
        for i2 in i1: # i2 is block of pauli strings
            for i3 in i2:  # i3 is pauli string
                cns = ps2nodes(i3.ps)
                pauli_single_gates(qc, pauli_map, i3.ps, left=True)
                # for i in cns:
                #     if i3.ps[i] == 'X':
                #         qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                #         # qc.h(pauli_map[i])
                #     elif i3.ps[i] == 'Y':
                #         qc.u(np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
                for i4 in range(len(cns)-1):
                    dummy_local_move(qc, graph, pauli_map, pauli_map[cns[i4]], pauli_map[cns[i4+1]])
                if len(cns) >= 1:
                    qc.rz(i3.real, pauli_map[cns[-1]])
                for i4 in range(len(cns)-1, 0, -1):
                    dummy_local_move(qc, graph, pauli_map, pauli_map[cns[i4-1]], pauli_map[cns[i4]])
                pauli_single_gates(qc, pauli_map, i3.ps, left=False)
                # for i in cns:
                #     if i3.ps[i] == 'X':
                #         # qc.h(pauli_map[i])
                #         qc.u(np.pi/2, 0, np.pi, pauli_map[i])
                #     elif i3.ps[i] == 'Y':
                #         # Y = 1/sqrt{2} [[1, i],[i, 1]]
                #         qc.u(-np.pi/2, -np.pi/2, np.pi/2, pauli_map[i])
    return qc

def qiskit_synthesis(ps_layers, coupling_map=None, arch='manhattan', initial_layout=None, time_parameter=1):
    from qiskit.aqua.operators.legacy import evolution_instruction
    from qiskit.quantum_info import Pauli
    from qiskit import transpile
    # in evolution_instruction, 1.0 means \pi, so, we need to assign time parameter other than 1.0
    # in assign_time_parameter, we assign time_parameter/3.14 to each pauli string.
    # assign_time_parameter(ps_layers, 1)
    nq = len(ps_layers[0][0][0])
    qc = QuantumCircuit(nq)
    # psl = []
    for i in ps_layers:
        for j in i:
            for k in j:
                qc.append(evolution_instruction([[1, Pauli.from_label(k.ps)]], 1, 1), qc.qubits)
    if coupling_map == None:
        coupling_map = load_coupling_map(arch)
    if initial_layout != None:
        return transpile(qc, basis_gates=['u', 'cx'], initial_layout=initial_layout, coupling_map=coupling_map, optimization_level=0)
    else:
        return transpile(qc, basis_gates=['u', 'cx'], coupling_map=coupling_map, optimization_level=0)