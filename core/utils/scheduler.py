from utils.hardware import pNode, pGraph
from utils.mst import UnionFind, kruskal_mst
from utils.floyd import floyd_warshall
from utils.tree import Tree

import numpy as np
import pdb


class Scheduler:
    def __init__(self, pauli_map, graph, qc):
        self.pauli_map = pauli_map
        self.graph = graph
        self.qc = qc
        
        self.reverse_pauli_map = [-1 for i in self.graph.data]
        for i, j in enumerate(self.pauli_map):
            # logical qubit i mapped to physical qubit j
            self.reverse_pauli_map[j] = i
        
        # apply floyd_warshall algorithm to calculate the distance
        self.distance = floyd_warshall(self.graph.G)
        self.tree = None
        
        # enable_cancel == True means you can only record the instructions not do any qubit routing
        self.enable_cancel = True
        
        self.instruction_list = []
        self.not_compiled_pointer = 0
    
    def physical_swap(self, physical_i, physical_j):
        assert self.graph.G[physical_i, physical_j] == 1
        logical_i, logical_j = self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j]
        self.pauli_map[logical_i], self.pauli_map[logical_j] = self.pauli_map[logical_j], self.pauli_map[logical_i]
        self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j] = \
            self.reverse_pauli_map[physical_j], self.reverse_pauli_map[physical_i]
    
    def MST_init(self, n_nodes):
        self.union_find = UnionFind(n_nodes)
    
    def MST(self, nodes, edges):
        edges = sorted(edges, key=lambda x: self.distance[self.pauli_map[x[0]]][self.pauli_map[x[1]]])
        mst_edges = kruskal_mst(edges, self.union_find, self.distance)
        return mst_edges
    
    def Tree_init(self, edges, root):
        self.tree = Tree(edges, root)
    
    def shortest_path(self, u, v):
        path = []
        p = u
        while p != v:
            for t in self.graph.data[p].adj:
                if self.distance[p][t] + self.distance[t][v] == self.distance[p][v]:
                    path.append((p, t))
                    p = t
                    break
        return path
    
    def clear_uncompiled_logical_instructions(self):
        for i in range(self.not_compiled_pointer, len(self.instruction_list)):
            (instruction, data) = self.instruction_list[i]
            if instruction.startswith('Logical_left_X'):
                self.qc.u(np.pi/2, 0, np.pi, self.pauli_map[data])
            elif instruction.startswith('Logical_left_Y'):
                self.qc.u(np.pi/2, -np.pi/2, np.pi/2, self.pauli_map[data])
            elif instruction.startswith('Logical_CNOT'):
                u, v = data
                p_u, p_v = self.pauli_map[u], self.pauli_map[v]
                path = self.shortest_path(p_u, p_v)
                for u, v in path[:-1]:
                    self.physical_swap(u, v)
                    self.qc.swap(u, v)
                self.qc.cx(path[-1][0], path[-1][1])
            elif instruction.startswith('Logical_RZ'):
                self.qc.rz(1, self.pauli_map[data])
            elif instruction.startswith('Logical_right_X'):
                self.qc.u(np.pi/2, 0, np.pi, self.pauli_map[data])
            elif instruction.startswith('Logical_right_Y'):
                self.qc.u(-np.pi/2, -np.pi/2, np.pi/2, self.pauli_map[data])
            else:
                raise Exception('Illegal instruction: ' + instruction)
        
        self.not_compiled_pointer = len(self.instruction_list)
    
    def add_instruction(self, instruction, data):
        if self.enable_cancel == False:
            self.instruction_list.append((instruction, data))
            return

        # check if there is any cancellation
        tmp = len(self.instruction_list) - 1
        if instruction.startswith('Logical_CNOT'):
            set_a = set(data)
        else:
            set_a = set([data])
        while tmp >= self.not_compiled_pointer:
            if self.instruction_list[tmp][0].startswith('Logical_CNOT'):
                set_b = set(self.instruction_list[tmp][1])
                if len(set_a & set_b) != 0:
                    break
            else:
                set_b = set([self.instruction_list[tmp][1]])
                if len(set_a & set_b) != 0:
                    break
            tmp = tmp - 1
        
        # tmp is the closest logical operators that not compiled and might be cancelled 
        if tmp >= self.not_compiled_pointer:
            if instruction.startswith('Logical_left_X'):
                if self.instruction_list[tmp][0] == 'Logical_right_X':
                    del self.instruction_list[tmp]
                    return
            elif instruction.startswith('Logical_left_Y'):
                if self.instruction_list[tmp][0] == 'Logical_right_Y':
                    del self.instruction_list[tmp]
                    return
            elif instruction.startswith('Logical_CNOT'):
                if self.instruction_list[tmp][0] == 'Logical_CNOT':
                    if self.instruction_list[tmp][1] == data:
                        del self.instruction_list[tmp]
                        return
            elif instruction.startswith('Logical_RZ'):
                pass
            elif instruction.startswith('Logical_right_X'):
                pass
            elif instruction.startswith('Logical_right_Y'):
                pass
            else:
                raise Exception('Illegal instruction: ' + instruction)
        
        self.instruction_list.append((instruction, data))
