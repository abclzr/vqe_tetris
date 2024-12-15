import sys, time
import sys
sys.path.append('../core/')
import argparse
from utils.parallel_bl import *
from qiskit import QuantumCircuit, transpile
from tools import *
from arch import *
import time, sys, os
from t_arch import *
import pdb
import random
from utils.synthesis_lookahead import synthesis_lookahead
from utils.synthesis_lookahead_bfs import synthesis_lookahead_bfs
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
import ast
from benchmark.mypauli import pauliString
from qiskit.quantum_info import Pauli, Statevector
from qiskit.circuit.library import PauliEvolutionGate

old_cwd = os.getcwd()
ctime = time.time

def load_oplist(mapper_name, mole_name):
    if mapper_name == 'random':
        fth = os.path.join('../core/data', 'random', f'random_{mole_name}.pickle')
    else:
        fth = os.path.join('../core/data', mapper_name, mole_name + '_UCCSD.pickle')
    with open(fth, 'rb') as f:
        entry = pickle.load(f)
    return entry

# Tetris Mahattan device method
def Tetris_lookahead_melbourne(parr, use_bridge, swap_coefficient=3, k=10):
    print('Tetris passes, Our schedule, Our synthesis, melbourne', flush=True)
    lnq = len(parr[0][0])
    length = lnq // 2 # `length' is a hyperparameter, and can be adjusted for best performance. Here we keep `length' fixed for simplicity.
    coup = load_coupling_map('melbourne')
    t0 = ctime()
    # a2 = gate_count_oriented_scheduling(parr)#, length=length, maxiter=30)
    # a2 = [block for blocks in a2 for block in blocks]
    a2 = parr
    qc, metrics, sorted_pauli_blocks = synthesis_lookahead_bfs(a2, arch='melbourne', use_bridge=use_bridge, swap_coefficient=swap_coefficient, k=k)
    print('Tetris, Time costed:', ctime()-t0, flush=True)
    print('qubit mapping: ', metrics['mapping'])

    return metrics, sorted_pauli_blocks, qc

def generate_random_pauli_string(n_qubits):
    """Generate a random Pauli string for n_qubits.
    
    Args:
        n_qubits (int): Number of qubits.
    
    Returns:
        str: A random Pauli string.
    """
    pauli_operators = ['I', 'X', 'Y', 'Z']
    return ''.join(random.choice(pauli_operators) for _ in range(n_qubits))


if __name__ == '__main__':
    k = 1
    
    #################################
    # Generating the Pauli strings. #
    #################################
    driver = PySCFDriver(atom="Li .0 .0 .0; H .0 .0 1.3", basis='sto3g')
    mapper = JordanWignerMapper()
    mole = 'LiH'
    mapper_name = 'jordan_wigner'

    problem = driver.run()

    ansatz = UCCSD(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper,
        initial_state=HartreeFock(
            problem.num_spatial_orbitals,
            problem.num_particles,
            mapper,
        ),
    )
    
    pauli_blocks = []
    for pauli_list in ansatz._operators:
        # block = ast.literal_eval(pauli_list._primitive._pauli_list.__str__())
        block = ast.literal_eval(pauli_list.paulis.__str__())
        block = [pauliString(ps) for ps in block]
        pauli_blocks.append(block)

    print("##################Tetris Compilation##################")
    print('UCCSD:', mole)
    metrics, sorted_pauli_blocks, tetris_circuit = Tetris_lookahead_melbourne(pauli_blocks, use_bridge=False)
    print(sorted_pauli_blocks)
    mapping = metrics['mapping']
    n_qubits = len(sorted_pauli_blocks[0][0])
    
    print('##################Correctness Checking##################')
    correct_circuit = QuantumCircuit(n_qubits)
    for block in sorted_pauli_blocks:
        for pauli_string in block:
            evo = PauliEvolutionGate(Pauli(pauli_string.ps), time=.5)
            # The leftmost pauli letter should be applied to the rightmost qargs.
            correct_circuit.append(evo, qargs=list(range(n_qubits))[::-1])
    correct_statevec = Statevector.from_instruction(correct_circuit)

    tetris_statevec = Statevector.from_instruction(tetris_circuit)
    num_trials = 100
    for _ in range(num_trials):
        observable = generate_random_pauli_string(n_qubits)
        value1 = tetris_statevec.expectation_value(Pauli(observable), qargs=mapping[:n_qubits])
        if abs(value1) < 1e-10:
            continue
        value2 = correct_statevec.expectation_value(Pauli(observable), qargs=list(range(n_qubits)))
        print(f"Measure on observable: {observable}")
        print(f"Correct: {value2}")
        print(f"Tetris: {value1}")
        