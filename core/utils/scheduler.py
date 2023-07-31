from utils.hardware import pNode, pGraph
import pdb
from enum import ENUM

class IR_TYPE(ENUM):
    Logical_1q = 1
    Logical_CNOT = 2
    Physical_1q = 3
    Physical_CNOT = 4

class Scheduler:
    def __init__(self, pauli_map, graph, qc):
        self.pauli_map = pauli_map
        self.graph = graph
        self.qc = qc
        
        self.reverse_pauli_map = [-1 for i in self.graph.data]
        for i, j in enumerate(self.pauli_map):
            # logical qubit i mapped to physical qubit j
            self.reverse_pauli_map[j] = i
        
        pass
    
    def physical_swap(self, physical_i, physical_j):
        assert self.graph.G[physical_i, physical_j] == 1
        logical_i, logical_j = self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j]
        self.pauli_map[logical_i], self.pauli_map[logical_j] = self.pauli_map[logical_j], self.pauli_map[logical_i]
        self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j] = \
            self.reverse_pauli_map[physical_j], self.reverse_pauli_map[physical_i]
    
    
        