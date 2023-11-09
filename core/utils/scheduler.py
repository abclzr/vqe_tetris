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
        self.is_ancilla = [True for i in self.graph.data]
        for p_q in pauli_map:
            self.is_ancilla[p_q] = False
        
        for i, j in enumerate(self.pauli_map):
            # logical qubit i mapped to physical qubit j
            self.reverse_pauli_map[j] = i
        
        logical_qubit_id = len(self.pauli_map)
        for i in range(len(self.reverse_pauli_map)):
            if self.reverse_pauli_map[i] == -1:
                self.reverse_pauli_map[i] = logical_qubit_id
                logical_qubit_id = logical_qubit_id + 1
        assert logical_qubit_id == len(self.reverse_pauli_map)
        
        self.pauli_map = [0 for i in self.reverse_pauli_map]
        for i, j in enumerate(self.reverse_pauli_map):
            self.pauli_map[j] = i
        
        # apply floyd_warshall algorithm to calculate the distance
        self.distance = floyd_warshall(self.graph.G)
        self.tree = Tree([], 0)
        
        # enable_cancel == True means you can only record the instructions not do any qubit routing
        self.enable_cancel = True
        
        self.instruction_list = []
        self.record = []
        self.not_compiled_pointer = 0
        self.total_logical_instruction = 0
        self.canceled_logical_instruction = 0
        self.total_swap_cnt = 0
        self.total_cx_cnt = 0
        self.total_bridge_cnt = 0
    
    def notify_ancilla(self, logical_i):
        physical_i = self.pauli_map[logical_i]
        self.is_ancilla[physical_i] = True
    
    def physical_swap(self, physical_i, physical_j):
        self.total_swap_cnt = self.total_swap_cnt + 1
        self.qc.swap(physical_i, physical_j)
        # assert self.graph.G[physical_i, physical_j] == 1
        logical_i, logical_j = self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j]
        self.pauli_map[logical_i] = physical_j
        self.pauli_map[logical_j] = physical_i
        self.reverse_pauli_map[physical_i], self.reverse_pauli_map[physical_j] = \
            self.reverse_pauli_map[physical_j], self.reverse_pauli_map[physical_i]
        self.is_ancilla[physical_i], self.is_ancilla[physical_j] = \
            self.is_ancilla[physical_j], self.is_ancilla[physical_i]
        return
    
    def find_centor(self, nodes):
        centor = -1
        min_total_distance = 0x7777777
        for c in range(len(self.reverse_pauli_map)):
            total_distance = 0
            for n in nodes:
                total_distance = total_distance + self.distance[self.pauli_map[n]][c]
            if total_distance < min_total_distance:
                min_total_distance = total_distance
                centor = c
        return centor
    
    def gather_root_tree(self, nodes, centor):
        close_to_far = sorted(nodes, key=lambda n: self.distance[self.pauli_map[n]][centor])
        connected_component = []
        edges = []
        for n in close_to_far:
            min_dis = 0x7777777
            closest_dest = -1
            for dest in connected_component + [centor]:
                if self.distance[self.pauli_map[n]][dest] < min_dis:
                    min_dis = self.distance[self.pauli_map[n]][dest]
                    closest_dest = dest
            path = self.shortest_path(self.pauli_map[n], closest_dest)
            for edge in path[:-1]:
                self.physical_swap(edge[0], edge[1])
            
            if len(path) == 0:
                # the node sits on the centor
                pass
            elif not self.reverse_pauli_map[path[-1][1]] in nodes:
                # the node doesn't sit on the centor and centor is empty
                if path[-1][1] != centor:
                    pdb.set_trace()
                assert path[-1][1] == centor
                self.physical_swap(path[-1][0], path[-1][1])
            else:
                edges.append((n, self.reverse_pauli_map[path[-1][1]]))
            # now n is attached to the connected component
            connected_component.append(self.pauli_map[n])
        return connected_component, edges
    
    def gather_leaf_tree(self, leaf_nodes, connected_component, n_paulistring, use_bridge, swap_coefficient=3):
        leaf_nodes = sorted(leaf_nodes, key=lambda leaf: min([self.distance[self.pauli_map[leaf]][c] for c in connected_component]))
        connected_leaf_physical = []
        edges = []
        for leaf in leaf_nodes:
            min_dis = 0x7777777
            closest_dest = -1
            for dest in connected_component:
                if (self.distance[self.pauli_map[leaf]][dest] - 1) * swap_coefficient + 2 * n_paulistring < min_dis:
                    min_dis = (self.distance[self.pauli_map[leaf]][dest] - 1) * swap_coefficient + 2 * n_paulistring
                    closest_dest = dest
            for dest in connected_leaf_physical:
                if (self.distance[self.pauli_map[leaf]][dest] - 1) * swap_coefficient + 2 < min_dis:
                    min_dis = (self.distance[self.pauli_map[leaf]][dest] - 1) * swap_coefficient + 2
                    closest_dest = dest
            
            path = self.shortest_path(self.pauli_map[leaf], closest_dest)
            assert len(path) > 0
            tmp = 0
            while tmp < len(path) and not path[tmp][1] in connected_component + connected_leaf_physical:
                tmp = tmp + 1
            closest_dest = path[tmp][1]
            path = path[:tmp]
            
            bridge_edges = []
            if use_bridge:
                while len(path) > 0 and self.is_ancilla[path[-1][1]]:
                    bridge_edges.append(path.pop())
                if bridge_edges != []:
                    # print(bridge_edges)
                    self.total_bridge_cnt = self.total_bridge_cnt + len(bridge_edges)
            # for edge in path[:-1]:
            #     if not edge[1] in connected_component + connected_leaf_physical:
            #         self.physical_swap(edge[0], edge[1])
            #     else:
            #         closest_dest = edge[1]
            #         break
            for edge in path:
                self.physical_swap(edge[0], edge[1])
            for edge in reversed(bridge_edges):
                edges.append((self.reverse_pauli_map[edge[0]], self.reverse_pauli_map[edge[1]]))
            
            if bridge_edges == []:
                edges.append((leaf, self.reverse_pauli_map[closest_dest]))
            else:
                # edges.append((leaf, self.reverse_pauli_map[bridge_edges[-1][0]]))
                edges.append((self.reverse_pauli_map[bridge_edges[0][1]], self.reverse_pauli_map[closest_dest]))
            
            connected_leaf_physical.append(self.pauli_map[leaf])
            for edge in bridge_edges:
                connected_leaf_physical.append(edge[0])
                connected_leaf_physical.append(edge[1])
        return edges

    
    def MST_init(self, n_nodes):
        self.union_find = UnionFind(n_nodes)
    
    def MST(self, nodes, edges):
        edges = sorted(edges, key=lambda x: self.distance[self.pauli_map[x[0]]][self.pauli_map[x[1]]])
        mst_edges = kruskal_mst(edges, self.union_find, self.distance)
        return mst_edges
    
    def accept_MST(self, edges):
        for u, v in edges:
            self.union_find.union(u, v)
    
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
                price = 0
            elif instruction.startswith('Logical_left_Y'):
                self.qc.u(np.pi/2, -np.pi/2, np.pi/2, self.pauli_map[data])
                price = 0
            elif instruction.startswith('Logical_CNOT'):
                u, v = data
                p_u, p_v = self.pauli_map[u], self.pauli_map[v]
                path = self.shortest_path(p_u, p_v)
                for u, v in path[:-1]:
                    self.physical_swap(u, v)
                self.qc.cx(path[-1][0], path[-1][1])
                self.total_cx_cnt += 1
                price = 3 * len(path) - 3 + 1
            elif instruction.startswith('Logical_SWAP'):
                u, v = data
                p_u, p_v = self.pauli_map[u], self.pauli_map[v]
                path = self.shortest_path(p_u, p_v)
                for u, v in path:
                    self.physical_swap(u, v)
                price = 3 * len(path)
            elif instruction.startswith('Logical_RZ'):
                self.qc.rz(1, self.pauli_map[data])
                price = 0
            elif instruction.startswith('Logical_right_X'):
                self.qc.u(np.pi/2, 0, np.pi, self.pauli_map[data])
                price = 0
            elif instruction.startswith('Logical_right_Y'):
                self.qc.u(-np.pi/2, -np.pi/2, np.pi/2, self.pauli_map[data])
                price = 0
            else:
                raise Exception('Illegal instruction: ' + instruction)
            self.record.append((instruction, data, price))
        
        self.not_compiled_pointer = len(self.instruction_list)
    
    def add_instruction(self, instruction, data):
        self.total_logical_instruction = self.total_logical_instruction + 1
        if self.enable_cancel == False:
            self.instruction_list.append((instruction, data))
            return

        # check if there is any cancellation
        tmp = len(self.instruction_list) - 1
        if instruction.startswith('Logical_CNOT') or instruction.startswith('Logical_SWAP'):
            set_a = set(data)
        else:
            set_a = set([data])
        while tmp >= self.not_compiled_pointer:
            if self.instruction_list[tmp][0].startswith('Logical_CNOT') or self.instruction_list[tmp][0].startswith('Logical_SWAP'):
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
            elif instruction.startswith('Logical_SWAP'):
                if self.instruction_list[tmp][0] == 'Logical_SWAP':
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
