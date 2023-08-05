from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper, BravyiKitaevSuperFastMapper, BravyiKitaevMapper
from benchmark.mypauli import pauliString

import pickle
import ast
import ipdb

mapper = JordanWignerMapper()
# mapper = ParityMapper()
# mapper = BravyiKitaevMapper()
# mapper = BravyiKitaevSuperFastMapper()

# driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.735", basis="sto-3g")
# driver = PySCFDriver(atom="B 0 0 0; H 0 1 1; H 1 0 1; H 1 1 0", basis="sto-3g")
driver = PySCFDriver(atom="Li .0 .0 .0; Cl .0 .0 -1.5", basis='sto3g')
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

# ipdb.set_trace()

pauli_blocks = []
for pauli_list in ansatz._operators:
    block = ast.literal_eval(pauli_list._primitive._pauli_list.__str__())
    block = [pauliString(ps) for ps in block]
    pauli_blocks.append(block)

with open('data/jordan_wigner/LiCl_UCCSD.pickle', 'wb') as f:
    pickle.dump(pauli_blocks, f)
