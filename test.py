from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper

from utils.mypauli import pauliString

import ipdb

mapper = JordanWignerMapper()
# driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.735", basis="sto-3g")
driver = PySCFDriver(atom="B 0 0 0; H 0 1 1; H 1 0 1; H 1 1 0", basis="sto-3g")
problem = driver.run()

def gene_uccsd_oplist(num_orbitals, num_particles):
    var_form = UCCSD(num_spatial_orbitals=num_orbitals, num_particles=num_particles, qubit_mapper=mapper)
    oplist = []
    for i in var_form._hopping_ops:
        t = []
        ps = i.to_dict()['paulis']
        for j in ps:
            t.append(pauliString(j['label'], real=j['coeff']['real'], imag=j['coeff']['imag']))
        oplist.append(t)
    return oplist


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

ipdb.set_trace()

with open('output.txt', 'w') as f:
    for pauli_list in ansatz._operators:
        f.write(pauli_list.paulis.__str__() + '\n')