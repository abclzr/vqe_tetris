# vqe_tetris

Install pyscf, qiskit 0.43.2, qiskit_nature, qiskit_optimization
```
pip install qiskit==0.43.2
pip install pytket
cd artifact_evaluation
export PYTHONPATH=../core
```
Demo (This demo is checking the correctness of tetris compiled UCCSD from LiH jordan-wigner mapper. It works under qiskit version 1.1.0 and not guaranteed to work under qiksit 0.43.2. You can skip this demo and go to "Run experiments".)
```
python3 demo_correctness.py
```
Run experiments
```
python3 run_all.py -test_scale=6
python3 calculate_duration.py
python3 show.py
python3 test_fidelity.py
```
