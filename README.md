# vqe_tetris

Install pyscf, qiskit 0.43.2, qiskit_nature, qiskit_optimization
```
pip install qiskit==0.43.2
pip install pytket
```
Run experiments
```
export PYTHONPATH=core
cd artifact_evaluation
python3 run_all.py -test_scale=1
```
