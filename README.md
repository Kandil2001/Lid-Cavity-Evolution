# Lid Cavity Evolution

**Lid Cavity Evolution** is an ambitious, growing repository dedicated to exploring and comparing algorithms for solving the classic lid-driven cavity problem in computational fluid dynamics (CFD). This project showcases the progression of solver implementations, starting from straightforward and basic algorithms, moving through vectorized and parallel approaches, and striving toward ever more efficient and elegant solutions.

## Goals

- To serve as a benchmark suite for lid-driven cavity solvers in MATLAB and Python.
- To provide educational insight into how algorithms evolve for speed, simplicity, and scalability.
- To enable clear comparisons of solving time, accuracy, and code complexity across approaches.

## Features

- **MATLAB Implementations:** Includes both loop-based and vectorized SIMPLE solvers.
- **Python Implementations:** Serial and parallel versions (coming soon).
- **Benchmarking:** Tools and scripts for measuring and visualizing solving time and performance.
- **Documentation:** Progressive explanations for each algorithm and its optimizations.

## Structure

```
matlab/
    SimpleLidcavity.m
    SimpleLidCavityVector.m
python/
    lid_cavity_serial.py         # To be added
    lid_cavity_parallel.py       # To be added
README.md
```

## Roadmap

- [x] MATLAB base solvers (loop-based & vectorized)
- [ ] Python serial implementation
- [ ] Python parallel implementation
- [ ] Benchmarking and performance comparison
- [ ] Advanced CFD extensions

---

**Contributions, ideas, and feedback are welcome as Lid Cavity Evolution grows!**
