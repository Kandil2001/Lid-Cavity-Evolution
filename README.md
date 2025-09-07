# Lid Cavity Evolution

A benchmark suite for comparing algorithms that solve the classic **lid-driven cavity** problem in Computational Fluid Dynamics (CFD). This project showcases the evolution of solver implementations from basic educational code to optimized, production-ready versions.

## 🚀 Current Features

- **MATLAB Solvers**:
  - `SimpleLidCavity.m` - A basic, loop-based implementation of the SIMPLE algorithm. Ideal for understanding the fundamentals.
  - `SimpleLidCavityVector.m` - A significantly faster, vectorized version of the same solver. Demonstrates the power of MATLAB's array operations.

## 🗺️ Roadmap & Coming Soon

- [ ] Python serial implementation (using NumPy)
- [ ] Python parallel implementation (using Dask or Multiprocessing)
- [ ] Automated benchmarking scripts to compare solver performance
- [ ] Validation against canonical data from Ghia et al. (1982)

## ⏱️ The Goal: Performance Comparison

The ultimate aim of this repository is to provide a clear, fair comparison of solving time and accuracy across different implementations and languages.

| Solver | Language | Paradigm | Avg. Time (s) | Grid Size | Re |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `SimpleLidCavity` | MATLAB | Serial (Loops) | TBD | 128x128 | 1000 |
| `SimpleLidCavityVector` | MATLAB | Serial (Vectorized) | TBD | 128x128 | 1000 |
| *More coming...* | | | | | |

## 🤝 Contributing

Ideas, contributions, and feedback are welcome! Feel free to open an issue to discuss new solver implementations, optimizations, or ideas for comparison.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
