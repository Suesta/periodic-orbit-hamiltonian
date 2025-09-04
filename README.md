# Periodic Orbit Computation in Hamiltonian Systems

Development of numerical methods in C for the computation of periodic orbits in Hamiltonian systems.

## 🚀 Description

This project implements advanced numerical routines to detect and refine periodic orbits in Hamiltonian dynamical systems. It was developed by Víctor Suesta at the age of 23 as part of an advanced scientific computing study. The core functionalities include:

- **Runge-Kutta-Fehlberg 7-8 integration (RKF78)** for high-precision trajectory propagation.
- **Variational equation flow integration**, used for stability analysis and Newton refinement.
- **QR decomposition routine**, essential for handling the monodromy matrix and the stability multipliers.
- **Newton’s method for multiple variables**, applied to locate and refine periodic orbits using shooting techniques.

All implementations are written in pure C, without external libraries, and optimized for clarity and scientific reproducibility.

## 📁 Project Structure

- `main.c`: Entry point and main simulation logic.
- `Pr*.c`: Different examples and parameter configurations.
- `FUNCTIONS.c`: Core numerical integration and Newton method functions.
- `FuncionsQR.h`: Custom header for the QR decomposition.
- `Instruccions per executar.txt`: Execution instructions (in Catalan).

## 🧪 How to Run

1. Compile the desired source file using `gcc`. Example:
   ```bash
   gcc -o orbit Pr4Ex2.c FUNCTIONS.c -lm


## 🔗 Author
Víctor Suesta  
[LinkedIn](https://www.linkedin.com/in/víctor-suesta-arribas-7b1250322/)
