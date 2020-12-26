# Cliffords

![Build Status](https://github.com/BBN-Q/Cliffords.jl/workflows/CI/badge.svg)

[![codecov](https://codecov.io/gh/BBN-Q/Cliffords.jl/branch/master/graph/badge.svg?token=Zz2pT31upJ)](https://codecov.io/gh/BBN-Q/Cliffords.jl)

This library allows for efficient calculation of Clifford circuits by tracking the evolution of X and Z generators (the so-called **tableau** representation). No special effort has been made to strictly minimize the number of bits needed to store each Clifford. Rather, the goal was clarity. One unique feature compared to other such utilities is that we also efficiently track the inverse operations. This is useful to, e.g., compute 'undo' gates in randomized benchmarking sequences.

## Usage

```
using Cliffords

# single-qubit Pauli operators
X * Y => iZ

Pauli([0 1; 1 0]) => +X

# multi-qubit Pauli operators
kron(X,X) * kron(Z,Z) => -YY

# Cliffords
Clifford([1 0 0 0;
          0 1 0 0;
          0 0 0 1;
          0 0 1 0]) == CNOT => true

# Clifford * Pauli
H * X => Z
H * Z => X

# Clifford * Clifford
CNOT21 = expand(CNOT, [2,1], 2)
CNOT * CNOT21 * CNOT == SWAP => true
```
