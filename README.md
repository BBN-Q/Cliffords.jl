# Cliffords

[![Build Status](https://travis-ci.org/BBN-Q/Cliffords.jl.svg?branch=master)](https://travis-ci.org/BBN-Q/Cliffords.jl)

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

## Installation 
For Julia v0.6 run:

```julia> Pkg.add("Cliffords")```

For v0.7 and above, run the following in the package manager to grab the master branch:

```(v1.1) pkg> add Cliffords#master```
