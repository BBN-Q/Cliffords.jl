using Cliffords, Base.Test

import Cliffords: IZ, ZI, XI, IX, YI, IY, XX, YY, ZZ, XY, XZ, YX, ZX,
	generators

# Paulis
@test X * Id == X
@test Id * X == X
@test Y * Id == Y
@test Z * Id == Z

@test X * X == Id
@test Y * Y == Id
@test Z * Z == Id

@test X * Y == im*Z
@test Y * Z == im*X
@test Z * X == im*Y
@test Y * X == -im*Z

# generators
@test generators(X) == [X]
@test generators(Z) == [Z]
@test generators(Y) == [im*X, Z]

@test generators(XX) == [XI, IX]
@test generators(ZZ) == [ZI, IZ]
@test generators(XY) == [im*XI, IX, IZ]
@test generators(YX) == [im*XI, ZI, IX]
@test generators(YY) == [-XI, ZI, IX, IZ]

XII = kron(X,Id,Id)
IXI = kron(Id,X,Id)
IIX = kron(Id,Id,X)
ZII = kron(Z,Id,Id)
IZI = kron(Id,Z,Id)
IIZ = kron(Id,Id,Z)

@test generators(kron(X,X,X)) == [XII, IXI, IIX]
@test generators(kron(Z,Z,Z)) == [ZII, IZI, IIZ]
@test generators(kron(X,Y,Z)) == [im*XII, IXI, IZI, IIZ]

@test Cliffords.Pauli([1 0; 0 1]) == Id
@test Cliffords.Pauli([0 1; 1 0]) == X
@test Cliffords.Pauli([0 -im; im 0]) == Y
@test Cliffords.Pauli([1 0; 0 -1]) == Z

@test convert(Matrix{Complex{Int}},Cliffords.Pauli(eye(4))) == eye(4)

# Cliffords

@test Clifford(1/sqrt(2) * [1 1; 1 -1]) == H
@test Clifford(eye(2)) == RI
@test Clifford(eye(4)) == kron(RI,RI)
@test Clifford([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) == CNOT
@test Clifford(diagm([1,1,1,-1])) == CZ
@test Clifford(expm(-im*pi/4*[1 0; 0 -1])) == S

# Cliffords * Paulis
@test H * Id == Id
@test H * X == Z
@test H * Z == X
@test H * Y == -Y

II = kron(Id,Id)

@test kron(RI, H) * II == II
@test kron(RI, H) * IX == IZ
@test kron(RI, H) * IZ == IX
@test kron(RI, H) * IY == -IY

@test CNOT * II == II
@test CNOT * IX == IX
@test CNOT * XI == XX
@test CNOT * YY == -XZ

# Cliffords \ Paulis
@test H \ Id == Id
@test H \ X == Z
@test H \ Z == X
@test H \ Y == -Y

II = kron(Id,Id)

@test kron(RI, H) \ II == II
@test kron(RI, H) \ IX == IZ
@test kron(RI, H) \ IZ == IX
@test kron(RI, H) \ IY == -IY

@test CNOT \ II == II
@test CNOT \ IX == IX
@test CNOT \ XI == XX
@test CNOT \ YY == -XZ

# Cliffords * Cliffords

@test H * S * H == Clifford([+X=>+X,+Z=>-Y],[+X=>+X,+Z=>+Y])
@test H * S * H == Clifford(expm(-im*pi/4*[0 1; 1 0]))
@test RY * RX * RY == SelfInverseClifford([+X=>+X,+Z=>-Z])
@test kron(RI, H) * CNOT * kron(RI, H) == SelfInverseClifford([+IZ=>+IZ,+ZI=>+ZI,+XI=>+XZ,+IX=>+ZX])
@test kron(RI, H) * CNOT * kron(RI, H) == CZ
@test S * S == RZ

# 3 CNOTs is a SWAP
CNOT21 = expand(CNOT, [2,1], 2)
@test CNOT * CNOT21 * CNOT == SWAP

# iSWAP
iSWAP = CZ*kron(S,S)*SWAP
@test iSWAP == Clifford([1 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1])

# Cliffords \ Cliffords
@test H \ H == H * H
@test CNOT \ CNOT == CNOT * CNOT
@test RY * RX * RY == RY \ (RX * RY)
@test S \ S == RI
@test S \ RZ == S
