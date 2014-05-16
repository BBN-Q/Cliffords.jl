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

# Cliffords * Paulis
# @test H * Id == Id # failing
@test H * X == Z
@test H * Z == X
@test H * Y == -Y

@test CNOT * IX == IX
@test CNOT * XI == XX
@test CNOT * YY == -XZ

# Cliffords * Cliffords

@test H * S * H == Clifford([+X=>+X,+Z=>-Y],[+X=>+X,+Z=>+Y])
@test RY * RX * RY == SelfInverseClifford([+X=>+X,+Z=>-Z])
@test kron(RI, H) * CNOT * kron(RI, H) == SelfInverseClifford([+IZ=>+IZ,+ZI=>+ZI,+XI=>+XZ,+IX=>+ZX])

# 3 CNOTs is a SWAP
CNOT21 = expand(CNOT, [2,1], 2)
@test CNOT * CNOT21 * CNOT == SelfInverseClifford([+IZ=>+ZI,+ZI=>+IZ,+XI=>+IX,+IX=>+XI])

# inverse circuits...

