using Cliffords, Base.Test

import Cliffords: IZ, ZI, XI, IX, YI, IY, XX, YY, ZZ, XZ, ZX

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

# Cliffords * Paulis

# Cliffords * Cliffords

@test H * S * H == Clifford([+X=>+X,+Z=>-Y],[+X=>+X,+Z=>+Y])
@test RY * RX * RY == SelfInverseClifford([+X=>+X,+Z=>-Z])
@test kron(RI, H) * CNOT * kron(RI, H) == SelfInverseClifford([+IZ=>+IZ,+ZI=>+ZI,+XI=>+XZ,+IX=>+ZX])
