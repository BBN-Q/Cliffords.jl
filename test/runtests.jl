using Cliffords, Base.Test
using Compat

import Cliffords: IZ, ZI, XI, IX, YI, IY, II, XX, YY, ZZ, XY, XZ, YX, ZX,
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

@test Id < X
@test X < Y
@test Y < Z

@test II < IX
@test IX < XX
@test XI < XX
@test IX < IY
@test XI < YI
@test XX < YY
@test XX < XY
@test YX > XX

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

@test Pauli([1 0; 0 1]) == Id
@test Pauli([0 1; 1 0]) == X
@test Pauli([0 -im; im 0]) == Y
@test Pauli([1 0; 0 -1]) == Z

@test complex(paulieye(2)) == eye(4)

# Cliffords

@test Clifford(1/sqrt(2) * [1 1; 1 -1]) == H
@test Clifford(eye(2)) == RI
@test Clifford(eye(4)) == kron(RI,RI)
@test Clifford([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) == CNOT
@test Clifford(diagm([1,1,1,-1])) == CZ
@test Clifford(expm(-im*pi/4*[1 0; 0 -1])) == S

function eq_upto_phase(A, B)
	idx = findfirst(A)
	rel_phase = (B[idx] == 0) ? 1.0 : (A[idx] / B[idx])
	return rel_phase * A ≈ B
end
@test eq_upto_phase(complex(RI), eye(2))
@test eq_upto_phase(complex(H), 1/sqrt(2) * [1 1; 1 -1])
@test eq_upto_phase(complex(CNOT), [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])

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

@test H * S * H == Clifford(@compat(Dict(+X=>+X,+Z=>-Y)),@compat Dict(+X=>+X,+Z=>+Y))
@test H * S * H == Clifford(expm(-im*pi/4*[0 1; 1 0]))
@test RY * RX * RY == SelfInverseClifford(@compat Dict(+X=>+X,+Z=>-Z))
@test kron(RI, H) * CNOT * kron(RI, H) == SelfInverseClifford(@compat Dict(+IZ=>+IZ,+ZI=>+ZI,+XI=>+XZ,+IX=>+ZX))
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

# Local Cliffords
@test localclifford(1) == RI

# Weight
for i=1:20
  v = rand(10)
  vX = (v .> .25) .* (v .< .5)
  vY = 2*(v .> .5) .* (v .< .75)
  vZ = 3*(v .> .75)
  r = rand()
  rS = (r > .25) + (r > .5) + (r > .75)
  @test weight(Pauli(vX + vY + vZ,rS)) == sum(v .> .25)
end

#
@test_approx_eq norm(complex(X)-[0 1; 1 0]) 0
@test_approx_eq norm(complex(Y)-[0 -1im;1im 0]) 0
@test_approx_eq norm(complex(Z)-[1 0; 0 -1]) 0
