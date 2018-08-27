# Cliffords

@test Clifford(1/sqrt(2) * [1 1; 1 -1]) == H
@test Clifford(eye(2)) == RI
@test Clifford(eye(4)) == kron(RI,RI)
@test Clifford([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) == CNOT
@test Clifford(diagm(0 => [1,1,1,-1])) == CZ
@test Clifford(exp(-im*pi/4*[1 0; 0 -1])) == S

function eq_upto_phase(A, B)
    idx = findfirst(x -> x != 0, A)
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

@test kron(RI, H) \ II == II
@test kron(RI, H) \ IX == IZ
@test kron(RI, H) \ IZ == IX
@test kron(RI, H) \ IY == -IY

@test CNOT \ II == II
@test CNOT \ IX == IX
@test CNOT \ XI == XX
@test CNOT \ YY == -XZ

# Cliffords * Cliffords

@test H * S * H == Clifford(Dict(+X=>+X,+Z=>-Y),Dict(+X=>+X,+Z=>+Y))
@test H * S * H == Clifford(exp(-im*pi/4*[0 1; 1 0]))
@test RY * RX * RY == SelfInverseClifford(Dict(+X=>+X,+Z=>-Z))
@test kron(RI, H) * CNOT * kron(RI, H) == SelfInverseClifford(Dict(+IZ=>+IZ,+ZI=>+ZI,+XI=>+XZ,+IX=>+ZX))
@test kron(RI, H) * CNOT * kron(RI, H) == CZ
@test S * S == RZ

# 3 CNOTs is a SWAP
CNOT21 = Cliffords.expand(CNOT, [2,1], 2)
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

# Cliffords -> Pauli
@test Clifford(Id) == RI
@test Clifford(X)  == RX
@test Clifford(Y)  == RY
@test Clifford(Z)  == RZ
@test Clifford(XX) == kron(RX,RX)
@test Clifford(XI) == kron(RX,RI)

# Clifford and Matrix interaction
@test typeof(promote(H, eye(2))) == Tuple{Array{Float64,2}, Array{Float64,2}}
@test typeof(promote(H, eye(ComplexF64, 2))) == Tuple{Array{ComplexF64,2}, Array{ComplexF64,2}}
