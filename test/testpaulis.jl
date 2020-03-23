# Paulis
@test X * Id == X
@test Id * X == X
@test Y * Id == Y
@test Z * Id == Z

@test X * X == Id
@test Y * Y == Id
@test Z * Z == Id

@test X * Y == im∘Z
@test Y * Z == im∘X
@test Z * X == im∘Y
@test Y * X == -im∘Z

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
@test generators(Y) == [im∘X, Z]

@test generators(XX) == [XI, IX]
@test generators(ZZ) == [ZI, IZ]
@test generators(XY) == [im∘XI, IX, IZ]
@test generators(YX) == [im∘XI, ZI, IX]
@test generators(YY) == [-XI, ZI, IX, IZ]

XII = kron(X,Id,Id)
IXI = kron(Id,X,Id)
IIX = kron(Id,Id,X)
ZII = kron(Z,Id,Id)
IZI = kron(Id,Z,Id)
IIZ = kron(Id,Id,Z)

@test generators(kron(X,X,X)) == [XII, IXI, IIX]
@test generators(kron(Z,Z,Z)) == [ZII, IZI, IIZ]
@test generators(kron(X,Y,Z)) == [im∘XII, IXI, IZI, IIZ]

@test Pauli([1 0; 0 1]) == Id
@test Pauli([0 1; 1 0]) == X
@test Pauli([0 -im; im 0]) == Y
@test Pauli([1 0; 0 -1]) == Z

@test complex(paulieye(2)) == eye(4)

# Factors
@test factor(XII)  == [ X,Id,Id]
@test factor(XII*IZI)  == [ X,Z,Id]
@test factor(-XII)  == [-X,Id,Id]
@test factor(-XII*IZI)  == [-X,Z,Id]
@test factor(-IXI)  == [-Id,X,Id]
@test factor(-IXI*IIZ)  == [-Id,X,Z]
@test factor(im∘XII)  == [im∘X,Id,Id]
@test factor(im∘XII*IZI)  == [im∘X,Z,Id]
@test factor(im∘IXI)  == [im∘Id,X,Id]
@test factor(im∘IXI*IIZ)  == [im∘Id,X,Z]

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

# Matrix conversions
@test complex(X) == [0 1; 1 0]
@test complex(Y) == [0 -1im; 1im 0]
@test complex(Z) == [1 0; 0 -1]

# Matrix arthimetic
@test iszero(X - [0 1; 1 0])
@test iszero(Y - [0 -1im;1im 0])
@test iszero(Z - [1 0; 0 -1])

# weight returns Int
@test isa(weight(Pauli(1)), Int)

@testset "Construct from matrix" begin
    xmat = Complex{Float64}[4.702883518586134e-7 + 0.0im 0.9999991609182552 + 0.0im; 0.9999990516115099 + 0.0im 5.793263216649176e-7 + 0.0im]
    @test Pauli(xmat) == X
    ymat = -im * Complex{Float64}[7.6121754928998e-7 + 0.0im 5.577787991615943e-7 - 0.9999984033382637im; 2.066866029117562e-7 + 0.9999984033382637im 8.745277847991395e-8 + 0.0im]
    @test Pauli(ymat) == -im ∘ Y
end

