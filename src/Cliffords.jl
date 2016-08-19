# Copyright 2014: Raytheon BBN Technologies
# Original authors: Blake Johnson and Marcus da Silva

VERSION >= v"0.5-" && __precompile__()

module Cliffords

import Base: convert, show, kron, abs, length, hash, isequal, vec, promote_rule,
    zero, inv, expand, ==, *, +, -, \, isless, ctranspose
import Iterators: product

export Clifford, SelfInverseClifford, expand,
       RI, RX, RY, RZ, H, S, CNOT, CZ, SWAP, cliffordeye

using Compat, FixedSizeArrays

include("Paulis.jl")

type Clifford{N}
    T::Dict{Pauli{N}, Pauli{N}}
    Tinv::Dict{Pauli{N}, Pauli{N}}
end

SelfInverseClifford(T) = Clifford(T, T)
length{N}(c::Clifford{N}) = N

=={N}(a::Clifford{N}, b::Clifford{N}) = (a.T == b.T)
isequal{N}(a::Clifford{N}, b::Clifford{N}) = (a == b) # for backward compatibility with Julia 0.2
hash{N}(c::Clifford{N}, h::UInt) = hash(c.T, h)

function convert{N}(::Type{Clifford{N}},U::Matrix)
    # N = round(Int,log(2,size(U,1)))
    T = Dict{Pauli{N},Pauli{N}}()
    Tinv = Dict{Pauli{N},Pauli{N}}()
    t = typeof(complex(U))
    n = round(Int, log(2,size(U,1)))
    ri = cliffordeye(n)
    for p in keys(ri.T)
        T[p] = U * p * U'
        Tinv[p] = U' * p * U
    end
    Clifford(T, Tinv)
end

function Clifford(U::Matrix)
    N = round(Int64,log(2,size(U,1)))
    convert(Clifford{N},U)
end

function convert{T}(::Type{Matrix{T}}, c::Clifford)
    d = 2^length(c)
    l = liou(c)

    # liou2choi
    rl = reshape(l, (d, d, d, d) )
    rl = permutedims(rl, [1,3,2,4])
    choi = reshape(rl, size(l))

    # unitary is the reshaped eigenvector corresponding to the non-zero eigenvalue
    v = eigfact(Hermitian(choi), d^2:d^2)[:vectors]
    m = reshape(v, d, d)
    return convert(Matrix{T}, m * sqrt(d))
end
complex(c::Clifford) = convert(Matrix{Complex128},c)

function liou(c::Clifford)
    d = 4^length(c)
    m = zeros(Complex128,d,d)
    for p in allpaulis(length(c))
        m += vec(c*p)*vec(p)'/sqrt(d)
    end
    m
end

promote_rule{T}(::Type{Clifford}, ::Type{Matrix{T}}) = Matrix{T}

const RI = SelfInverseClifford(@compat Dict(Z => Z, X => X))
const H = SelfInverseClifford(@compat Dict(Z => X, X => Z))
const S = Clifford(@compat(Dict(Z => Z, X => Y)), @compat Dict(Z => Z, X => -Y))
const CNOT = SelfInverseClifford(@compat Dict(ZI => ZI, XI => XX, IZ => ZZ, IX => IX))
const CZ   = SelfInverseClifford(@compat Dict(ZI => ZI, XI => XZ, IZ => IZ, IX => ZX))
const SWAP = SelfInverseClifford(@compat Dict(ZI => IZ, XI => IX, IZ => ZI, IX => XI))
const RX = SelfInverseClifford(@compat Dict(Z => -Z, X => X))
const RY = SelfInverseClifford(@compat Dict(Z => -Z, X => -X))
const RZ = SelfInverseClifford(@compat Dict(Z => Z, X => -X))

function *{N}(a::Clifford{N}, b::Clifford{N})
    T = Dict{Pauli{N},Pauli{N}}()
    for p = keys(b.T)
        T[p] = a * (b * p)
    end
    Tinv = Dict{Pauli{N},Pauli{N}}()
    for p = keys(a.Tinv)
        Tinv[p] = b \ (a \ p)
    end
    Clifford(T, Tinv)
end

function \{N}(a::Clifford{N}, b::Clifford{N})
    T = Dict{Pauli{N},Pauli{N}}()
    for p = keys(b.T)
        T[p] = a \ (b * p)
    end
    Tinv = Dict{Pauli{N},Pauli{N}}()
    for p = keys(a.Tinv)
        Tinv[p] = b \ (a * p)
    end
    Clifford(T, Tinv)
end

function *(c::Clifford, p::Pauli)
    if isid(p)
        return p
    end
    # rewrite p in terms of generators (X and Z)
    G = generators(p)
    r = paulieye(length(p))
    for g in G
        r *= phase(g) * c.T[abs(g)]
    end
    return r
end

function \(c::Clifford, p::Pauli)
    if isid(p)
        return p
    end
    G = generators(p)
    r = paulieye(length(p))
    for g in G
        r *= phase(g) * c.Tinv[abs(g)]
    end
    return r
end

*(c::Clifford, m::Matrix) = *(promote(c,m)...)
*(m::Matrix, c::Clifford) = *(promote(m,c)...)

inv(c::Clifford) = Clifford(c.Tinv, c.T)
ctranspose(c::Clifford) = inv(c)

function expand{N}(c::Clifford{N}, subIndices, n)
    T = Dict{Pauli{n},Pauli{n}}()
    for (k,v) in c.T
        T[expand(k, subIndices, n)] = expand(v, subIndices, n)
    end
    Tinv = Dict{Pauli{n}, Pauli{n}}()
    for (k,v) in c.Tinv
        Tinv[expand(k, subIndices, n)] = expand(v, subIndices, n)
    end
    # add trivial mapping for missing subspaces
    for dim in setdiff(1:n, subIndices)
        T[expand(X, dim, n)] = expand(X, dim, n)
        T[expand(Z, dim, n)] = expand(Z, dim, n)
        Tinv[expand(X, dim, n)] = expand(X, dim, n)
        Tinv[expand(Z, dim, n)] = expand(Z, dim, n)
    end
    Clifford(T, Tinv)
end

function kron(a::Clifford, b::Clifford)
    n = length(a) + length(b)
    expand(a, [1:length(a);], n) * expand(b, [length(a)+1:n;], n)
end

zero(::Type{Clifford}) = RI
cliffordeye(n) = expand(RI, [1], n)

include("C1.jl")

end
