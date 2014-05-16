# Copyright 2014: Raytheon BBN Technologies
# Original authors: Blake Johnson and Marcus da Silva

module Cliffords

import Base: kron, length
export Clifford, SelfInverseClifford, expand,
	RI, RX, RY, RZ, H, S, CNOT

include("Paulis.jl")

immutable Clifford
	T::Dict{Pauli, Pauli}
	Tinv::Dict{Pauli, Pauli}
end

SelfInverseClifford(T) = Clifford(T, T)
length(c::Clifford) = length(first(keys(c.T)))

==(a::Clifford, b::Clifford) = (a.T == b.T) && (a.Tinv == b.Tinv)

# TODO
# function Clifford(label::String, U::Matrix)

# end

const RI = SelfInverseClifford([Z => Z, X => X])
const H = SelfInverseClifford([Z => X, X => Z])
const S = Clifford([Z => Z, X => Y], [Z => Z, X => -Y])
const CNOT = SelfInverseClifford([ZI => ZI, XI => XX, IZ => ZZ, IX => IX])
const RX = SelfInverseClifford([Z => -Z, X => X])
const RY = SelfInverseClifford([Z => -Z, X => -X])
const RZ = SelfInverseClifford([Z => Z, X => -X])

function *(a::Clifford, b::Clifford)
	T = (Pauli=>Pauli)[]
	for p = keys(b.T)
		T[p] = a * (b * p)
	end
	Tinv = (Pauli=>Pauli)[]
	for p = keys(a.Tinv)
		Tinv[p] = b \ (a \ p)
	end
	Clifford(T, Tinv)
end

function *(c::Clifford, p::Pauli)
    if isid(p)
        return p
    else
	# rewrite p in terms of generators (X and Z)
	G = generators(p)
	r = paulieye(length(p))
	for g in G
		r *= phase(g) * c.T[abs(g)]
	end
	return r
    end
end

function \(c::Clifford, p::Pauli)
    if isid(p)
        return p
    else
	G = generators(p)
	r = paulieye(length(p))
	for g in G
		r *= phase(g) * c.Tinv[abs(g)]
	end
	return r
    end
end

function expand(c::Clifford, subIndices, n)
	T = (Pauli=>Pauli)[]
	for (k,v) in c.T
		T[expand(k, subIndices, n)] = expand(v, subIndices, n)
	end
	Tinv = (Pauli=>Pauli)[]
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
	expand(a, [1:length(a)], n) * expand(b, [length(a)+1:n], n)
end

end