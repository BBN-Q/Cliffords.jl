import Base: convert, show, kron, abs, length, hash
export Id, X, Y, Z

immutable Pauli
    v::Vector{Uint8} # 0 = I, 1 = X, 2 = Z, 3 = Y
    s::Uint8 # 0 = +1, 1 = +i, 2 = -1, 3 = -i (im^s)
    function Pauli(v::Vector{Uint8}, s::Uint8)
        if any(v .> 3)
            error("Invalid Pauli operator.")
        else
            new(v,mod(s,4))
        end
    end
end

Pauli(v::Vector, s = 0) = Pauli(uint8(v), uint8(s))
Pauli(v::Integer, s = 0) = Pauli([v], s)

==(a::Pauli, b::Pauli) = (a.v == b.v && a.s == b.s)
hash(a::Pauli) = hash(convert(String, a))

function convert(::Type{String}, p::Pauli)
    phases = ["+","i","-","-i"]
    paulis = "IXZY"
    phases[p.s+1] * join([paulis[i+1] for i in p.v])
end

function show(io::IO, p::Pauli)
    print(io,convert(String,p))
end

function levicivita(a::Uint8, b::Uint8)
    # an unusual Levi-Civita pseudo-tensor for the (1,2,3) = (X,Z,Y) convention
    if (a,b) == (1,3) || (a,b) == (3,2) || (a,b) == (2,1)
        1
    elseif (a,b) == (1,2) || (a,b) == (2,3) || (a,b) == (3,1)
        3
    else
        0
    end
end

function levicivita(a::Vector{Uint8}, b::Vector{Uint8})
    lc = sum(map(levicivita, a, b))
    mod(lc, 4)
end

*(a::Pauli, b::Pauli) = Pauli(a.v $ b.v, a.s + b.s + levicivita(a.v, b.v))

function *(n::Number, p::Pauli)
    ns = findfirst(n .== [1, im, -1, -im]) - 1
    @assert(ns >= 0, "Multiplication only supported for +/- 1, +/- im")
    Pauli(p.v, p.s + ns)
end
*(p::Pauli, n::Number) = n * p
+(p::Pauli) = p
-(p::Pauli) = Pauli(p.v, p.s + 2)
abs(p::Pauli) = Pauli(p.v, 0)
phase(p::Pauli) = im^p.s
length(p::Pauli) = length(p.v)

kron(a::Pauli, b::Pauli) = Pauli([a.v, b.v], a.s + b.s)

function expand(a::Pauli, subIndices, n)
    v = zeros(n)
    for (ct, i) in enumerate(subIndices)
        v[i] = a.v[ct]
    end
    Pauli(v, a.s)
end

function generators(a::Pauli)
    G = Pauli[]
    s = phase(a)
    for (idx, p) in enumerate(a.v)
        if p == 0 # I, skip it
            continue
        elseif p == 1 || p == 2 # X or Z
            push!(G, expand(Pauli(p), [idx], length(a.v)))
        else # Y
            s *= im
            push!(G, expand(X, [idx], length(a.v)))
            push!(G, expand(Z, [idx], length(a.v)))
        end
    end
    G[1] *= s # give phase to first generator
    return G
end

function paulieye(n)
    expand(Id, [1], n)
end

# 1-qubit Paulis
const Id = Pauli(0)
const X = Pauli(1)
const Y = Pauli(3)
const Z = Pauli(2)

# 2-qubit Paulis
labelOpPairs = [("I", Id), ("X", X), ("Y", Y), ("Z", Z)]
for (a, ao) in labelOpPairs, (b, bo) in labelOpPairs
    @eval const $(symbol(a*b)) = kron($ao, $bo)
end
