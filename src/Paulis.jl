import Base.complex

export Pauli, Id, X, Y, Z, allpaulis, paulieye, weight, complex

# Paulis's are represented by a vector of numbers (0-3) corresponding to
# single-qubit Paulis, along with a phase parameter.
immutable Pauli
    v::Vector{UInt8} # 0 = I, 1 = X, 2 = Z, 3 = Y
    s::UInt8 # 0 = +1, 1 = +i, 2 = -1, 3 = -i (im^s)
    function Pauli(v::Vector{UInt8}, s::UInt8)
        new(v,mod(s,4))
    end
end

Pauli(v::Vector, s = 0) = Pauli(convert(Vector{UInt8}, v), convert(UInt8, s))
Pauli(v::Integer, s = 0) = Pauli([v], s)
Pauli(m::Matrix) = convert(Pauli, m)

weight(p::Pauli) = sum( p.v .> 0 )

show(io::IO, p::Pauli) = print(io,convert(String,p))

==(a::Pauli, b::Pauli) = (a.v == b.v && a.s == b.s)
isequal(a::Pauli, b::Pauli) = (a == b)
hash(a::Pauli, h::UInt) = hash(a.v, hash(a.s, h))
isid(a::Pauli) = isempty(find(a.v))

function isless(a::Pauli, b::Pauli)
    # canonical total order defined by weight and then "lexicographic":
    # Id < X < Y < Z
    if weight(a) != weight(b)
        return weight(a) < weight(b)
    else
        return lex_tuple(a) < lex_tuple(b)
    end
end

function lex_tuple(p::Pauli)
    const pauli_lex = (1, 2, 4, 3)
    tuple((pauli_lex[x+1] for x in p.v)...)
end

function convert(::Type{String}, p::Pauli)
    phases = ["+","i","-","-i"]
    paulis = "IXZY"
    phases[p.s+1] * join([paulis[i+1] for i in p.v])
end

function convert{T}(::Type{Matrix{Complex{T}}}, p::Pauli)
    const mats = @compat Dict(
        0x00 => eye(Complex{T},2),
        0x01 => Complex{T}[0 1; 1 0],
        0x02 => Complex{T}[1 0; 0 -1],
        0x03 => Complex{T}[0 -im; im 0])

    return phase(p)*reduce(kron,[mats[x] for x in p.v])
end

complex(p::Pauli) = convert(Matrix{Complex128},p)

function convert(::Type{Pauli}, m::Matrix)
    d = size(m,1)
    n = log(2,size(m,1))
    for p in allpaulis(n)
        overlap = trace(m*convert(typeof(complex(m)),p)) / d
        if isapprox(abs(overlap),1,atol=d*eps(Float64))
            return (round(real(overlap))+im*round(imag(overlap)))*p
        elseif !isapprox(abs(overlap),0,atol=d*eps(Float64))
            println(m, overlap,isapprox(abs(overlap),0))
            error("Trying to convert non-Pauli matrix to a Pauli object")
        end
    end
end

promote_rule{T<:Real}(::Type{Pauli}, ::Type{Matrix{T}}) = Matrix{Complex{T}}
promote_rule{T<:Complex}(::Type{Pauli}, ::Type{Matrix{T}}) = Matrix{T}

function levicivita(a::UInt8, b::UInt8)
    # an unusual Levi-Civita pseudo-tensor for the (1,2,3) = (X,Z,Y) convention
    if (a,b) == (1,3) || (a,b) == (3,2) || (a,b) == (2,1)
        0x01
    elseif (a,b) == (1,2) || (a,b) == (2,3) || (a,b) == (3,1)
        0x03
    else
        0x00
    end
end
levicivita(x::@compat Tuple{UInt8,UInt8}) = levicivita(x...)
levicivita(a::Vector{UInt8}, b::Vector{UInt8}) = mapreduce(levicivita, +, zip(a,b))

# with our Pauli representation, multiplication is the sum (mod 4), or equivalently, the
# XOR of the bits
*(a::Pauli, b::Pauli) = Pauli(a.v $ b.v, a.s + b.s + levicivita(a.v, b.v))

const phases_ = [1, im, -1, -im]

function *(n::Number, p::Pauli)
    ns = findfirst(n .== phases_) - 1
    @assert(ns >= 0, "Multiplication only supported for +/- 1, +/- im")
    Pauli(p.v, p.s + ns)
end
*(p::Pauli, n::Number) = n * p
*(p::Pauli, u::Matrix) = *(promote(p, u)...)
*(u::Matrix, p::Pauli) = *(promote(u, p)...)

+(p::Pauli) = p
-(p::Pauli) = Pauli(p.v, p.s + 2)

abs(p::Pauli) = Pauli(p.v, 0)
phase(p::Pauli) = phases_[p.s+1]

length(p::Pauli) = length(p.v)
vec(p::Pauli) = vec(convert(Matrix{Complex{Int}}, p))

kron(a::Pauli, b::Pauli) = Pauli([a.v; b.v], a.s + b.s)

function expand(a::Pauli, subIndices, n)
    v = zeros(n)
    for (ct, i) in enumerate(subIndices)
        v[i] = a.v[ct]
    end
    Pauli(v, a.s)
end

function generators(a::Pauli)
    G = Pauli[]
    all(a.v .== 0) && return abs(a)
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
    @eval const $(Symbol(a*b)) = kron($ao, $bo)
end

function allpaulis(n)
    if n <= 0
        error("Need at least 1 qubit")
    elseif n==1
        return [Id,X,Y,Z]
    else
        return map(x->kron(x[1],x[2]),product([Id,X,Y,Z],allpaulis(n-1)))
    end
end
