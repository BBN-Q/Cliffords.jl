import Base: complex, getindex

export Pauli, Id, X, Y, Z, allpaulis, paulieye, weight, complex

# deal with symbols added or removed from Base in Julia 0.5/0.6
if isdefined(Base, :∘)
    import Base.∘
else
    export ∘
end
if isdefined(Base, :factor)
    import Base.factor
else
    export factor
end

# Paulis's are represented by an immutable vector of numbers (0-3) corresponding to
# single-qubit Paulis, along with a phase parameter.
struct Pauli{N}
    v::SVector{N,UInt8} # 0 = I, 1 = X, 2 = Z, 3 = Y
    s::UInt8 # 0 = +1, 1 = +i, 2 = -1, 3 = -i (im^s)
    Pauli{N}(v, s) where N = new(map(x->mod(x,0x4),v),mod(s,0x4))
end

Pauli(v::SVector{N,UInt8}, s) where {N} = Pauli{N}(v,s)
Pauli(v::MVector{N,UInt8}, s) where {N} = Pauli{N}(SVector{N}(v), s)
Pauli(v::Integer, s = 0) = Pauli{1}(SVector(v % UInt8), s)
Pauli(v::Vector, s = 0) = Pauli{length(v)}(SVector{length(v),UInt8}(v), s)
#Pauli(m::Matrix) = convert(Pauli{isqrt(size(m,1))}, m)

function Pauli(m::AbstractMatrix; tolfac=eps(Float64))
    d = size(m,1)
    n = log(2,size(m,1))
    for p in allpaulis(n)
        overlap = tr(m*p) / d
        if isapprox(abs(overlap),1,atol = tolfac * d)
            return (round(real(overlap))+im*round(imag(overlap)))∘p
        elseif !isapprox(abs(overlap),0,atol=d*eps(Float64))
            error("Trying to convert non-Pauli matrix to a Pauli object")
        end
    end
end

weight(p::Pauli) = Int(sum( p.v .> 0 ))

show(io::IO, p::Pauli) = print(io,convert(AbstractString,p))

==(a::Pauli, b::Pauli) = (a.v == b.v && a.s == b.s)
isequal(a::Pauli, b::Pauli) = (a == b)
hash(a::Pauli, h::UInt) = hash(a.v, hash(a.s, h))
isid(a::Pauli) = all(p == 0 for p in a.v)

"""
factor(p::Pauli)

Given an `N` qubit Pauli operation `p`, `factor` returns a list of
single-qubit pauli operations corresponding to the tensor product
factors.

"""
function factor(a::Pauli)
    factors = [ Pauli(a.v[i]) for i in 1:length(a) ]
    factors[1] = phase(a)∘factors[1]
    return factors
end

function getindex(a::Pauli,i)
    if i==1
        return phase(a)∘Pauli(a.v[i])
    else
        return Pauli(a.v[i])
    end
end

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
    pauli_lex = (1, 2, 4, 3)
    tuple([pauli_lex[x+1] for x in p.v]...)
end

function convert(::Type{AbstractString}, p::Pauli)
    phases = ["+","i","-","-i"]
    paulis = "IXZY"
    phases[p.s+1] * join([paulis[i+1] for i in p.v])
end

function convert(::Type{Matrix{Complex{T}}}, p::Pauli) where T
    mats = Dict(
        0x00 => Complex{T}[1 0; 0 1],
        0x01 => Complex{T}[0 1; 1 0],
        0x02 => Complex{T}[1 0; 0 -1],
        0x03 => Complex{T}[0 -im; im 0])

    return phase(p)*reduce(kron,[mats[x] for x in p.v])
end

complex(p::Pauli) = convert(Matrix{ComplexF64},p)

function convert(::Type{Pauli{N}}, m::Matrix) where N
    return Pauli(m)
end

promote_rule(::Type{Pauli{N}}, ::Type{Matrix{T}}) where {T<:Real,N} = Matrix{Complex{T}}
promote_rule(::Type{Pauli{N}}, ::Type{Matrix{T}}) where {T<:Complex,N} = Matrix{T}

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
levicivita(a, b) = reduce(+, levicivita.(a, b))

# with our Pauli representation, multiplication is the sum (mod 4), or equivalently, the
# XOR of the bits
*(a::Pauli{1}, b::Pauli{1}) = Pauli(a.v[1] ⊻ b.v[1],
                                    mod(a.s + b.s + levicivita(a.v[1], b.v[1]),4))
*(a::Pauli{2}, b::Pauli{2}) = Pauli{2}(SVector{2,UInt8}(a.v[1] ⊻ b.v[1], a.v[2] ⊻ b.v[2]),
                                       mod(a.s + b.s + levicivita(a.v, b.v),4))
*(a::Pauli{N}, b::Pauli{N}) where {N} = Pauli{N}(SVector{N,UInt8}(ntuple(i -> a.v[i] ⊻ b.v[i], N)),
                                                 mod(a.s + b.s + levicivita(a.v, b.v),4))

const phases_ = [1, im, -1, -im]
const phaseDict_ = Dict(1 => 0x0, im => 0x1, -1 => 0x2, -im => 0x3)
# special "multiplication" operator that returns a Pauli
function ∘(n::Number, p::Pauli)
    Pauli(p.v, mod(p.s + phaseDict_[n],4))
end

# unary operators return Paulis
+(p::Pauli) = p
-(p::Pauli) = Pauli(p.v, p.s + 0x02)

# standard binary operators are defined to return Matrix
*(n::Number, p::Pauli) = n * complex(p)
*(p::Pauli, n::Number) = n * p
*(p::Pauli, u::Matrix) = *(promote(p, u)...)
*(u::Matrix, p::Pauli) = *(promote(u, p)...)

+(p::Pauli, u::Matrix) = +(promote(p, u)...)
+(u::Matrix, p::Pauli) = +(promote(p, u)...)
+(a::Pauli, b::Pauli) = complex(a) + complex(b)

-(p::Pauli, u::Matrix) = -(promote(p, u)...)
-(u::Matrix, p::Pauli) = -(promote(p, u)...)
-(a::Pauli, b::Pauli) = complex(a) - complex(b)

abs(p::Pauli) = Pauli(p.v, 0)
phase(p::Pauli) = phases_[p.s+1]

length(p::Pauli{N}) where {N} = N
vec(p::Pauli) = vec(convert(Matrix{Complex{Int}}, p))

kron(a::Pauli{N}, b::Pauli{M}) where {N,M} = Pauli{N+M}([a.v; b.v], a.s + b.s)

expand(a::Pauli{1}, index::Number, ::Type{Val{1}}) = a
function expand(a::Pauli{1}, index::Number, ::Type{Val{N}}) where N
    v = SVector{N,UInt8}(ntuple(i -> i == index ? a.v[1] : 0x0, N))
    Pauli(v, a.s)
end

function expand(a::Pauli, subIndices::Vector, n)
    v = zeros(n)
    for (ct, i) in enumerate(subIndices)
        v[i] = a.v[ct]
    end
    Pauli(v, a.s)
end

function generators(a::Pauli{1})
    if abs(a) == Y
        return Pauli{1}[im*phase(a)∘X,Z]
    else
        return Pauli{1}[a]
    end
end

function generators(a::Pauli{N}) where N
    if isid(a)
        return abs(a)
    end
    G = Pauli{N}[]
    s = phase(a)
    for (idx, p) in enumerate(a.v)
        if p == 0 # I, skip it
            continue
        elseif p == 1 || p == 2 # X or Z
            push!(G, expand(Pauli(p), idx, Val{N}))
        else # Y
            s *= im
            push!(G, expand(X, idx, Val{N}))
            push!(G, expand(Z, idx, Val{N}))
        end
    end
    # give phase to first generator
    G[1] = s ∘ G[1]
    return G
end

function paulieye(n)
    expand(Id, 1, Val{n})
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
