import Base.complex

export Pauli, Id, X, Y, Z, allpaulis, paulieye, weight, complex

# Paulis's are represented by an immutable vector of numbers (0-3) corresponding to
# single-qubit Paulis, along with a phase parameter.
immutable Pauli{N}
    v::Vec{N,UInt8} # 0 = I, 1 = X, 2 = Z, 3 = Y
    s::UInt8 # 0 = +1, 1 = +i, 2 = -1, 3 = -i (im^s)
    # function Pauli{N}(v::Vec{N,UInt8}, s::UInt8)
    #     new(v,mod(s,4))
    # end
end

Pauli{N}(v::Vec{N,UInt8}, s::Integer) = Pauli(v, convert(UInt8,s))
Pauli(v::Vector, s = 0) = Pauli{length(v)}(Vec{length(v),UInt8}(v), convert(UInt8, s))
Pauli(v::Integer, s = 0) = Pauli([v], s)
Pauli(m::Matrix) = convert(Pauli, m)

weight{N}(p::Pauli{N}) = sum( p.v .> 0 )

show{N}(io::IO, p::Pauli{N}) = print(io,convert(UTF8String,p))

=={N,M}(a::Pauli{N}, b::Pauli{M}) = (a.v == b.v && a.s == b.s)
isequal{N,M}(a::Pauli{N}, b::Pauli{M}) = (a == b)
hash{N}(a::Pauli{N}, h::UInt) = hash(a.v, hash(a.s, h))
isid{N}(a::Pauli{N}) = isempty(find(a.v))

function convert{N}(::Type{UTF8String}, p::Pauli{N})
    phases = ["+","i","-","-i"]
    paulis = "IXZY"
    phases[p.s+1] * join([paulis[i+1] for i in p.v])
end

function convert{T,N}(::Type{Matrix{Complex{T}}}, p::Pauli{N})
    const mats = @compat Dict(
        0x00 => eye(Complex{T},2),
        0x01 => Complex{T}[0 1; 1 0],
        0x02 => Complex{T}[1 0; 0 -1],
        0x03 => Complex{T}[0 -im; im 0])

    return phase(p)*reduce(kron,[mats[x] for x in p.v])
end

complex{N}(p::Pauli{N}) = convert(Matrix{Complex128},p)

function Pauli(m::Matrix)
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

function convert{N}(::Type{Pauli{N}}, m::Matrix)
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

promote_rule{T<:Real,N}(::Type{Pauli{N}}, ::Type{Matrix{T}}) = Matrix{Complex{T}}
promote_rule{T<:Complex,N}(::Type{Pauli{N}}, ::Type{Matrix{T}}) = Matrix{T}

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
levicivita(a, b) = mapreduce(levicivita, +, zip(a,b))

# with our Pauli representation, multiplication is the sum (mod 4), or equivalently, the 
# XOR of the bits
*(a::Pauli, b::Pauli) = Pauli(map($,a.v,b.v), a.s + b.s + levicivita(a.v, b.v))

const phases_ = [1, im, -1, -im]

function *{N}(n::Number, p::Pauli{N})
    ns = findfirst(n .== phases_) - 1
    @assert(ns >= 0, "Multiplication only supported for +/- 1, +/- im")
    Pauli(p.v, p.s + ns)
end
*{N}(p::Pauli{N}, n::Number) = n * p
*{T,N}(p::Pauli{N}, u::Matrix{T}) = *(promote(p, u)...)
*{T,N}(u::Matrix{T}, p::Pauli{N}) = *(promote(u, p)...)

+{N}(p::Pauli{N}) = p
-{N}(p::Pauli{N}) = Pauli(p.v, p.s + 0x02)

abs{N}(p::Pauli{N}) = Pauli(p.v, 0)
phase{N}(p::Pauli{N}) = phases_[p.s+1]

length{N}(p::Pauli{N}) = length(p.v)
vec{N}(p::Pauli{N}) = vec(convert(Matrix{Complex{Int}}, p))

kron{N,M}(a::Pauli{N}, b::Pauli{M}) = Pauli{N+M}(Vec{N+M,UInt8}([Vector{UInt8}(a.v); Vector{UInt8}(b.v)]), a.s + b.s)

function expand(a::Pauli, subIndices, n)
    v = zeros(n)
    for (ct, i) in enumerate(subIndices)
        v[i] = a.v[ct]
    end
    Pauli(v, a.s)
end

function generators{N}(a::Pauli{N})
    G = Pauli[]
    if all(map(Bool,a.v .== 0)) # hack because .== in FixedSizeArrays is broken
        return abs(a) 
    end
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

function allpaulis(n)
    if n <= 0
        error("Need at least 1 qubit")
    elseif n==1
        return [Id,X,Y,Z]
    else
	return map(x->kron(x[1],x[2]),product([Id,X,Y,Z],allpaulis(n-1)))
    end
end
