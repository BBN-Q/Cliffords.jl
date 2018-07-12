# single qubit Cliffords

export localclifford

C1 = Vector{Clifford{1}}(24)

# identity
C1[1]  = RI

# pi/2 and pi rotations about X,Y,Z
C1[2]  = exp(-im*pi/4*X)
C1[3]  = im * X
C1[4]  = exp(im*pi/4*X)
C1[5]  = exp(-im*pi/4*Y)
C1[6]  = im * Y
C1[7]  = exp(im*pi/4*Y)
C1[8]  = exp(-im*pi/4*Z)
C1[9]  = im * Z
C1[10] = exp(im*pi/4*Z)

# Hadamard class
C1[11] = exp(-im*pi/2/sqrt(2) * (X+Y))
C1[12] = exp(-im*pi/2/sqrt(2) * (X-Y))
C1[13] = exp(-im*pi/2/sqrt(2) * (X+Z)) # standard Hadamard
C1[14] = exp(-im*pi/2/sqrt(2) * (X-Z))
C1[15] = exp(-im*pi/2/sqrt(2) * (Y+Z))
C1[16] = exp(-im*pi/2/sqrt(2) * (Y-Z))

# Axis exchange class
C1[17] = exp(-1im*pi/3/sqrt(3) * (X+Y+Z))
C1[18] = exp(-2im*pi/3/sqrt(3) * (X+Y+Z))
C1[19] = exp(-1im*pi/3/sqrt(3) * (X-Y+Z))
C1[20] = exp(-2im*pi/3/sqrt(3) * (X-Y+Z))
C1[21] = exp(-1im*pi/3/sqrt(3) * (X+Y-Z))
C1[22] = exp(-2im*pi/3/sqrt(3) * (X+Y-Z))
C1[23] = exp(-1im*pi/3/sqrt(3) * (-X+Y+Z))
C1[24] = exp(-2im*pi/3/sqrt(3) * (-X+Y+Z))

rC1 = Dict{Clifford{1},UInt}()

for (idx,v) in enumerate(C1)
    rC1[v] = idx
end

localclifford(i::Int) = C1[i]
localclifford(v::Vector) = kron(C1[v]...)

localcliffordindex(c::Clifford) = rC1[c]
