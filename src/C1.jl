# single qubit Cliffords

const pX = [0 1; 1 0]
const pY = [0 -im; im 0]
const pZ = [1 0; 0 -1]

C1 = zeros(Clifford, 24)

# identity
C1[1]  = RI

# pi/2 and pi rotations about X,Y,Z
C1[2]  = expm(-im*pi/4*pX)
C1[3]  = im * pX
C1[4]  = expm(im*pi/4*pX)
C1[5]  = expm(-im*pi/4*pY)
C1[6]  = im * pY
C1[7]  = expm(im*pi/4*pY)
C1[8]  = expm(-im*pi/4*pZ)
C1[9]  = im * pZ
C1[10] = expm(im*pi/4*pZ)

# Hadamard class
C1[11] = expm(-im*pi/2/sqrt(2) * (pX+pY))
C1[12] = expm(-im*pi/2/sqrt(2) * (pX-pY))
C1[13] = expm(-im*pi/2/sqrt(2) * (pX+pZ)) # standard Hadamard
C1[14] = expm(-im*pi/2/sqrt(2) * (pX-pZ))
C1[15] = expm(-im*pi/2/sqrt(2) * (pY+pZ))
C1[16] = expm(-im*pi/2/sqrt(2) * (pY-pZ))

# Axis exchange class
C1[17] = expm(-1im*pi/3/sqrt(3) * (pX+pY+pZ))
C1[18] = expm(-2im*pi/3/sqrt(3) * (pX+pY+pZ))
C1[19] = expm(-1im*pi/3/sqrt(3) * (pX-pY+pZ))
C1[20] = expm(-2im*pi/3/sqrt(3) * (pX-pY+pZ))
C1[21] = expm(-1im*pi/3/sqrt(3) * (pX+pY-pZ))
C1[22] = expm(-2im*pi/3/sqrt(3) * (pX+pY-pZ))
C1[23] = expm(-1im*pi/3/sqrt(3) * (-pX+pY+pZ))
C1[24] = expm(-2im*pi/3/sqrt(3) * (-pX+pY+pZ))
