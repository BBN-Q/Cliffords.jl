using Cliffords, Test, LinearAlgebra

import Cliffords: IZ, ZI, XI, IX, YI, IY, II, XX, YY, ZZ, XY, XZ, YX, ZX,
        generators

# avoid annoying deprecation warnings about eye()
import LinearAlgebra: eye
eye(n::Integer) = Matrix{Float64}(I, n, n)
eye(T::Type, n::Integer) = Matrix{T}(I, n, n)

include("testpaulis.jl")
include("testcliffords.jl")
