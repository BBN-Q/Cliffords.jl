using Cliffords, Test, LinearAlgebra

import Cliffords: IZ, ZI, XI, IX, YI, IY, II, XX, YY, ZZ, XY, XZ, YX, ZX,
        generators

# avoid an annoying deprecation warning
import LinearAlgebra: eye
eye(n::Integer) = Matrix{Float64}(I, n, n)

include("testpaulis.jl")
include("testcliffords.jl")
