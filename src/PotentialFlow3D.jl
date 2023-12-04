module PotentialFlow3D

    using StaticArrays
    using LinearAlgebra
    using UnPack

    export VortexPoint3, VortexLineSegment, unitsegment, VortexLoop, segment, unitloop

    export cross, dot, norm

    include("types.jl")
    include("basicoperations.jl")
    include("induction.jl")

end
