module PotentialFlow3D

    using StaticArrays

    export VortexPoint3, VortexLineSegment, unitsegment, VortexLoop, segment, unitloop
    
    export cross

    include("types.jl")
    include("basicoperations.jl")

end
