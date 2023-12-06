module PotentialFlow3D

    using StaticArrays
    using LinearAlgebra
    using UnPack

    export VortexPoint3, VortexLineSegment, unitsegment, VortexLoop, segment, unitloop

    export cross, dot, norm, velocity, midpoint, normal, loopcenter

    export surfacegrid_to_panels, influence_matrix, panelrhs, assign_panel_strengths

    include("types.jl")
    include("basicoperations.jl")
    include("induction.jl")
    include("panels.jl")

end
