import Base: +, -, vec, length
import LinearAlgebra: cross, dot, norm

const DEFAULTBLOB = 1.0e-3
const DEFAULTINF = 1.0e12

struct VortexPoint3{T}
    coords :: SVector{3,T}
    VortexPoint3(x::SVector{3,T}) where {T} = new{T}(map(xi -> isinf(xi) ? sign(xi)*DEFAULTINF : xi, x)) 
end
Base.show(io::IO, p::VortexPoint3) = print(io, "point($(p.coords))")

VortexPoint3(p::AbstractVector{T}) where {T} = VortexPoint3(SVector{3,T}(p))


(+)(p1::VortexPoint3,p2::VortexPoint3) = p1.coords.+p2.coords
(-)(p1::VortexPoint3,p2::VortexPoint3) = p1.coords.-p2.coords
(-)(p::VortexPoint3) = VortexPoint3(-p1.coords)
(-)(x::AbstractVector{T},p::VortexPoint3{T}) where {T} = SVector{3,T}(x .- p.coords)

cross(p1::VortexPoint3,p2::VortexPoint3) = cross(p1.coords,p2.coords)
dot(p1::VortexPoint3,p2::VortexPoint3) = dot(p1.coords,p2.coords)

struct VortexLineSegment{XT,GT,BT}
    xstart :: VortexPoint3{XT}
    xend :: VortexPoint3{XT}
    Γ :: GT
    σ :: BT
end

"""
    VortexLineSegment(xstart::VortexPoint3,xend::VortexPoint3,Γ::Real; [blobradius::Real = 1e-3])

Create a line segment from `xstart` to `xend`, with strength `Γ` and blob radius `blobradius`
"""
VortexLineSegment(xstart::VortexPoint3,xend::VortexPoint3,Γ;blobradius::Real=DEFAULTBLOB) = VortexLineSegment(xstart,xend,Γ,blobradius)

unitsegment(xstart::VortexPoint3,xend::VortexPoint3;blobradius::Real=DEFAULTBLOB) = VortexLineSegment(xstart,xend,1.0;blobradius=blobradius)

function Base.show(io::IO, p::VortexLineSegment)
    println(io, "Vortex line segment from $(p.xstart) -> $(p.xend)")
    println(io, "  circulation: $(p.Γ)")
    println(io, "  blob radius: $(p.σ)")
end

"""
    vec(v::VortexLineSegment)

Return a vector `v.xend - v.xstart` for line segment `v`.
"""
vec(v::VortexLineSegment) = v.xend - v.xstart

"""
    cross(v1::VortexLineSegment,v2::VortexLineSegment)

Return the cross product between line segments `v1` and `v2`, where
the line segment vectors extend from start point to end point.
"""
cross(v1::VortexLineSegment,v2::VortexLineSegment) = cross(vec(v1),vec(v2))

"""
    dot(v1::VortexLineSegment,v2::VortexLineSegment)

Return the dot product between line segments `v1` and `v2`, where
the line segment vectors extend from start point to end point.
"""
dot(v1::VortexLineSegment,v2::VortexLineSegment) = dot(vec(v1),vec(v2))

"""
    norm(v::VortexLineSegment)

Return the length of line segment `v`
"""
norm(v::VortexLineSegment) = norm(vec(v))

### Vortex loop ###

struct VortexLoop{XT,GT,BT}
    x :: Vector{VortexPoint3{XT}}
    Γ :: GT
    σ :: BT
end

"""
    VortexLoop(x::Vector{VortexPoint3},Γ::Real[;blobradius=1e-3])

Create a loop of singular vorticity, with vertices `x`, of circulation `Γ`.
"""
function VortexLoop(x::Vector{T},Γ::Real;blobradius=DEFAULTBLOB) where {T<:VortexPoint3}
    n = length(x)
    VortexLoop(x,Γ,blobradius)
end

VortexLoop(x::Vector{T},Γ::Real;kwargs...) where T <: AbstractVector{S} where S <: Real = VortexLoop(VortexPoint3.(x),Γ;kwargs...)


function Base.show(io::IO, v::VortexLoop)
    println(io, "Vortex loop with $(length(v.x)) points:")
    println(io, "  circulation: $(v.Γ)")
    println(io, "  blob radius: $(v.σ)")
end

"""
    unitloop(x::Vector{Vector{Real}}[,blobradius=1e-3])

Create a vortex loop with unit strength with vertices `x`.
"""
unitloop(x::Vector{T};blobradius=DEFAULTBLOB) where {T<:VortexPoint3} = VortexLoop(x,1.0,blobradius)
unitloop(x::Vector{T};kwargs...) where T <: AbstractVector = VortexLoop(x,1.0;kwargs...)


"""
    segment(v::VortexLoop,j::Int)

Return the `j`th line segment on vortex loop `v`.
"""
segment(v::VortexLoop,j::Int) = VortexLineSegment(v.x[j],v.x[mod(j,length(v.x))+1],v.Γ,v.σ)

"""
    length(v::VortexLoop)

Number of segments in loop `v`
"""
length(v::VortexLoop) = length(v.x)