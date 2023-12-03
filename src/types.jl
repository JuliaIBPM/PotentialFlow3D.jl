import Base: +, -

const DEFAULTBLOB = 1.0e-3

struct VortexPoint3{T}
    coords :: SVector{3,T}
end
Base.show(io::IO, p::VortexPoint3) = print(io, "point($(p.coords))")

VortexPoint3(p::AbstractVector{T}) where {T} = VortexPoint3(SVector{3,T}(p))


(+)(p1::VortexPoint3,p2::VortexPoint3) = VortexPoint3(p1.coords.+p2.coords)
(-)(p1::VortexPoint3,p2::VortexPoint3) = VortexPoint3(p1.coords.-p2.coords)
(-)(p::VortexPoint3) = VortexPoint3(-p1.coords)
(-)(x::AbstractVector{T},p::VortexPoint3{T}) where {T} = x .- p.coords

cross(p1::VortexPoint3,p2::VortexPoint3)= cross(p1.coords,p2.coords)

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


mutable struct VortexLoop{XT,GT,BT}
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

VortexLoop(x::Vector{T},Γ::Real;kwargs...) where T <: AbstractVector{S} where S <: Real = VortexLoop(VortexPoint3{S}.(x),Γ;kwargs...)


function Base.show(io::IO, v::VortexLoop)
    println(io, "Vortex loop with $(length(v.x)) points:")
    println(io, "  circulation: $(v.Γ)")
    println(io, "  blob radius: $(v.σ)")
end

unitloop(x::Vector{T};blobradius=DEFAULTBLOB) where {T<:VortexPoint3} = VortexLoop(x,1.0,blobradius)
unitloop(x::Vector{T};kwargs...) where T <: AbstractVector = VortexLoop(x,1.0;kwargs...)


"""
    segment(v::VortexLoop,j::Int)

Return the `j`th line segment on vortex loop `v`.
"""
segment(v::VortexLoop,j::Int) = VortexLineSegment(v.x[j],v.x[mod(j,length(v.x))+1],v.Γ,v.σ) 