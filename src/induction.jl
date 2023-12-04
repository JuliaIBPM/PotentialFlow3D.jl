## Operations to evaluate velocity, vector potential, etc

function velocity(x::AbstractVector{T},vl::VortexLoop) where {T<:Real}
    vel = zeros(T,3)
    for j in 1:length(vl)
        vel += velocity(x,segment(vl,j))
    end
    return SVector{3,T}(vel)

end

function velocity(x::AbstractVector{T},vl::VortexLineSegment) where {T<:Real}
    @unpack σ, Γ = vl
    lensq = dot(vl,vl)
    r1 = x - vl.xstart
    r2 = x - vl.xend

    r1norm = sqrt(dot(r1,r1)+σ^2)
    r1dir = r1/r1norm

    r2norm = sqrt(dot(r2,r2)+σ^2)
    r2dir = r2/r2norm

    b = cross(r1,r2)
    bsq = dot(b,b) #+ σ*lensq

    # f = r1/|r1| - r2/|r2|
    f = r1dir - r2dir

    g = -dot(vec(vl),f)

    boverbsq = _boverbsq(b,bsq,lensq,Val(g==0.0))

    return -0.25/π*Γ*boverbsq*g


end


velocity(x::VortexPoint3,vl) = velocity(x.coords,vl)

_boverbsq(b::AbstractVector{T},bsq,lensq,::Val{true}) where {T} = SVector{3,T}([0.0,0.0,0.0])

_boverbsq(b::AbstractVector{T},bsq,lensq,::Val{false}) where {T} = b/(bsq + eps(T)^2*lensq)
