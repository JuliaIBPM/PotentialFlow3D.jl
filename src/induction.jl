## Operations to evaluate velocity, vector potential, etc

const HTHRESH = 0.9

"""
    velocity(x::AbstractVector,v::Vector{<:VortexLoop}) -> SVector{3}

Return the velocity induced at point `x` by a vector of vortex loops.
"""
function velocity(x::AbstractVector{T},vl::Vector{<:VortexLoop}) where {T<:Real}
    vel = zeros(T,3)
    for v in vl
        vel += velocity(x,v)
    end
    return SVector{3,T}(vel)
end

"""
    velocity(x::AbstractVector,v::VortexLoop) -> SVector{3}

Return the velocity induced at point `x` by vortex loop `v`.
"""
function velocity(x::AbstractVector{T},vl::VortexLoop) where {T<:Real}
    vel = zeros(T,3)
    for j in 1:length(vl)
        vel += velocity(x,segment(vl,j))
    end
    return SVector{3,T}(vel)

end

"""
    velocity(x::AbstractVector,v::VortexLineSegment) -> SVector{3}

Return the velocity induced at point `x` by vortex line segment `v`.
"""
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


function vectorpotential(x::AbstractVector{T},vl::VortexLineSegment) where {T<:Real}
    @unpack σ, Γ = vl
    
    vvec = vec(vl)
    τ = vvec/norm(vvec)
    
    r1 = x - vl.xstart
    r2 = x - vl.xend

    r1norm = sqrt(dot(r1,r1)+σ^2)
    r1dir = r1/r1norm

    r2norm = sqrt(dot(r2,r2)+σ^2)
    r2dir = r2/r2norm

    cosβ1 = dot(r1dir,τ)
    cosβ2 = dot(r2dir,τ)
    println([cosβ1,cosβ2])

    return Γ*τ/(8π)*log((1-cosβ2)*(1+cosβ1)/((1+cosβ2)*(1-cosβ1)))


end

function vectorpotential(vs::VortexLineSegment{XT},vt::VortexLineSegment{XT}) where {XT}
    σssq, σtsq, σsσt = vs.σ^2, vt.σ^2, vs.σ*vt.σ
    Γs, Γt = vs.Γ, vt.Γ
    
    ls = vec(vs)
    lt = vec(vt)

    cosα = dot(ls,lt)
    lens = norm(vs)

    σsσt = 1.0e-5^2*lens^2


    cosα == 0 && return zero(XT)

    lent = norm(vt)
    es = ls/lens # unit vector along source segment

    cosα /= lens*lent

    sinαsq = 1.0 - cosα^2
    sinαsq = abs(sinαsq) < 100eps(1.0) ? 0.0 : sinαsq

    r00 = vt.xstart - vs.xstart  # target start - source start
    r01 = r00 + lt               # target end - source start
    r10 = r00 - ls               # target start - source end
    r11 = r01 - ls               # target end - source end

    # Lengths of relative position vectors
    lenr = zeros(Float64,2,2)
    lenr[1,1] = sqrt(dot(r00,r00) + σsσt)
    lenr[2,1] = sqrt(dot(r01,r01) + σsσt)
    lenr[1,2] = sqrt(dot(r10,r10) + σsσt)
    lenr[2,2] = sqrt(dot(r11,r11) + σsσt)

    et = copy(es)

    # Projections of relative position vectors onto source axis
    esr = zeros(Float64,2,2)
    esr[1,1] = dot(es,r00)
    esr[2,1] = esr[1,1] + lent*cosα
    esr[1,2] = esr[1,1] - lens
    esr[2,2] = esr[2,1] - lens

    # Integration limits I = (I1(ulim(2,2))-I1(ulim(1,2))) -
    #                        (I0(ulim(2,1))-I0(ulim(1,1)))

    # Half-space check. If any of these are greater than a critical
    # value (~0.9), then they are too close to the positive axis of the source
    # segment and we will use the other branch of the solution there.
    hspace = esr./lenr
    hmax = maximum(hspace)

    ulim = hmax > HTHRESH ? lenr + esr : lenr - esr

    sgns = [1 -1; -1 1]
    if sinαsq == 0.0
        # Parallel segments
        inta = lenr .+ (lenr .- ulim).*log.(ulim)
        
        int = sum(inta.*sgns)/(4π)

        return int
    end
    sinα = sqrt(sinαsq)

    et = lt/lent # Unit vector along target segment
    Ls = I - es*es' # Perpendicular projectors 
    Lt = I - et*et'
    d = [Lt*r00, Lt*r10]

    lend = norm.(d) # Perpendicular distances from ends of source to target line axis
    ns = Ls*et/sinα # Orthogonal coordinate system
    nt = Lt*es/sinα
    eperp = cross(es,ns)

    dβ = [dot(es,di) for di in d] # Projections of orthogonal components on source axis
    csq = lend.^2 .* sinαsq .- dβ.^2
    csq[abs.(csq) .< 100eps(1.0)] .= 0.0
    c = sqrt.(csq)
    fact = cosα/sinαsq

    # Projections of relative position vectors onto source axis
    etr = zeros(Float64,2,2)
    etr[1,1] = dot(et,r00)
    etr[2,1] = dot(et,r01)
    etr[1,2] = dot(et,r10)
    etr[2,2] = dot(et,r11)
    
    # lam+- = sqrt(u^2 +- 2|d|βu - c^2)
    hsign = 1 - 2(hmax <= 0.9)
    ulimt = lenr .+ hsign.*etr
    lam = cosα.*lenr .+ hsign*etr
    dβ .*= -hsign

    inta = (cosα.*ulim .- lam).*log.(ulim)
   
    tanterm = atan.(-c' .+ dβ'.*ulim,lam.*c')
    inta .+= dβ'.*log.(ulimt) .+ c'.*tanterm

    int = fact*sum(inta.*sgns)/(4π)
    return int

end




#=


%% Compute the integral.
cnt = 0;
if (omalfsq > 0),
    
    
    
    osgn = 1;
    for jj = 2:-1:1
        % jj represents the point on source axis (1 = lower point, 2 = upper point)
        isgn = 1;
        dbj = dbeta(jj);
        cj = c(jj);
        for ii = 2:-1:1
            % ii represents the limit on target axis (1 = lower , 2 = upper)
            u = ulim(ii,jj);
            ut = ulimt(ii,jj);
            lamij = lam(ii,jj);
            tanterm = atan2(-cj^2+dbj*u,cj*lamij); % Fix sign on dbj??            
            
            %logu = 0.5*log(u^2+sigsq);
            %logut = 0.5*log(ut^2+sigsq);
            logu = log(u);
            logut = log(ut);
            
            
            temp = (alf*u-lamij)*logu+dbj*logut+cj*tanterm;
            I = I + isgn*osgn*fact*temp;
            isgn = -isgn;
        end
        osgn = -osgn;
    end

end

I = I/(4*pi);

=#