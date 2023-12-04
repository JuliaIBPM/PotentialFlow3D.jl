## Operations to evaluate velocity, vector potential, etc

function segmentvelocity(x::AbstractVector{T},vl::VortexLineSegment) where {T<:Real}
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


segmentvelocity(x::VortexPoint3,vl::VortexLineSegment) = segmentvelocity(x.coords,vl)

_boverbsq(b::AbstractVector{T},bsq,lensq,::Val{true}) where {T} = SVector{3,T}([0.0,0.0,0.0])

_boverbsq(b::AbstractVector{T},bsq,lensq,::Val{false}) where {T} = b/(bsq + eps(T)^2*lensq)

#=
rj = xeval(je,:) - x{1};
            %xjnorm = norm(rj);
            xjnorm = sqrt(rj*rj' + sigsq);
            xjdir = rj/xjnorm;
            
            % Loop through each segment formed by vertex pair
            for j = 1:nvert
                
                
                xjp1j  = x{j+1} - x{j};
                lensq = xjp1j*xjp1j';
                
                rjp1 = xeval(je,:) - x{j+1};
                xjp1norm = sqrt(rjp1*rjp1' + sigsq);
                xjp1dir = rjp1/xjp1norm;
                
                % b = rj x rj+1
                bj(1) = rj(2)*rjp1(3)-rj(3)*rjp1(2);
                bj(2) = rj(3)*rjp1(1)-rj(1)*rjp1(3);
                bj(3) = rj(1)*rjp1(2)-rj(2)*rjp1(1);
                bjsq = bj*bj'+sigsq*lensq;
           
                % f = rj/|rj| - rj+1/|rj+1|
                fj = xjdir - xjp1dir;
                                
                gj = -xjp1j*fj';
                if (gj == 0),
                    % Evaluation point lies on axis of this segment, but not on
                    % segment itself
                    boverbsq = zeros(1,3);
                    
                elseif (bjsq < eps),
                    % Point lies on this segment
                    bjsq = bjsq + eps*eps*lensq;
                    boverbsq = bj/bjsq;
                    
                else
                    boverbsq = bj/bjsq;
                    
                end
                
                % Velocity
                vseg = -boverbsq*gj;
                vtemp = vtemp + vseg;
                
                % Update for next vertex
                rj = rjp1;
                xjdir = xjp1dir;
                
                
            end
=#