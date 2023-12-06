## Basic operations on vortex elements

"""
    normal(v::VortexLoop)

Return the unit normal for the vortex loop `v`, averaged over the area inscribed by the loop.
The direction of the normal is determined by the order of the vertices, assuming a right-hand rule.
"""
function normal(vl::VortexLoop{XT}) where {XT}
    sum = zeros(3)
    for j in 1:length(vl)
        segj = segment(vl,j)
        sum .+= cross(midpoint(segj),vec(segj))
    end
    return SVector{3,XT}(sum/norm(sum))
end 

"""
    center(v::VortexLoop)

Return the center of the vortex loop, based on an arc-average of the points.
"""
function loopcenter(vl::VortexLoop{XT}) where {XT}
    sumn = zeros(3)
    sumd = 0.0
    for j in 1:length(vl)
        segj = segment(vl,j)
        sumn .+= norm(segj).*midpoint(segj)
        sumd += norm(segj)
    end
    return SVector{3,XT}(sumn/sumd)
end

                                                                

