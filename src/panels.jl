function surfacegrid_to_panels(xv,yv,zfunc)
    zv = zfunc.(xv,yv')
    panels = VortexLoop[]
    for j in 1:length(yv)-1, i in 1:length(xv)-1
       p1 = [xv[i],  yv[j],  zv[i,  j]]
       p2 = [xv[i+1],yv[j],  zv[i+1,j]]
       p3 = [xv[i+1],yv[j+1],zv[i+1,j+1]]
       p4 = [xv[i],  yv[j+1],zv[i,  j+1]]
       push!(panels,VortexLoop([p1,p2,p3,p4],1.0) )
    end
    return panels
end

function influence_matrix(panels::Vector{<:VortexLoop})
    A = [dot(velocity(loopcenter(panels[i]),panels[j]),normal(panels[i])) for i in eachindex(panels), j in eachindex(panels)]
end

"""
    panelrhs(vort::Vector{<:VortexLoop},panels::Vector{<:VortexLoop})

Return -(vw + vw_img)⋅n evaluated at the panel centers.
"""
function panelrhs(vort::Vector{<:VortexLoop},panels::Vector{<:VortexLoop})
    b = [-dot(velocity(loopcenter(panels[i]),vort),normal(panels[i])) for i in eachindex(panels)]
end

function loops2loops(target::Vector{<:VortexLoop},source::Vector{<:VortexLoop})
    return [dot(velocity(loopcenter(target[i]),source),normal(target[i])) for i in eachindex(target)]
end

function assign_panel_strengths(panels::Vector{<:VortexLoop},Γ::Vector{<:Real})
    new_panels = VortexLoop[]
    for (j,p) in enumerate(panels)
       push!(new_panels,VortexLoop(p.x,Γ[j]) )
    end
    return new_panels
end