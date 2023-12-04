@testset "Velocity" begin
    x1 = VortexPoint3([0.0,0.0,0.0])
    x2 = VortexPoint3([Inf,0.0,0.0])

    r = 5.0*rand()
    x = VortexPoint3([0.0,r,0.0])

    seg = VortexLineSegment(x1,x2,1.0)

    vel = PotentialFlow3D.segmentvelocity(x,seg)
    @test vel ≈ [0.0,0.0,1.0/(4π*r)]

    x1 = VortexPoint3([0.0,0.0,-Inf])
    x2 = VortexPoint3([0.0,0.0,Inf])
    Γ = rand()
    seg = VortexLineSegment(x1,x2,Γ)

    r = 5.0*rand()
    x = VortexPoint3([r,0.0,rand()])
    vel = PotentialFlow3D.segmentvelocity(x,seg)
    @test vel ≈ [0.0,Γ/(2π*r),0.0]

end