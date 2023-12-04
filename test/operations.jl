@testset "Velocity" begin
    x1 = VortexPoint3([0.0,0.0,0.0])
    x2 = VortexPoint3([Inf,0.0,0.0])

    r = 5.0*rand()
    x = [0.0,r,0.0]

    seg = VortexLineSegment(x1,x2,1.0)

    vel = velocity(x,seg)
    @test vel ≈ [0.0,0.0,1.0/(4π*r)]

    x1 = VortexPoint3([0.0,0.0,-Inf])
    x2 = VortexPoint3([0.0,0.0,Inf])
    Γ = rand()
    seg = VortexLineSegment(x1,x2,Γ)

    r = 5.0*rand()
    x = [r,0.0,rand()]
    vel = velocity(x,seg)
    @test vel ≈ [0.0,Γ/(2π*r),0.0]

    x = [0.0,0.0,0.0]
    vel = velocity(x,seg)
    @test vel ≈ [0.0,0.0,0.0]

    x1 = VortexPoint3(rand(3))
    x2 = VortexPoint3(rand(3))
    seg = VortexLineSegment(x1,x2,Γ)
    vel = velocity(x1,seg)
    @test vel ≈ [0.0,0.0,0.0]

    x1, x2, x3, x4 = rand(3), rand(3), rand(3), rand(3)
    vl1 = VortexLoop([x1,x2,x3,x4],1.0)

    x = [0.0,1.0,0.0]

    vel = velocity(x,vl1)

    vel1 = velocity(x,VortexLineSegment(VortexPoint3(x1),VortexPoint3(x2),1.0))
    vel2 = velocity(x,VortexLineSegment(VortexPoint3(x2),VortexPoint3(x3),1.0))
    vel3 = velocity(x,VortexLineSegment(VortexPoint3(x3),VortexPoint3(x4),1.0))
    vel4 = velocity(x,VortexLineSegment(VortexPoint3(x4),VortexPoint3(x1),1.0))
    @test vel1+vel2+vel3+vel4 ≈ vel
end