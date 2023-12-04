using PotentialFlow3D.StaticArrays

@testset "Vortex types" begin
    
    x = rand(3)
    p = VortexPoint3(x)
    @test p.coords == SVector{3}(x)

    x1, x2, x3, x4 = rand(3), rand(3), rand(3), rand(3)

    vl1 = VortexLoop([x1,x2,x3,x4],1.0)

    vl2 = unitloop([x1,x2,x3,x4])

    @test vl1.x == vl2.x && vl1.Γ == vl2.Γ && vl1.σ == vl2.σ

    @test x1 - vl1.x[2] isa SVector

    x1 = VortexPoint3([Inf,-Inf,5.0])
    @test x1.coords[1] == PotentialFlow3D.DEFAULTINF && x1.coords[2] == -PotentialFlow3D.DEFAULTINF && x1.coords[3] == 5.0

end

@testset "Basic operations" begin

    x1 = Float64[0,0,0]
    x2 = Float64[1,0,0]
    x3 = Float64[1,1,0]
    x4 = Float64[0,1,0]

    vl1 = VortexLoop([x1,x2,x3,x4],1.0)

    r1 = vec(segment(vl1,1))
    r2 = -vec(segment(vl1,4))
    @test cross(r1,r2) == SVector{3}(Float64[0,0,1])

    @test dot(r1,r2) == 0.0


end
