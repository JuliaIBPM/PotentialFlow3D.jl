@testset "Vortex types" begin
    
    x = rand(3)
    p = VortexPoint3(x)

    x1, x2, x3, x4 = rand(3), rand(3), rand(3), rand(3)

    vl1 = VortexLoop([x1,x2,x3,x4],1.0)

    vl2 = unitloop([x1,x2,x3,x4])

    @test vl1.x == vl2.x && vl1.Γ == vl2.Γ && vl1.σ == vl2.σ
end
