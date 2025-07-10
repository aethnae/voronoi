@testset "Voronoi" begin
    println("Starting $(Tests) Voronoi tests.")
    bbox=[Vertex(0.0,0.0), Vertex(1.0,0.0), Vertex(1.0,1.0), Vertex(0.0,1.0)]
    D = Delaunay()
    
    for i = 1:Tests
        p = round(Vertex(rand(Float64),rand(Float64), mod1(i,2)))
        D = insert_point!(p, D)

        #V = voronoi(D,bbox)        hat +10% Abweichung
        V = voronoi_2(D)
        A = areas(V)
        S = sum(keys(A)) do v
            A[v]
        end

        @test isapprox(S, 1; atol = 0.01)
    end
end