@testset "Voronoi" begin
    println("Starting $(Tests) Voronoi tests.")

    D = Delaunay()
    
    for _ = 1:Tests
        p = round(Vertex(rand(Float64),rand(Float64)))
        D = insert_point!(p, D)

        V = voronoi_2(D)
        A = areas_2(V)
        S = sum(keys(A)) do v
            A[v]
        end

        @test isapprox(S, 1; atol = 0.01)
    end
end