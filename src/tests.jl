include("delaunay.jl")

using Random, Test

Random.seed!(42)

@testset "Delaunay.jl" begin
	for _ = 1:100
		a = Vertex(rand(Float64), rand(Float64))
		b = Vertex(rand(Float64), rand(Float64))
		c = Vertex(0,0)
		p = Vertex(rand(Float64), rand(Float64))

		if a == b
			continue
		end

		L = is_left(p,a,b)

		if a.y == b.y
			@test (a.x < b.x) == (p.y > a.y)
		elseif a.x == b.x
			@test (a.y < b.y) == (p.x > a.x)
		else
			m = (b.y - a.y) / (b.x - a.x)
			y_g = m*(p.x - a.x) + a.y
			Lcomp = (p.y > y_g) == (b.x > a.x)
			@test Lcomp == L
		end
	end
end
