include("logic.jl")

using Random, Test, .Logic, .Logic.DCEL

Random.seed!(42)

#=								IS_LEFT WORKS.
@testset "is_left" begin
	for _ = 1:100
		a = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)
		b = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)
		p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)

		if a == b
			continue
		end

		L = is_left(p,a,b)
		# Alternate messy calculation using y=mx+b
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
=#

#=									INSERT_POINT WITHOUT FLIP WORKS.
@testset "Structural Test" begin
	D = Delaunay(Set{Triangle}())

	for num = 1:100
		p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)

		println("Attempting to insert:", p)
		D,_ = insert_point_no_flip!(p, D)

		println("Testing if triangles were added")
		@test length(D.triangles) == 1+2*num
		@test length(filter(T -> T.edge.origin == p || T.edge.next.origin == p || T.edge.next.next.origin == p, D.triangles)) > 2

		println("Testing if all triangles are circular connected")
		for T in D.triangles
			e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
			@test e2.prev === e1
			@test e3.prev === e2
			@test e3.next === e1
		end

		println("Testing if half-edges have correct origins")
		for T in D.triangles
			for e in (T.edge, T.edge.next, T.edge.prev)
				@test e isa Border || e.twin.origin === e.next.origin
			end
		end
	end
end
=#

@testset "Flip Test" begin
	D = Delaunay(Set{Triangle}())

	for num = 1:10
		p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)

		println("Attempting to insert:", p)
		D = insert_point!(p, D)
		#println("After insertion:", D)

		println("Testing if triangles were added")
		@test length(D.triangles) == 1+2*num
		@test length(filter(T -> T.edge.origin == p || T.edge.next.origin == p || T.edge.next.next.origin == p, D.triangles)) > 2

		println("Testing if all triangles are circular connected")
		for T in D.triangles
			e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
			@test e2.prev === e1
			@test e3.prev === e2
			@test e3.next === e1
		end

		println("Testing if half-edges have correct origins")
		for T in D.triangles
			for e in (T.edge, T.edge.next, T.edge.prev)
				@test e isa Border || e.twin.origin === e.next.origin
			end
		end

		println("Testing if Delaunay condition is upheld")
		for T in D.triangles
			for e in (T.edge, T.edge.next, T.edge.prev)
				@test is_delaunay(e)
			end
		end

	end
end