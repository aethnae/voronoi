include("logic.jl")

using Random, Test, .Logic, .Logic.DCEL

Random.seed!(42)

Tests = 100
Digits = 3

function randomVertex()
	return Vertex(round(rand(Float64)*10^Digits)/(10^Digits), round(rand(Float64)*10^Digits)/(10^Digits))
end

@testset "is_left" begin
	println("Starting $(Tests) is_left tests.")
	for _ = 1:Tests
		a = randomVertex()
		b = randomVertex()
		p = randomVertex()

		if a == b
			continue
		end

		L = is_left(p,a,b)
		Lcomp = nothing

		# Alternate messy calculation using y=mx+b
		if a.y == b.y
			Lcomp = (a.x < b.x) == (p.y > a.y)
		elseif a.x == b.x
			Lcomp = (a.y < b.y) == (p.x > a.x)
		else
			m = (b.y - a.y) / (b.x - a.x)
			y_g = m*(p.x - a.x) + a.y
			Lcomp = (p.y > y_g) == (b.x > a.x)
		end

		if Lcomp != L
			@test "Incorrect result for p=$(p) relative to $(a)->$(b)" == nothing
		end
	end
end

@testset "Full Logic Test" begin
	println("Starting $(Tests) full logic tests.")

	D = Delaunay()

	for num = 1:Tests
		p = randomVertex()
		D = insert_point!(p, D)

		if length(D.triangles) != 1 + 2*num
			@test "Expected 1+2*$(num) triangles, got $(length(D.triangles))" == nothing
		end
		p_found = length(filter(T -> T.edge.origin == p || T.edge.next.origin == p || T.edge.next.next.origin == p, D.triangles))
		if p_found < 3
			@test "Expected at least 3 triangles containing inserted point, found $(p_found)" == nothing
		end

		for T in D.triangles
			e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
			if e2.prev !== e1 || e2.prev !== e1 || e3.prev !== e2 || e3.next !== e1
				@test "Triangle $(T) is not connected correctly after insert" == nothing
			end
		end

		for T in D.triangles
			twin_fail = filter(e -> !(e isa Border) && e.twin.origin !== e.next.origin, (T.edge, T.edge.next, T.edge.prev))
			if !isempty(twin_fail)
				@test "Twins $(twin_fail[1]) and $(twin_fail[1].twin) not connected correctly after insert" == nothing
			end
		end

		for T in D.triangles
			delaunay_fail = filter(e -> !is_delaunay(e), (T.edge, T.edge.next, T.edge.prev))
			if !isempty(delaunay_fail)
				@test "Triangle $(T) breaks Delaunay condition at edge $(delaunay_fail[1])" == nothing
			end
		end

	end
end