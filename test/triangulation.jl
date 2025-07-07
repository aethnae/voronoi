@testset "Triangulation, low level logic" begin
	println("Starting $(Tests) low level logic (is_left) tests.")
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

		@test Lcomp == L
		if Lcomp != L
			println("Incorrect result for p=$(p) relative to $(a)->$(b).")
		end
	end
	println("Finished $(Tests) low level logic (is_left) tests.")
end

@testset "Triangulation, high level logic" begin
	println("Starting $(Tests) full logic tests.")

	D = Delaunay()

	for num = 1:Tests
		p = randomVertex()
		D = insert_point!(p, D)

		@test length(D.triangles) == 1 + 2*num
		if length(D.triangles) != 1 + 2*num 
			println("Expected 1+2*$(num) triangles, got $(length(D.triangles))")
		end

		p_found = filter(T -> T.edge.origin == p || T.edge.next.origin == p || T.edge.next.next.origin == p, D.triangles)
		@test length(p_found) >= 3 
		if !(length(p_found) >= 3 )
			println("Expected at least 3 triangles containing inserted point, found only $(p_found)")
		end

		for T in D.triangles
			e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
			@test e2.prev === e1 && e2.prev === e1 && e3.prev === e2 && e3.next === e1 
			if !(e2.prev === e1 && e2.prev === e1 && e3.prev === e2 && e3.next === e1 )
				println("Triangle $(T) is not connected correctly after insert")
			end
		end

		for T in D.triangles
			twin_fail = filter(e -> !(e isa Border) && e.twin.origin !== e.next.origin, (T.edge, T.edge.next, T.edge.prev))
			@test isempty(twin_fail) 
			if !(isempty(twin_fail) )
				println("Twins $(twin_fail[1]) and $(twin_fail[1].twin) not connected correctly after insert")
			end
		end

		for T in D.triangles
			delaunay_fail = filter(e -> !is_delaunay(e), (T.edge, T.edge.next, T.edge.prev))
			@test isempty(delaunay_fail) 
			if !(isempty(delaunay_fail) )
				println("Triangle $(T) breaks Delaunay condition at edge $(delaunay_fail[1])")
			end
		end
	end
	println("Finished $(Tests) full logic tests.")
end