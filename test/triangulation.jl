@testset "Triangulation" begin

	function test_low_level_logic(Tests::Int = 100, Digits::Int = 3)
		println("Starting $(Tests) low level logic (is_left) tests.")
		for _ = 1:Tests
			a = randomVertex(Digits)
			b = randomVertex(Digits)
			p = randomVertex(Digits)

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

			@assert Lcomp == L "
				Incorrect result for p=$(p) relative to $(a)->$(b)."
		end
		println("Finished $(Tests) low level logic (is_left) tests.")
	end

	function test_high_level_logic(Tests::Int = 100, Digits::Int = 3)
		println("Starting $(Tests) full logic tests.")

		D = Delaunay()

		for num = 1:Tests
			p = randomVertex(Digits)
			D = insert_point!(p, D)

			@assert length(D.triangles) == 1 + 2*num "
				Expected 1+2*$(num) triangles, got $(length(D.triangles))"

			p_found = filter(T -> T.edge.origin == p || T.edge.next.origin == p || T.edge.next.next.origin == p, D.triangles)
			@assert length(p_found) >= 3 "
				Expected at least 3 triangles containing inserted point, found only $(p_found)"

			for T in D.triangles
				e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
				@assert e2.prev === e1 && e2.prev === e1 && e3.prev === e2 && e3.next === e1 "
					Triangle $(T) is not connected correctly after insert"
			end

			for T in D.triangles
				twin_fail = filter(e -> !(e isa Border) && e.twin.origin !== e.next.origin, (T.edge, T.edge.next, T.edge.prev))
				@assert isempty(twin_fail) "
					Twins $(twin_fail[1]) and $(twin_fail[1].twin) not connected correctly after insert"
			end

			for T in D.triangles
				delaunay_fail = filter(e -> !is_delaunay(e), (T.edge, T.edge.next, T.edge.prev))
				@assert isempty(delaunay_fail) "
					Triangle $(T) breaks Delaunay condition at edge $(delaunay_fail[1])"
			end
		end
		println("Finished $(Tests) full logic tests.")
	end
end