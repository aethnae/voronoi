Tests = 100
Digits = 3

@testset "Triangulation, low level logic" begin
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
		p = randomVertex(Digits)
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

@testset "Area calculation test" begin
	println("Inserting points...")
	D = Delaunay()
	for _ = 1:Tests
		p = randomVertex(Digits)
		D = insert_point!(p, D)
	end

	println("Done. Converting to graph...")
	neighbors = Dict{Vertex, Set{Vertex}}()
	for T in D.triangles
		for e in (T.edge, T.edge.next, T.edge.prev)
			v = e.origin
			if !(v in (BottomLeft, BottomRight, TopLeft))
				neigh = haskey(neighbors, v) ? neighbors[v] : Set{Vertex}()
				push!(neigh, e.next.origin)
				push!(neigh, e.prev.origin)
				neighbors[v] = neigh
			end
		end
	end
	N = length(neighbors)

	println("$(N) vertices found. Generating polygons...")
	polygons = Dict{Vertex, Vector{Vertex}}()
	for v in keys(neighbors)
		neigh = sort([x for x in neighbors[v]], by = (p -> atan(p.y, p.x)))

		# Conversion to anticlockwise perpendicular bisectors x=m+nt
		perp = Vector{Tuple{Vertex,Vertex}}()
		for p in neigh
			if p == TopLeft
				push!(perp, (Vertex(0.0, 1.0), Vertex(-1.0, 0.0)))
				push!(perp, (Vertex(0.0, 1.0), Vertex(0.0, -1.0)))
			elseif p == BottomLeft
				push!(perp, (Vertex(0.0, 0.0), Vertex(0.0, -1.0)))
				push!(perp, (Vertex(0.0, 0.0), Vertex(1.0, 0.0)))
			elseif p == BottomRight
				push!(perp, (Vertex(1.0, 0.0), Vertex(1.0, 0.0)))
				push!(perp, (Vertex(1.0, 0.0), Vertex(0.0, 1.0)))
			else
				push!(perp, (0.5*(p+v), Vertex(p.y-v.y, v.x-p.x)))
			end
		end
		
		# intersections of perpendicular bisectors
		corners = Vector{Vertex}()
		for i in 1:length(perp)
			j = i % length(perp) + 1
			println("Center $(v) Intersect $(perp[i]), $(perp[j])")

			mij = perp[i][1]-perp[j][1]
			ni, nj = perp[i][2], perp[j][2]
			denom = ni.x * nj.y - ni.y * nj.x
			if denom == 0.0
				continue
			end
			ti = (nj.x * mij.y - nj.y * mij.x) / denom
			tj = (ni.x * mij.y - ni.y * mij.x) / denom
			@test ti >= 0
			@test tj <= 0
			u1 = perp[i][1] + ti * ni
			u2 = perp[i][2] + tj * nj
			@test u1.x == u2.x && u1.y == u2.y
			push!(corners, u1)
		end
		@test length(corners) > 2
		polygons[v] = corners
	end

	println("Calculating area...")
	a = 0.0
	for v in keys(polygons)
		av = 0.0
		corners = polygons[v]
		for i in length(corners)
			j = i % length(corners) + 1
			av = av + (corners[i].y+corners[j].y)*(corners[i].x-corners[j].y) / 2
		end
		println("$(av) for Vertex $(v)")
		a = a + av
	end
	println("$(a) in total")
end