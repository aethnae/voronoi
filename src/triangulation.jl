module Triangulation

include("dcel.jl")
using .DCEL, LinearAlgebra

export insert_point!, test_low_level_logic, test_high_level_logic

#================================ LOGIC ===============================================#

"""
	is_delaunay(ab::Edge)::Bool

Checks if the Delaunay condition is upheld for this edge.
Borders always fulfill the Delaunay condition. For a half-edge
AB in triangle ABC, the point D of triangle BAD must not be in
the circumcircle of ABC.
"""
function is_delaunay(ab::Edge)::Bool
	if ab isa Border
		return true
	end
	ba = ab.twin
	a, b, c, d = ab.origin, ab.next.origin, ab.prev.origin, ba.prev.origin

	return det([a.x a.y a.x^2+a.y^2 1; b.x b.y b.x^2+b.y^2 1; c.x c.y c.x^2+c.y^2 1; d.x d.y d.x^2+d.y^2 1]) <= 0
end

"""
	is_left(p::Vertex, a::Vertex, b::Vertex)::Bool

Uses dot product to determine if a point p is left of a
directed edge e. If the point lies exactly on the edge,
it is not considered left of the edge.
"""
function is_left(p::Vertex, a::Vertex, b::Vertex)::Bool
	return (p.x - a.x) * (a.y - b.y) + (p.y - a.y) * (b.x - a.x) > 0
end


"""
    find_triangle(p::Vertex, D::Delaunay)::Triangle

Finds the triangle containing the given vertex. The vertex
is to the left of each directed edge and does not lie on any
edge. If such a triangle does not exist, nothing is returned.
"""
function find_triangle(p::Vertex, D::Delaunay)::Union{Triangle,Nothing}
	for triangle in D.triangles
		ab = triangle.edge
		a, b, c = ab.origin, ab.next.origin, ab.prev.origin
		if is_left(p,a,b) && is_left(p,b,c) && is_left(p,c,a)
			return triangle
		end
	end
	return nothing
end

"""
    flip!(ab::HalfEdge, D::Delaunay):Tuple{Delaunay, Edge, Edge}

Should be called when the following conditions are met:
- AB is in triangle ABC,
- BA is in triangle BAP,
- P is inside circumcircle of ABC.

Then it replaces triangles ABC, BAP with APC, BCP.
Returns updated Delaunay structure and new edges AP,PB.
"""
function flip!(ab::HalfEdge, D::Delaunay)::Tuple{Delaunay, Edge, Edge}
	# Gather edges and triangles
	bc, ca = ab.next, ab.prev
	ba, ap, pb = ab.twin, ab.twin.next, ab.twin.prev
	abc, bap = ab.face, ba.face
	p, c = pb.origin, ca.origin

	# Remove old triangles before constructing new ones
	D = D - abc - bap

	# Construct apc, bcp
	pc, cp = HalfEdges(p,c)
	apc, bcp = Triangle(ap, pc, ca), Triangle(bc, cp, pb)

	return D + apc + bcp, ap, pb
end

"""
	recursive_flip!(ab::Edge, D::Delaunay)::Delaunay

Flips edge AB if it breaks the Delaunay condition. If flipped,
recursively checks AP and PB.
"""
function recursive_flip!(ab::Edge, D::Delaunay)::Delaunay
	if is_delaunay(ab)
		return D
	end

	D, ap, pb = flip!(ab, D)
	D = recursive_flip!(ap, D)
	D = recursive_flip!(pb, D)
	return D
end

function insert_point_no_flip!(p::Vertex, D::Delaunay)::Tuple{Delaunay,Edge,Edge,Edge}
	@assert !isempty(D.triangles) "This should not be empty!"

	abc = find_triangle(p, D)
	@assert abc != nothing "Point $(p) is not inside any triangle!"

	# Gather old objects
	ab, bc, ca = abc.edge, abc.edge.next, abc.edge.prev
	a,b,c = ab.origin, bc.origin, ca.origin

	# Remove old triangle before constructing new ones
	D = D - abc

	# Construct edges and triangles
	(ap, pa), (bp, pb), (cp, pc) = HalfEdges(a,p), HalfEdges(b,p), HalfEdges(c,p)
	abp, bcp, cap = Triangle(ab, bp, pa), Triangle(bc, cp, pb), Triangle(ca, ap, pc)

	return D + abp + bcp + cap, ab, bc, ca
end

function insert_point!(p::Vertex, D::Delaunay)::Delaunay
	D,ab,bc,ca = insert_point_no_flip!(p,D)
	D = recursive_flip!(ab,D)
	D = recursive_flip!(bc,D)
	D = recursive_flip!(ca,D)
	return D
end

#=================================== TESTS ============================================#
using Random

function randomVertex(d::Int)
	return Vertex(round(rand(Float64)*10^d)/(10^d), round(rand(Float64)*10^d)/(10^d))
end

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
#===============================================================================#

end