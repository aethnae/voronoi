include("diagram.jl")
using .DCEL, LinearAlgebra

#=
Neuer Punkt p
	1. Finde Dreieck ABC, in dem p liegt
	2. Ersetze ABC durch ABP, BCP, CAP
	3. Prüfe, ob drei Umkreise keine weiteren Ecken enthalten
	4. Edge-Flip rekursiv

Schritt 3: Umkreis-Test
=#

"""
	check_circumference(e::HalfEdge)::Bool

Checks if the circumcircle of the triangle formed by the half-edge `e` and its adjacent half-edge contains any other vertices.
Returns True if the circumcircle contains no other vertices.
"""
function check_circumference(e::HalfEdge)::Bool
	a = e.origin
	b = e.next.origin
	c = e.prev.origin

	adj = e.twin
	d = adj.prev.origin

	return det([a.x a.y a.x^2+a.y^2 1; b.x b.y b.x^2+b.y^2 1; c.x c.y c.x^2+c.y^2 1; d.x d.y d.x^2+d.y^2 1]) > 0
end

"""
	is_left(p::Vertex, e::HalfEdge)::Bool

Uses dot product to determine if a point p is left of a
directed edge e. If the point lies exactly on the edge,
it is not considered left of the edge.
"""
function is_left(p::Vertex, a::Vertex, b::Vertex)::Bool
	return (p.x - a.x) * (a.y - b.y) + (p.y - a.y) * (b.x - a.x) > 0
end


"""
    find_triangle(p::Vertex, D::Delaunay)::Face

Finds the triangle containing the given vertex. The vertex
is to the left of each directed edge and does not lie on any
edge. If such a triangle does not exist, nothing is returned.
"""
function find_triangle(p::Vertex, D::Delaunay)::Union{Face,Nothing}
	for triangle in D.triangles
		ab = triangle.halfedge
		a, b, c = ab.origin, ab.next.origin, ab.prev.origin
		if is_left(p,a,b) && is_left(p,b,c) && is_left(p,c,a)
			return triangle
		end
	end
	return nothing
end

"""
    flip!(ab::HalfEdge, D::Delaunay):Tuple{Delaunay, HalfEdge, HalfEdge}

Should be called when the following conditions are met:
- AB is in triangle ABC,
- BA is in triangle BAP,
- P is inside circumcircle of ABC.

Then it replaces triangles ABC, BAP with APC, BCP.
Returns updated Delaunay structure and new edges AP,PB.
"""
function flip!(ab::HalfEdge, D::Delaunay):Tuple{Delaunay, HalfEdge, HalfEdge}
	# Gather edges and triangles
	bc, ca = ab.next, ab.prev
	ba, ap, pb = ab.twin, ab.twin.next, ab.twin.prev
	abc, bap = ab.face, ba.face
	p, c = pb.origin, ca.origin

	# Construct apc, bcp
	pc, cp, apc, bcp = HalfEdge(p, cp, ca, ap, apc), HalfEdge(c, pc, pb, bc, bcp),
		Dreieck(pc), Dreieck(cp)

	# Replace
	D.triangles = delete!(D.triangles, abc)
	D.triangles = delete!(D.triangles, bap)
	D.triangles = push!(D.triangles, apc)
	D.triangles = push!(D.triangles, bcp)

	return D, ap, pb
end

"""
	recursive_flip!(ab::HalfEdge, D::Delaunay)::Delaunay

Checks if edge AB in triangle ABC should be flipped.
If yes: calls `flip!`, then recursively checks AP and PB.
"""
function recursive_flip!(ab::HalfEdge, D::Delaunay)::Delaunay
	if check_circumference(ab)
		return D
	end
	D, ap, pb = flip!(ab, D)
	D = recursive_flip!(ap, D)
	D = recursive_flip!(pb, D)

	return D
end

function insert_point!(p::Vertex, D::Delaunay)::Delaunay
	if isempty(D.triangles)
		# Construct large triangle
		
	end

	abc = find_triangle(p, D)
	@assert abc != nothing "Point $(p) is not inside any triangle!"

	# Gather old objects
	ab, bc, ca = abc.halfedge, abc.halfedge.next, abc.halfedge.prev
	a,b,c = ab.origin, bc.origin, ca.origin

	# Construct edges AP, BP, CP, their twins, and ABP, BCP, CAP
	ap, bp, cp, pa, pb, pc, abp, bcp, cap =
		HalfEdge(a, pa, pc, ca, cap),
		HalfEdge(b, pb, pa, ab, abp),
		HalfEdge(c, pc, pb, bc, bcp),
		HalfEdge(p, ap, ab, bp, abp),
		HalfEdge(p, bp, bc, cp, bcp),
		HalfEdge(p, cp, ca, ap, cap),
		Triangle(ab), Triangle(bc), Triangle(ca)

	D.triangles = delete!(D.triangles, abc)
	D.triangles = push!(D.triangles, abp)
	D.triangles = push!(D.triangles, bcp)
	D.triangles = push!(D.triangles, cap)

	D = recursive_flip!(ab,D)
	D = recursive_flip!(bc,D)
	D = recursive_flip!(ca,D)

	return D
end