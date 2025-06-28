module Logic

include("dcel.jl")
using .DCEL, LinearAlgebra

export circumcenter, is_delaunay, is_left, find_triangle, flip!, recursive_flip!, insert_point!, insert_point_no_flip!

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


"""
	circumcenter(T::Triangle)::Vertex

Calculates circumcenter coordinates of a given triangle.
"""
function circumcenter(T::Triangle)::Vertex
	D = 2*(a.x * (b.y-c.y) + b.x * (c.y-a.y) + c.x * (a.y - b.y))
    X = ((a.x^2 + a.y^2)*(b.y - c.y) + (b.x^2 + b.y^2)*(c.y - a.y) + (c.x^2 + c.y^2)*(a.y - b.y)) / D
    Y = ((a.x^2 + a.y^2)*(c.x - b.x) + (b.x^2 + b.y^2)*(a.x - c.x) + (c.x^2 + c.y^2)*(b.x - a.x)) / D
    return Vertex(X,Y)
end

end