include("diagram.jl")
using .DCEL

#=
Neuer Punkt p
	1. Finde Dreieck ABC, in dem p liegt
	2. Ersetze ABC durch ABP, BCP, CAP
	3. Prüfe, ob drei Umkreise keine weiteren Ecken enthalten
	4. Edge-Flip rekursiv

Schritt 3: Umkreis-Test
=#

"""
	is_left(p::Vertex, e::HalfEdge, D::Delaunay)::Bool

Uses dot product to determine if a point p is left of a
directed edge e. If the point lies exactly on the edge,
it is not considered left of the edge.
"""
function is_left(p::Vertex, e::HalfEdge, D::Delaunay)::Bool
	a,b = e.origin, e.next.origin
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
		e1 = triangle.halfedge
		e = e1.next
		while e != e1 && is_left(p, e, D)
			e = e.next
		end
		if e == e1 && is_left(p, e, D)
			return triangle
		end
	end
	return nothing
end

function flip!(e::HalfEdge, D::Delaunay):Tuple{Delaunay, HalfEdge, HalfEdge}
	abc, bap = e.face, e.twin.face
	bc, ca = e.next, e.next.next
	ba, ap, pb = e.twin, e.twin.next, e.twin.next.next

	# Construct apc, bcp
	pc, cp, apc, bcp = HalfEdge(p, cp, ca, ap, apc), HalfEdge(c, pc, pb, bc, bcp),
		Dreieck(pc), Dreieck(cp)

	# Replace
	newTriangles = delete!(D.triangles, abc)
	newTriangles = delete!(newTriangles, bap)
	newTriangles = push!(newTriangles, apc)
	newTriangles = push!(newTriangles, bcp)
	D.triangles = newTriangles

	return D, pc, cp
end

function recursive_flip!(e::HalfEdge, D::Delaunay)::Nothing
	
end

function insert_point!(p::Vertex, D::Delaunay)::Nothing
end