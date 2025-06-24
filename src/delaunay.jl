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
	is_left(p::Vertex, e::HalfEdge, D::DCELStruct)::Bool

Uses dot product to determine if a point p is left of a
directed edge e. If the point lies exactly on the edge,
it is not considered left of the edge.
"""
function is_left(p::Vertex, e::HalfEdge, D::DCELStruct)::Bool
	a = D.vertices[e.origin]
	b = D.vertices[D.halfedges[e.next].origin]
	n = Vertex(a.y-b.y, b.x-a.x)
	ap = Vertex(p.x-a.x, p.y-a.y)
	return (p.x - a.x) * (a.y - b.y) + (p.y - a.y) * (b.x - a.x) > 0
end


"""
    find_triangle(p::Vertex, D::DCELStruct)::Face

Finds the triangle containing the given vertex. The vertex
is to the left of each directed edge and does not lie on any
edge. If such a triangle does not exist, nothing is returned.
"""
function find_triangle(p::Vertex, D::DCELStruct)::Union{Face,Nothing}
	for face in D.faces
		e1 = D.halfedges[face.outer_component]
		e = D.halfedges[e1.next]
		while e != e1 && is_left(p, e, D)
			e = D.halfedges[e.next]
		end
		if e == e1 && is_left(p, e, D)
			return face
		end
	end
	return nothing
end

function flip!(e::HalfEdge, D::DCELStruct):Nothing
	ab, ba = e, D.halfedges[e.twin]
	abc, bap = D.faces[ab.incident_face]
end

function recursive_flip!(e::HalfEdge, D::DCELStruct)::Nothing
	
end

function insert_point!(p::Vertex, D::DCELStruct)::Nothing
end
