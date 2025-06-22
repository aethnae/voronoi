# DCEL data structures for Voronoi diagram representation

module DCEL

export Vertex, HalfEdge, Face, DCELStruct

struct Vertex   # Represents a point in the plane, placed by a player.
    x::Float64
    y::Float64
end

struct HalfEdge
    origin::Int           # Index of vertice from which half-edge emanates.
    twin::Int             # Index from halfedges array (the dual edge in the opposite direction).
    next::Int             # Next halfedge of traversal sequence (given by it's index).
    prev::Int             # Previous halfedge in the traversal sequence (given by it's index).
    incident_face::Int    # Index of face that borders half-edge to the left (when traversing).
end

struct Face
    site::Int             # Index of generator point placed by a player.
    outer_component::Int  # Index into halfedges array (one of the bounding edges).
end

struct DCELStruct
    vertices::Vector{Vertex}
    halfedges::Vector{HalfEdge}
    faces::Vector{Face}
end

end
