# DCEL data structures for Voronoi diagram representation

module DCEL

export Vertex, HalfEdge, Face, Dreieck, Delaunay

abstract type Face end

struct Vertex   # Represents a point in the plane, placed by a player.
    x::Float64
    y::Float64
end

mutable struct HalfEdge
    origin::Vertex          # Origin vertex of the half-edge
    twin::HalfEdge          # Twin half-edge
    next::HalfEdge          # Next half-edge in the face
    prev::HalfEdge          # Previous half-edge in the face
    face::Face              # Face to which this half-edge belongs
end

mutable struct Dreieck <: Face
    halfedge::HalfEdge      # One of the half-edges of the face
end

struct Delaunay
    triangles::Set{Dreieck}  # Set of triangles in the Delaunay triangulation
end

end
