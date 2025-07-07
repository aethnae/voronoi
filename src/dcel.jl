# DCEL data structures for Voronoi diagram representation

import Base.+, Base.-, Base.:*, Base.==, Base.round

export Vertex, Edge, Border, HalfEdge, HalfEdges, Triangle, Delaunay

#========================================================================#
abstract type Edge end      # Directed edge which can be a bound
abstract type Face end
#========================================================================#
struct Vertex   # Represents a point in the plane, placed by a player.
    x::Float64
    y::Float64
    player::Union{Int, Nothing}

    Vertex(X::Float64, Y::Float64) = new(X,Y,nothing)
    Vertex(X::Float64, Y::Float64, P::Int) = new(X,Y,P)
end

BottomLeft = Vertex(0.0,-2.0)
BottomRight = Vertex(3.0,0.0)
TopLeft = Vertex(0.0,3.0)
OuterVertices = (TopLeft, BottomLeft, BottomRight)
is_inside(v::Vertex) = v.x >= 0.0 && v.y >= 0.0 && v.x <= 1.0 && v.y <= 1.0


Base.show(io::IO, v::Vertex) = print(io, "($(v.x), $(v.y))")

a::Vertex + b::Vertex = Vertex(a.x + b.x, a.y + b.y)
a::Vertex - b::Vertex = Vertex(a.x - b.x, a.y - b.y)
a::Vertex * b::Vertex = a.x*b.x + a.y*b.y
k::Float64 * b::Vertex = Vertex(k * b.x, k * b.y)
round(a::Vertex;digits=10,base=2) = Vertex(round(a.x,digits=digits,base=base), round(a.y,digits=digits,base=base),a.player)

#========================================================================#
mutable struct Border <: Edge
    origin::Vertex              # Origin vertex
    next::Union{Edge,Nothing}   # Next edge in the face
    prev::Union{Edge,Nothing}   # Previous edge in the face
    face::Union{Face,Nothing}   # Face to which this edge belongs

    Border(o::Vertex) = new(o, nothing, nothing, nothing)
    Border(o,n,p,f) = new(o,n,p,f)
end

Base.show(io::IO, b::Border) = print(io, "Border[$(b.origin)->$(b.next.origin)]")

mutable struct HalfEdge <: Edge
    origin::Vertex                  # Origin vertex
    twin::Union{HalfEdge,Nothing}   # Twin half-edge
    next::Union{Edge,Nothing}       # Next edge in the face
    prev::Union{Edge,Nothing}       # Previous edge in the face
    face::Union{Face,Nothing}       # Face to which this half-edge belongs

    HalfEdge(o::Vertex) = new(o, nothing, nothing, nothing, nothing)
    HalfEdge(o,t,n,p,f) = new(o,t,n,p,f)
end

function HalfEdges(a::Vertex,b::Vertex)::Tuple{HalfEdge, HalfEdge}
    ab, ba = HalfEdge(a), HalfEdge(b)
    ab.twin, ba.twin = ba, ab
    return ab, ba
end

function Base.show(io::IO, e::HalfEdge)
    v = e.next !== nothing ? e.next.origin : nothing
    print(io, "[$(e.origin)->$(v)]")
end

e1::Edge == e2::Edge = e1.origin == e2.origin && e1.next.origin == e2.next.origin

#========================================================================#
mutable struct Triangle <: Face
    edge::Union{Edge,Nothing}       # One of the edges of the face, can be a half-edge

    function Triangle(xy::Edge, yz::Edge, zx::Edge)
        xyz = new(xy)
        xy.next, yz.next, zx.next = yz, zx, xy
        xy.prev, yz.prev, zx.prev = zx, xy, yz
        xy.face = yz.face = zx.face = xyz
        return xyz
    end
end

T1::Triangle == T2::Triangle = T1.edge in (T2.edge, T2.edge.next, T2.edge.prev)

Base.show(io::IO, T::Triangle) = print(io, "Tri{$(T.edge), $(T.edge.next), $(T.edge.prev)}")

#========================================================================#
mutable struct Delaunay
    triangles::Set{Triangle}  # Set of triangles in the Delaunay triangulation

    Delaunay() = new(Set([Triangle(
        Border(TopLeft),
        Border(BottomLeft),
        Border(BottomRight)
    )]))
end

function Base.show(io::IO, D::Delaunay)
    for T in D.triangles
        println(io, T)
    end
end

D::Delaunay + T::Triangle = begin
    D.triangles = push!(D.triangles, T)
    return D
end

D::Delaunay - T::Triangle = begin
    D.triangles = delete!(D.triangles, T)
    return D
end
