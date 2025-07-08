# DCEL data structures for Voronoi diagram representation

import Base.+, Base.-, Base.:*, Base.==, Base.round

export Vertex, Edge, Border, HalfEdge, HalfEdges, Triangle, Delaunay

#========================================================================#
"""
    abstract type Edge

Represents an edge in the Delaunay triangulation.

# Subtypes
- @Border Edge outside of polygon space, used for initialization.
- @HalfEdge Edge created when adding new vertices.
"""
abstract type Edge end      # Directed edge which can be a bound
abstract type Face end
#========================================================================#

"""
    struct Vertex

Represents a point in the plane, placed by a player.

# Attributes
- `x::Float64` The x-coordinate of the vertex, must be in [0,1].
- `y::Float64` The y-coordinate of the vertex, must be in [0,1].
- `player::Union{Int, Nothing}` The player number, must be in {1,2}. For testing purposes,
a `nothing` value is permitted.

# Constructors
- `Vertex(x,y,p)` Constructs a vertex with `player=p`.
- `Vertex(x,y)` Constructs a vertex with `player=nothing`.

# Functions
- `a::Vertex + b::Vertex` adds the coordinates into a vertex with `player=nothing`.
- `a::Vertex - b::Vertex` subtracts the coordinates into a vertex with `player=nothing`.
- `a::Vertex * b::Vertex` calculates the dot product, x*x + y*y.
- `k::Float64 * b::Vertex` scales the coordinates by a factor.
"""
struct Vertex
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

"""
    is_inside(v::Vertex)

Determines if the vertex coordinates place it inside [0,1]².
"""
is_inside(v::Vertex) = v.x >= 0.0 && v.y >= 0.0 && v.x <= 1.0 && v.y <= 1.0


Base.show(io::IO, v::Vertex) = print(io, "($(v.x), $(v.y))")

a::Vertex + b::Vertex = Vertex(a.x + b.x, a.y + b.y)
a::Vertex - b::Vertex = Vertex(a.x - b.x, a.y - b.y)
a::Vertex * b::Vertex = a.x*b.x + a.y*b.y
k::Float64 * b::Vertex = Vertex(k * b.x, k * b.y)

"""
    round(a::Vertex; digits=10, base=2)

Rounds the x- and y-coordinates of a vertex to the given precision using @round.
"""
round(a::Vertex;digits=10,base=2) = Vertex(round(a.x,digits=digits,base=base), round(a.y,digits=digits,base=base),a.player)

#========================================================================#
"""
    mutable struct Border <: Edge

An outer border edge used for triangulation initialization.
Edges enclose their triangle in anticlockwise order. Due to the
double-linked references, attributes may be initialized as `nothing`.

# Attributes
- `origin::Vertex` The origin vertex of this directed edge.
- `next::Union{Edge,Nothing}` The next edge of the triangle.
- `prev::Union{Edge,Nothing}` The preceding edge of the triangle.
- `face::Union{Face,Nothing}` The triangle face containing this edge.

# Constructors
- `Border(o::Vertex)` only sets the origin vertex. Other attributes are `nothing`.
- `Border(o,n,p,f)` sets all attributes.
"""
mutable struct Border <: Edge
    origin::Vertex              # Origin vertex
    next::Union{Edge,Nothing}   # Next edge in the face
    prev::Union{Edge,Nothing}   # Previous edge in the face
    face::Union{Face,Nothing}   # Face to which this edge belongs

    Border(o::Vertex) = new(o, nothing, nothing, nothing)
    Border(o,n,p,f) = new(o,n,p,f)
end

Base.show(io::IO, b::Border) = print(io, "Border[$(b.origin)->$(b.next.origin)]")

"""
    mutable struct HalfEdge <: Edge

An edge created when adding new points to the triangulation. Each
half-edge is only contained within a single triangle. Adjacent triangles
contain a `twin` of this edge pointing in the opposite direction.
Edges enclose their triangle in anticlockwise order. Due to the
double-linked references, attributes may be initialized as `nothing`.

# Attributes
- `origin::Vertex` The origin vertex of this directed edge.
- `twin::Union{HalfEdge,Nothing}` The twin half-edge.
- `next::Union{Edge,Nothing}` The next edge of the triangle.
- `prev::Union{Edge,Nothing}` The preceding edge of the triangle.
- `face::Union{Face,Nothing}` The triangle face containing this edge.

# Constructors
- `HalfEdge(o::Vertex)` only sets the origin. Other attributes are `nothing`.
- `HalfEdge(o,t,n,p,f)` sets all attributes.
"""
mutable struct HalfEdge <: Edge
    origin::Vertex                
    twin::Union{HalfEdge,Nothing} 
    next::Union{Edge,Nothing}     
    prev::Union{Edge,Nothing}     
    face::Union{Face,Nothing}     

    HalfEdge(o::Vertex) = new(o, nothing, nothing, nothing, nothing)
    HalfEdge(o,t,n,p,f) = new(o,t,n,p,f)
end

"""
    HalfEdges(a::Vertex,b::Vertex)::Tuple{HalfEdge, HalfEdge}

Using vertices `a` and `b`, constructs two twin half-edges AB and BA.
"""
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
"""
    mutable struct Triangle <: Face

A triangle of the triangulation. Convenience object to connect edges.

# Attributes
- `edge::Union{Edge,Nothing}` One of the edges of the triangle.

# Constructors
- `Triangle(xy::Edge, yz::Edge, zx::Edge)` constructs a triangle XYZ from given
edges XY, YZ and ZX. This connects their `next`, `prev` and `face` attributes.

# Functions
- `T1::Triangle == T2::Triangle` determines if the triangles share any edge.
"""
mutable struct Triangle <: Face
    edge::Union{Edge,Nothing}

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
"""
    mutable struct Delaunay

The main triangulation object.

# Attributes
- `triangles::Set{Triangle}` Set of triangles in the Delaunay triangulation.

# Constructors
- `Delaunay()` Constructs a pseudo-empty triangulation, containing only one large
triangle encompassing the entire polygon.

# Functions
- `D::Delaunay + T::Triangle` adds the triangle to the triangulation.
- `D::Delaunay - T::Triangle` removes the triangle from the triangulation.
"""
mutable struct Delaunay
    triangles::Set{Triangle} 

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
