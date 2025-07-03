

export Point2D,
  point2d,
  compute!,
  cross_product,
  center_of_area

#=
##### Point2D
=#

"""
  Point2D

A 2D point with user-defined units and a radius.

# Example
```julia
p = Point2D(1.0, 2.0, 0.5)
```
"""
struct Point2D
    x::UserUnit
    y::UserUnit
    r::Float64 
end

"""
    point2d(x::UserUnit, y::UserUnit, r::Float64 = 1.0)

Create a `Point2D` instance with the given coordinates and radius.
```julia
point2d(UserUnit(1.0), UserUnit(2.0), 0.5)
```
"""
function point2d(x::UserUnit, y::UserUnit, r::Float64 = 1.0)
  return Point2D(x, y, r)
end

Base.convert(::Type{Float64}, u::UserUnit) = Float64(u)

"""
    orientation(p1::Point2D, p2::Point2D, p3::Point2D)
Determine the orientation of the triplet of points (p1, p2, p3). The output is
positive if the points are in counter-clockwise order, negative if they are in
clockwise order.

# Example
```julia
p1 = Point2D(0, 0, 0.01)
p2 = Point2D(1, 0, 0.01)
p3 = Point2D(0, 1, 0.01)
orientation(p1, p2, p3)  !== orientation(p1, p3, p2)  
```
"""
function orientation(p1::Point2D, p2::Point2D, p3::Point2D)
  v1  = Point2D(p2.x - p1.x, p2.y - p1.y, 0)
  v2  = Point2D(p3.x - p1.x, p3.y - p1.y, 0)
  return cross_product(v1, v2)
end

"""
    cross_product(p1::Point2D, p2::Point2D)
Calculate the cross product of two 2D points interpreted as vectors.

# Example
```julia
p1 = Point2D(1, 0, 0.01)
p2 = Point2D(0, 1, 0.01)
cross_product(p1, p2)  # Should return 0.5
```
"""
function cross_product(p1::Point2D, p2::Point2D)
  x1 = convert(Float64, p1.x)
  y1 = convert(Float64, p1.y)
  x2 = convert(Float64, p2.x)
  y2 = convert(Float64, p2.y)

  return (x1 * y2 - x2 * y1) / 2
end

#=

##### Area

=#
"""
    Area2D
A 2D area defined by a vector of `Point2D` points. The area is computed using a
simple triangulation.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
a = Area2D(points)
```
"""
mutable struct Area2D
    v::Vector{Point2D}
    value::Float64
    function Area2D(points::Vector{Point2D})
        if length(points) < 3
            error("An area must have at least 3 points")
        end
        return new(points, -1.0)
    end
end

"""
    area(points::Vector{Point2D})
Create an `Area2D` instance from a vector of `Point2D` points. Sorts the points
in counter-clockwise order starting from the first point.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
a = area(points)
```
"""
function area(points::Vector{Point2D})
  pstart = points[1]
  points_temp = points[2:end]
  sort!(points_temp, lt = (p1, p2) -> orientation(pstart, p1, p2) > 0)
  return Area2D(Point2D[pstart; points_temp])
end

"""
    compute!(area::Area2D)

Compute the area of the `Area2D` instance. If the area has already been computed, it returns the cached value.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
a = area(points)
compute!(a) 
```
"""
function compute!(area::Area2D) 
  area.value >= 0.0 && return area.value

  points = area.v

  ret = 0.0
  for (p1, p2) in zip(points, circshift(points, -1))
    ret += cross_product(p1, p2)
  end
  area.value = ret
  return ret
end

"""
    center_of_area(area::Area2D)

Compute the center of the area defined by the `Area2D` instance. The center is the average of the coordinates of all points in the area.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
a = area(points)
center = center_of_area(a)  # Should return Point2D(0.5, 0.5, 0.01)
```
"""
function center_of_area(area::Area2D)
  # Compute the center of the area
  x_sum = sum(p.x for p in area.v)
  y_sum = sum(p.y for p in area.v)
  n = length(area.v)
  return Point2D(x_sum / n, y_sum / n, 0.01)  # radius is arbitrary
end


#=

##### Graph

=#

"""
    Graph2D
A 2D graph structure containing points, edges, and areas. 

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
edges = [(points[1], points[2]), (points[2], points[3]), (points[3], points[4]), (points[4], points[1])]
areas = [area(points)]
G = Graph2D(points, edges, areas)
```
"""
mutable struct Graph2D
    points::Vector{Point2D}
    edges::Vector{Tuple{Point2D, Point2D}}  
    areas::Vector{Area2D}
    function Graph2D(points::Vector{Point2D} = Point2D[], edges::Vector{Tuple{Point2D, Point2D}} = Tuple{Point2D,Point2D}[], areas::Vector{Area2D} = Area2D[])
        new(points, edges, areas)
    end
end

"""
    which_point(G::Graph2D, x::UserUnit, y::UserUnit)
Find the point in the graph that contains the point (x, y). Returns the point if found, otherwise returns `nothing`.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
G = Graph2D(points)
which_point(G, UserUnit(0.5), UserUnit(0.5))  # Should return Point2D(0, 0, 0.01)
```
"""
function which_point(G::Graph2D, x::UserUnit, y::UserUnit)
  # Find the point that contains the point (x,y)
  for p in G.points
    if abs(p.x - x) < p.r && abs(p.y - y) < p.r
      return p
    end
  end

  return nothing
end

function Base.getindex(G::Graph2D, i::Int)
  return G.points[i]
end

"""
    Base.push!(G::Graph2D, p::Point2D)
Add a point to the graph.
# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01)]
G = Graph2D(points)
push!(G, Point2D(0.5, 0.5, 0.01))  # Adds a new point
```
"""
function Base.push!(G::Graph2D, p::Point2D)
  # Add a point to the graph
  push!(G.points, p)
end

"""
    Base.push!(G::Graph2D, edge::Tuple{Point2D, Point2D})
Add an edge to the graph. The edge must be a tuple of two points that are already in the graph.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01)]
G = Graph2D(points)
push!(G, (Point2D(0, 0, 0.01), Point2D(1, 0, 0.01)))
```
"""
function Base.push!(G::Graph2D, edge::Tuple{Point2D, Point2D})
  # Add an edge to the graph
  if length(edge) == 2 && edge[1] in G.points && edge[2] in G.points
    push!(G.edges, edge)
  else
    error("Edge must be a tuple of two points that are already in the graph")
  end
end

"""
    Base.push!(G::Graph2D, area::Area2D)
Add an area to the graph. The area must be an instance of `Area2D`.

# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01), Point2D(1, 1, 0.01), Point2D(0, 1, 0.01)]
a = area(points)
G = Graph2D(points)
push!(G, a) 
```
"""
function Base.push!(G::Graph2D, area::Area2D)
  push!(G.areas, area)
end

"""
    delete_point(G::Graph2D, p::Point2D)
Delete a point from the graph. This will also remove any edges connected to
this point and any areas that contain this point.
Warning! This operation is not efficient.
# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01)]
G = Graph2D(points)
delete_point(G, Point2D(0, 0, 0.01))  # Removes the point and connected edges
```
"""
function delete_point(G::Graph2D, p::Point2D)
  # Remove a point from the graph
  idx = findfirst(x -> x == p, G.points)
  isnothing(idx) && return
  deleteat!(G.points, idx)
  # Remove edges connected to this point
  filter!(e -> e[1] != p && e[2] != p, G.edges)
  # Remove areas that contain this point
  G.areas = filter(a -> !(p in a.v), G.areas)
end

"""
    delete_all(G::Graph2D)
Delete all points, edges, and areas from the graph. This will clear the graph completely.
# Example
```julia
points = [Point2D(0, 0, 0.01), Point2D(1, 0, 0.01)]
G = Graph2D(points)
delete_all(G)  # Clears the graph
```
"""
function delete_all(G::Graph2D)
  # Clear all points and edges from the graph
  G.points = Point2D[]
  G.edges = Tuple{Point2D, Point2D}[]
  G.areas = Area2D[]
end
