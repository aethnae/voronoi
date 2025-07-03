export voronoi, sort_vertices_ccw!, polygon_area, areas

"""
	circumcenter(T::Triangle)::Vertex

Calculates circumcenter coordinates of a given triangle.
"""
function circumcenter(T::Triangle)::Vertex
	ab, bc, ca = T.edge, T.edge.next, T.edge.prev
	a,b,c = ab.origin, bc.origin, ca.origin
	D = 2*(a.x * (b.y-c.y) + b.x * (c.y-a.y) + c.x * (a.y - b.y))
    X = ((a.x^2 + a.y^2)*(b.y - c.y) + (b.x^2 + b.y^2)*(c.y - a.y) + (c.x^2 + c.y^2)*(a.y - b.y)) / D
    Y = ((a.x^2 + a.y^2)*(c.x - b.x) + (b.x^2 + b.y^2)*(a.x - c.x) + (c.x^2 + c.y^2)*(b.x - a.x)) / D
    return Vertex(X,Y)
end

"""
	voronoi(D::Delaunay)::Tuple{Dict{Vertex,Vector{Vertex}},Dict{Vertex,Vector{Vertex}}}

Gives the Voronoi to a Delauney
Outputs: a dict V where each center of a Voronoi polygon is mapped to its corners.
         a dict A where each corner of a Voronoi polygon is mapped to its connected corners
"""
function voronoi(D::Delaunay)
    V = Dict{Vertex, Set{Vertex}}() # the centers of Voronoi-polygons and their edges
    A = Dict{Vertex, Set{Vertex}}() # adjaceny list which Voronoi edges are connected

    # collect all inner points
    pts = [e.origin for T in D.triangles for e in (T.edge, T.edge.next, T.edge.prev) if !(e isa Border)]
    pts = unique(pts)

    # only one point -> cell is the whole board
    if length(pts) == 1
        V = Dict(pts[1] => Set{Vertex}())
        A = Dict{Vertex, Set{Vertex}}()
        return V, A
    end

    # two points -> board is devided in the middle
    if length(pts) == 2
        p, q = pts[1], pts[2]
        mid = Vertex((p.x + q.x)/2, (p.y + q.y)/2)
        V = Dict(p => Set([mid]), q => Set([mid]))
        A = Dict(mid => Set([p, q]), p => Set([mid]), q => Set([mid]))
        return V, A
    end

    # Sonderfall: drei Punkte -> gemeinsamer Umkreismittelpunkt
    if length(pts) == 3
        # Finde das Dreieck aus genau diesen drei Punkten
        T0 = findfirst(T -> begin
            vs = (T.edge.origin, T.edge.next.origin, T.edge.prev.origin)
            all(v -> any(u -> u === v, pts), vs)
        end, D.triangles)
        c = circumcenter(T0)
        V = Dict{Vertex, Set{Vertex}}(
            pts[1] => Set([c]),
            pts[2] => Set([c]),
            pts[3] => Set([c])
        )
        A = Dict{Vertex, Set{Vertex}}(
            c      => Set(pts),
            pts[1] => Set([c]),
            pts[2] => Set([c]),
            pts[3] => Set([c])
        )
        return V, A
    end

    # get the centers of every triangle in Delauney
    centers = Dict{Triangle,Vertex}()
    for T in D.triangles
        # nur Dreiecke ohne Border-Kante
        if !any(e -> e isa Border, (T.edge, T.edge.next, T.edge.prev))
            centers[T] = circumcenter(T)
        end
    end
    println(centers)
    # connect the centers
    for T in D.triangles
        c1 = centers[T]
        for e in (T.edge, T.edge.next, T.edge.prev)
            if !(e isa Border)
                he = e::HalfEdge
                if he.twin !== nothing
                    T2 = he.twin.face
                    c2 = centers[T2]
        
                    # in V: he.origin is voronoi-center
                    pts = get!(V, he.origin, Set{Vertex}())
                    push!(pts, c1)
                    push!(pts, c2)
        
                    # in A: c1 <-> c2 in Voronoi-diagram
                    s1 = get!(A, c1, Set{Vertex}())
                    push!(s1, c2)
                    s2 = get!(A, c2, Set{Vertex}())
                    push!(s2, c1)
                end
            end
        end
    end
    return V, A
end 

"""
    sort_vertices_ccw!(verts::Vector{Tuple{Float64,Float64}})

Sorts a list of 2D points in counterclockwise order around their centroid.

Input: vector of tuples (x, y) where x and y are Float64.
Output: vector of tuples (x, y) sorted counterclockwise
"""
function sort_vertices_ccw!(verts::Vector{Tuple{Float64,Float64}})
    # Compute centroid
    cx = mean(x for (x, _) in verts)
    cy = mean(y for (_, y) in verts)

    # Sort by angle around centroid using atan2
    sort!(verts, by = v -> atan2(v[2] - cy, v[1] - cx)) # atan2(y,x) variant of arctan, gives the angle between the point (x,y) and the pos x-axis 

    return sorted_v
end

"""
    polygon_area(verts::Vector{Tuple{Float64,Float64}})

Calculates the area of a Polygon with n corners with shoelace formula

Input: vector of tuples (x, y), where x and y are Float64, sorted counterclockwise.
Output: Area of the Polygon
"""
function polygon_area(verts::Vector{Tuple{Float64,Float64}})
    n = length(verts)
    A = 0.0
    for i in 1:n
        x1, y1 = verts[i]
        x2, y2 = verts[mod1(i+1, n)] # mod1(a,n) gives n if a % n = 0, ensures first point is also last  
        A += (x1 * y2) - (x2 * y1)
    end
    return 0.5 * abs(A)
end

"""
    areas(V::Dict{Vertex, Vector{Vertex}})::Dict{Int, Float64}

Calculates the Areas of the players  

Input: Voronoi dict 
Output: Dict with the area of the polygons belonging to indiviual players
"""
function areas(V::Dict{Vertex, Vector{Vertex}})::Dict{Int, Float64}
    Areas = Dict{Int, Float64}()

    # sort the corners for a Voronoi polygon counterclockwise, calculate its area and add it to the players area
    for center in keys(V)
        V[center] = sort_vertices_ccw!(V[center])
        Areas[center.player] = get!(Areas, center.player, 0) + polygon_area(V[center])
    end
    return Areas
end