
export voronoi, sort_vertices_ccw!, polygon_area, areas, intersect_ray_bbox, filter_internal_triangles, distance, rotate_right

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
    filter_internal_triangles(D::Delaunay)::Set{Triangle}

Filter out the triangles that have one corner in common with the initial Triangle in Delauney
"""
function filter_internal_triangles(D::Delaunay)::Set{Triangle}
    inner_triangles = Set{Triangle}()
    S1 = Vertex(-10.0, -10.0) # initial triangle
    S2 = Vertex(20.0, -10.0)
    S3 = Vertex(0.0, 20.0)
    S = [S1, S2, S3]
    for T in D.triangles
        e = T.edge
        e_next = T.edge.next
        e_prev = T.edge.prev
        if !(e.origin in S || e_next.origin in S ||e_prev.origin in S)
            push!(inner_triangles, T)
        end 
    end
    return inner_triangles
end

"""
    intersect_ray_bbox(origin::Vertex, dir::Tuple{Float64,Float64}, bbox::Vector{Vertex})

Calculates the intersection of a ray starting in bounding box(bbox)in direction of dir with the edges of the box  
"""
function intersect_ray_bbox(origin::Vertex, dir::Tuple{Float64,Float64}, bbox::Vector{Vertex})
    xs = [p.x for p in bbox]
    ys = [p.y for p in bbox]
    xmin = minimum(xs)
    xmax = maximum(xs)
    ymin = minimum(ys)
    ymax = maximum(ys)

    ox, oy = origin.x, origin.y
    println("ox: $(ox), oy: $(oy)")
    dx, dy = dir
    ts = Float64[]
    # calculate t for every edge of bbox, where ox+dx*t = x_edge oder oy+dy*t = y_edge
    if dx != 0
        push!(ts, (xmin - ox)/dx)
        push!(ts, (xmax - ox)/dx)
    end
    if dy != 0
        push!(ts, (ymin - oy)/dy)
        push!(ts, (ymax - oy)/dy)
    end
    # only positive t
    ts = filter(t -> t > 0, ts)
    tmin = minimum(ts)
    return Vertex(ox + dx*tmin, oy + dy*tmin)
end

"""
    rotate_right(h::HalfEdge)::Tuple{Float64,Float64}

Takes a HalfEdge and returns a direction, that is HalfEdge rotated by 90° counterclockwise
"""
function rotate_right(h::HalfEdge)::Tuple{Float64,Float64}
    a = h.origin
    b = h.next.origin
    v = (b.x-a.x, b.y-a.y)
    rv = (v[2], -v[1])
    return rv
end

"""
    distance(a::Vertex, b::Vertex)::Float64

calculates the distance between two vertices
"""
function distance(a::Vertex, b::Vertex)::Float64
    return sqrt(abs2(b.x-a.x)+ abs2(b.y-a.y))
end

"""
	voronoi(D::Delaunay)::Tuple{Dict{Vertex,Vector{Vertex}},Dict{Vertex,Vector{Vertex}}}

Gives the Voronoi to a Delaunay
Outputs: a dict V where each center of a Voronoi polygon is mapped to its corners.
		 a dict A where each corner of a Voronoipolygon is mapped to its connected corners
"""
function voronoi(D::Delaunay, bbox::Vector{Vertex})::Dict{Vertex, Vector{Vertex}}
    V = Dict{Vertex, Vector{Vertex}}() # the centers of Voronoi-polygons and their edges
	
    ld = bbox[1]
    rd = bbox[2]
    ru = bbox[3]
    lu = bbox[4]

    S1 = Vertex(-10.0, -10.0)
    S2 = Vertex(20.0, -10.0)
    S3 = Vertex(0.0, 20.0)

    # collect all inner points
    pts = [e.origin for T in D.triangles for e in (T.edge, T.edge.next, T.edge.prev)]
    pts = unique(pts)
    pts = collect(setdiff(pts, [S1, S2, S3])) 
    #println("Inner points: $(pts)")

    # only one point -> cell is the whole board
    if length(pts) == 1
        V = Dict(pts[1] => [ld, rd, ru, lu]) 
        return V
    end 

    # two points -> board is devided in the middle
    if length(pts) == 2
        p, q = pts[1], pts[2]
        mid = Vertex((p.x + q.x)/2, (p.y + q.y)/2)
        v = (q.x-p.x, q.y-p.y) # direction p to q
        rv = (v[2], -v[1])
        lv = (-v[2],v[1])
        far1 = intersect_ray_bbox(mid, rv, bbox)
        far2 = intersect_ray_bbox(mid,lv, bbox)
        cellp = Vector{Vertex}()
        cellq = Vector{Vertex}()
        for corner in bbox
            corner_distance_p = distance(corner,p)
            corner_distance_q = distance(corner,q)
            if corner_distance_p < corner_distance_q  # corner belongs to p cell
                push!(cellp, corner)
            elseif corner_distance_p > corner_distance_q  # corner belongs to q cell
                push!(cellq, corner)
            else
                push!(cellp, corner)  #corner belongs to both cells (is= far1 or far2)
                push!(cellq, corner)
            end
        end

        # far1 and far2 are in both cells
        push!(cellp, far1)
        push!(cellp, far2)
        push!(cellq, far1)
        push!(cellq, far2)

        cellp = sort_vertices_ccw!(unique(cellp))
        cellq = sort_vertices_ccw!(unique(cellq))

        V = Dict(p => cellp, q => cellq)
        return V
    end

    # get the centers of every triangle 
    centers = Dict{Triangle,Vertex}()
    inner_triangles = filter_internal_triangles(D) # only triangles without border edge
    #println("Inner triagnles: $(inner_triangles)")
    for T in inner_triangles
        centers[T] = circumcenter(T)
    end
    #println("triangle centers: $(centers)")

    # inserects of rays to the board corners
    intersects = Vector{Vertex}()
    for T in inner_triangles
        for e in (T.edge,T.edge.next,T.edge.prev)
            he = e::HalfEdge
            # twin face not internal or not existent
            if !(he.twin.face in inner_triangles) #||he.twin === nothing 
                c1 = centers[T]
                dir = rotate_right(he) # direction to the board borders
                #println("normal direction to $(he): $(dir)")
                far = intersect_ray_bbox(c1, dir, bbox)
                #println("intersect to direction $(dir): $(far)")
                push!(intersects, far)
                #s = get!(V, p, Set{Vertex}()); push!(s, far)
            end
        end
    end

    # border points to cell centers
    # board corners 
    for corner in bbox
        min_point, min_dist = pts[1], 2*maximum([c.x for c in bbox]) #bigger than max dist inside board
        for pt in pts
            #println("Min-dist: $(min_dist)")
            #println("Corner points dist: $(distance(corner, pt))")
            if distance(corner, pt) < min_dist
                min_point = pt
                min_dist = distance(corner, pt)
            end
        end
        #println("closest fpoint to $(corner): $(min_point)")
        s = get!(V, min_point, Vector{Vertex}()); push!(s,corner) # corner belongs to the inner point its closest to
    end
    println("V after board corners: $(V)")

    # intersects: each intersect is a cell corner of two inner points, the ones it's closest to
    for far in intersects
        point_distances = Vector{Tuple{Float64,Vertex}}()
        for pt in pts
            push!(point_distances, (distance(far, pt),pt))
        end
        point_distances = sort!(point_distances, by = x -> x[1])
        min_point1 = point_distances[1][2]
        min_point2 = point_distances[2][2]
        #println("points closest to $(far): $(min_point1) and $(min_point2)")
        s = get!(V, min_point1, Vector{Vertex}()); push!(s,far)
        s = get!(V, min_point2, Vector{Vertex}()); push!(s,far)
    end
    #println("V after intersects: $(V)")

    # connect the centers
    for T in inner_triangles
        c1 = centers[T]
        for e in (T.edge, T.edge.next, T.edge.prev)
            he = e::HalfEdge
            if length(inner_triangles)==1
                s = get!(V, he.origin, Vector{Vertex}())
                push!(s, c1)
            end
            if he.twin !== nothing && he.twin.face in inner_triangles 
                T2 = he.twin.face
                c2 = centers[T2]
    
                # in V: he.origin is voronoi-center
                s = get!(V, he.origin, Vector{Vertex}())
                push!(s, c1)
                push!(s, c2)
            end
        end
    end

    for pt in keys(V)
        V[pt] = sort_vertices_ccw!(unique(V[pt])) # sort the corners of each Voronoi polygon counterclockwise
    end
    return V
end 

"""
	sort_vertices_ccw!(verts::Vector{Tuple{Float64,Float64}})

Sorts a list of 2D points in counterclockwise order around their centroid.

Input: vector of tuples (x, y) where x and y are Float64
Output: vector of tuples (x, y) sorted counterclockwise
"""
function sort_vertices_ccw!(verts::Vector{Vertex})
    n = length(verts)
    if n == 0
        return verts  # nothing to sort
    end
    # Compute centroid
    cx = sum(p.x for p in verts) / n
    cy = sum(p.y for p in verts) / n

    # Sort by angle around centroid 
    sort!(verts, by = v -> atan(v.y - cy, v.x - cx)) # atan(y,x) variant of arctan, gives the angle between the point (x,y) and the pos x-axis 

    return verts
end

"""
	polygon_area(verts::Vector{Tuple{Float64,Float64}})::Float64

Calculates the area of a Polygon with n corners with shoelace formula

Input: vector of tuples (x, y), where x and y are Float64, sorted counterclockwise.
Output: Area of the Polygon
"""
function polygon_area(verts::Vector{Vertex})::Float64
    n = length(verts)
    A = 0.0
    for i in 1:n
        x1, y1 = verts[i].x, verts[i].y
        #println("x1: $(x1), y1: $(y1)")
        x2, y2 = verts[mod1(i+1, n)].x, verts[mod1(i+1, n)].y # mod1(a,n) gives n if a % n = 0, ensures first point is also last  
        A += (x1 * y2) - (x2 * y1)
        #println("Area in step $(i): $(A)")
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
	# calculate its area and add it to the players area
	for center in keys(V)
		Areas[center.player] = get!(Areas, center.player, 0) + polygon_area(V[center])
	end
	return Areas
end

