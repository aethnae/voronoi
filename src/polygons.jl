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
function intersect_ray_bbox(origin::Vertex, dir::Tuple{Float64,Float64}, bbox::Vector{Vertex})::Union{Vertex, Nothing}
    ox, oy = origin.x, origin.y
    dx, dy = dir

    xs = [p.x for p in bbox]
    ys = [p.y for p in bbox]
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    ts = Float64[]

    # Nur t-Werte sammeln, bei denen der Strahl die Box schneiden KANN
    if dx != 0
        tx1 = (xmin - ox) / dx
        tx2 = (xmax - ox) / dx
        if tx1 > 0
            y_at_tx1 = oy + dy * tx1
            if ymin ≤ y_at_tx1 ≤ ymax
                push!(ts, tx1)
            end
        end
        if tx2 > 0
            y_at_tx2 = oy + dy * tx2
            if ymin ≤ y_at_tx2 ≤ ymax
                push!(ts, tx2)
            end
        end
    end

    if dy != 0
        ty1 = (ymin - oy) / dy
        ty2 = (ymax - oy) / dy
        if ty1 > 0
            x_at_ty1 = ox + dx * ty1
            if xmin ≤ x_at_ty1 ≤ xmax
                push!(ts, ty1)
            end
        end
        if ty2 > 0
            x_at_ty2 = ox + dx * ty2
            if xmin ≤ x_at_ty2 ≤ xmax
                push!(ts, ty2)
            end
        end
    end

    if isempty(ts)
        return nothing  # Kein gültiger Schnittpunkt
    end

    tmin = minimum(ts)
    return Vertex(ox + dx * tmin, oy + dy * tmin)
end

"""
    is_in_box(p::Vertex, bbox::Vector{Vertex})::Bool

Checks if the given Vertex is inside or on the border ob bbox
"""
function is_in_box(p::Vertex, bbox::Vector{Vertex})::Bool
    xs = [p.x for p in bbox]
    ys = [p.y for p in bbox]
    xmin = minimum(xs)
    xmax = maximum(xs)
    ymin = minimum(ys)
    ymax = maximum(ys)

    return xmin <= p.x <= xmax && ymin <= p.y <= ymax
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
    in_triangle(p::Vertex, T::Triangle)::Bool

checks if the point p is a corner of the triangle T
"""
function in_triangle(p::Vertex, T::Triangle)::Bool
    return p==T.edge.origin || p==T.edge.next.origin || p==T.edge.prev.origin
end

"""
    get_points(T::Triangle)::Vector{Vertex}

returns the corners of a triangle in a vetctor
"""
function get_points(T::Triangle)::Vector{Vertex}
    a = T.edge.origin
    b = T.edge.next.origin
    c = T.edge.prev.origin
    return [a, b, c]
end

"""
    get_common_points(T1::Triangle, T2::Triangle)::Union{Vector{Vertex},Nothing}

calculates the common corners of two triangles. Returns nothing if there are none
"""
function get_common_points(T1::Triangle, T2::Triangle)::Union{Vector{Vertex},Nothing}
    points1 = get_points(T1)
    points2 = get_points(T2)
    common_points = [p for p in points1 if p in points2]
    if isempty(common_points)
        return nothing
    else
        return common_points
    end
end

"""
	voronoi(D::Delaunay)::Tuple{Dict{Vertex,Vector{Vertex}},Dict{Vertex,Vector{Vertex}}}

Gives the Voronoi to a Delaunay
Outputs: a dict V where each center of a Voronoi polygon is mapped to its corners.
		 a dict A where each corner of a Voronoipolygon is mapped to its connected corners
"""
function voronoi(D::Delaunay, bbox::Vector{Vertex})::Dict{Vertex, Vector{Vertex}}
    V = Dict{Vertex, Vector{Vertex}}() # the centers of Voronoi-polygons and their edges
	inner_triangles = Set{Triangle}()
    centers = Dict{Triangle,Vertex}()

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
        #println("V: $(V)")
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
        #println("V: $(V)")
        return V
    end

    
    inner_triangles = filter_internal_triangles(D) # only triangles without border edge
    # get the centers of every triangles
    #println("All triangles: $(D.triangles)")
    #println("Inner triangles: $(inner_triangles)")
    for T in D.triangles
        centers[T] = circumcenter(T)
    end
    #println("triangle centers: $(centers)")

    # find the triangles that share an edge
    connected_Tri = Dict{Triangle,Vector{Triangle}}()
    for T in D.triangles
        for e in (T.edge, T.edge.next, T.edge.prev)
            he = e
            #println("Current edge: $(e) in $(T)")
            if !(he isa Border) && he.twin !== nothing && he.twin.face in D.triangles 
                T2 = he.twin.face
                k = get!(connected_Tri, T, Vector{Triangle}()); push!(k, T2)
            end
        end
    end

    #=
    for p in inner_triangles
        println("Center to $(p): $(centers[p])")
    end =#

    
    # for 3 points add edge to cells
    if length(inner_triangles)==1 
        for T in inner_triangles
            c = centers[T]
            if is_in_box(c, bbox) 
                for p in (T.edge.origin, T.edge.next.origin, T.edge.prev.origin)
                    # in V: p is cell center with edge c
                    s = get!(V, p, Vector{Vertex}())
                    push!(s, c)
                end
            end 
        end     
    end

    # inserects of connecting lines between centers
    for T in D.triangles
        c1 = centers[T]
        for T2 in connected_Tri[T]
            c2 = centers[T2]

            # case: c1 out and c2 in
            if !is_in_box(c1,bbox) && is_in_box(c2,bbox)
                dir = (c1.x-c2.x, c1.y-c2.y)
                far = intersect_ray_bbox(c2, dir, bbox)
                common_points = get_common_points(T, T2)
                for p in common_points
                    if p in pts
                        push!(get!(V, p, Vector{Vertex}()), far)
                    end
                end
            elseif is_in_box(c1,bbox) && !is_in_box(c2,bbox) # case: c2 out and c1 in
                dir = (c2.x-c1.x, c2.y-c1.y)
                far = intersect_ray_bbox(c1, dir, bbox)
                common_points = get_common_points(T, T2)
                for p in common_points
                    if p in pts
                        push!(get!(V, p, Vector{Vertex}()), far)
                    end
                end
            elseif !is_in_box(c1,bbox) && !is_in_box(c2,bbox) # case: both out
                dir1 = (c2.x-c1.x, c2.y-c1.y)
                dir2 = (-dir1[1],-dir1[2])
                far1 = intersect_ray_bbox(c1, dir1, bbox)
                far2 = intersect_ray_bbox(c2, dir2, bbox)
                common_points = get_common_points(T, T2)
                for p in common_points
                    if p in pts && far1 !== nothing
                        push!(get!(V, p, Vector{Vertex}()), far1)
                    end
                    if p in pts && far2 !== nothing
                        push!(get!(V, p, Vector{Vertex}()), far2)
                    end
                end
            else # both are in
                points1 = get_points(T)
                points2 = get_points(T2)
                all_points = [points1; points2]
                for p in unique(all_points)
                    if p in pts && in_triangle(p,T)
                        push!(get!(V, p, Vector{Vertex}()), c1)
                    elseif p in pts && in_triangle(p,T2)     
                        push!(get!(V, p, Vector{Vertex}()), c2)
                    end
                end
            end
        end
    end

    #=for pt in keys(V)
        V[pt] = sort_vertices_ccw!(unique(V[pt])) # sort the corners of each Voronoi polygon counterclockwise
    end
    println("V after connected centers: $(V)")=#

    for T in inner_triangles
        c1 = centers[T] 
        for e in (T.edge, T.edge.next, T.edge.prev)
            he = e::HalfEdge
            #twin face not internal 
            if !(he.twin.face in inner_triangles) #|| he.twin === nothing
                dir = rotate_right(he) # direction to the board borders
                dir_ = (-dir[1],-dir[2])
                #println("normal direction to $(he): $(dir)")
                if is_in_box(c1, bbox) 
                    far = intersect_ray_bbox(c1, dir, bbox)
                    #println("intersect to edge $(he): $(far)")
                    push!(get!(V, he.origin, Vector{Vertex}()), far)
                    push!(get!(V, he.next.origin, Vector{Vertex}()), far)
                else # lines perpendicular to edges left to c1 have two intersects with borders
                    if is_left(c1, he.origin, he.next.origin)
                        mid_he = Vertex((he.origin.x + he.next.origin.x)/2, (he.origin.y + he.next.origin.y)/2)
                        far1 = intersect_ray_bbox(mid_he, dir, bbox)
                        far2 = intersect_ray_bbox(mid_he, dir_, bbox)
                        #println("intersects to edge $(he): $(far1) and $(far2)")
                        if  far1 !== nothing
                            push!(get!(V, he.origin, Vector{Vertex}()), far1)
                            push!(get!(V, he.next.origin, Vector{Vertex}()), far1)
                        end
                        if far2 !== nothing
                            push!(get!(V, he.origin, Vector{Vertex}()), far2)
                            push!(get!(V, he.next.origin, Vector{Vertex}()), far2)
                        end
                    end
                    #println("center NOT in box: Intersect for $(he) is $(far)")
                end
                #println("intersect to direction $(dir): $(far)")
            end
        end
    end

    #=for pt in keys(V)
        V[pt] = sort_vertices_ccw!(unique(V[pt])) # sort the corners of each Voronoi polygon counterclockwise
    end
    println("V after intersects Inner Tri: $(V)") =#

    # border points to cell centers
    # board corners 
    for corner in bbox
        min_point, min_dist = pts[1], 2*maximum([c.x for c in bbox]) #bigger than max dist inside board
        for pt in pts
            if distance(corner, pt) < min_dist
                min_point = pt
                min_dist = distance(corner, pt)
            end
        end
        #println("closest fpoint to $(corner): $(min_point)")
        s = get!(V, min_point, Vector{Vertex}()); push!(s,corner) # corner belongs to the inner point its closest to
    end

    for pt in keys(V)
        V[pt] = sort_vertices_ccw!(unique(V[pt])) # sort the corners of each Voronoi polygon counterclockwise
    end
    #println("V after board corners: $(V)")
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


#============================================================================================
============================================================================================#
walls = (Vertex(0.,0.), Vertex(1.,0.)),
        (Vertex(1.,0.), Vertex(0.,1.)),
        (Vertex(1.,1.), Vertex(-1.,0.)),
        (Vertex(0.,1.), Vertex(0.,-1.))

function intersect_wall_between(v1::Vertex, v2::Vertex)::Union{Vertex,Nothing}
    intersection = nothing
    foreach(walls) do wall
        # solve v1 + (v2-v1)t = wall
        Mat = [v2.x-v1.x   -wall[2].x;
               v2.y-v1.y   -wall[2].y]
        d = det(Mat)

        # check if points are not parallel to wall
        if isapprox(d, 0.0; atol=1e-3)
            return
        end

        # check if intersection is between the two points
        T = inv(Mat) * [wall[1].x-v1.x;   wall[1].y-v1.y]
        if T[1] <= 0.0 || T[1] >= 1.0
            return
        end

        i = round(v1 + T[1]*(v2-v1))
        if is_inside(i) && (intersection isa Nothing || T[1] < intersection[2])
            intersection = (i, T[1])
        end
    end 
    return intersection isa Nothing ? nothing : intersection[1]
end

function is_in_polygon(V::Vertex, P::Vector{Vertex})::Bool
    for i in 1:length(P)
        if !is_left(V, P[i], P[i % length(P) + 1])
            return false
        end
    end
    return true
end

function voronoi_2(D::Delaunay)::Dict{Vertex, Vector{Vertex}}
    Polygons = Dict{Vertex,Vector{Vertex}}()

    # Collect unbounded polygon vertices.
    foreach(D.triangles) do T
        C = circumcenter(T)
        foreach((T.edge.origin, T.edge.next.origin, T.edge.prev.origin)) do v
            if is_inside(v)
                vertices = haskey(Polygons, v) ? Polygons[v] : Vector{Vertex}()
                push!(vertices, C)
                Polygons[v] = vertices
            end
        end
    end

    # Sort polygon vertices anticlockwise.
    foreach(keys(Polygons)) do v
        polygon = sort([x for x in Set(Polygons[v])], by = (p -> atan(p.y-v.y, p.x-v.x)))
        Polygons[v] = polygon
    end

    # Calculate bounded polygon vertices.
    foreach(keys(Polygons)) do v
        polygon = Polygons[v]
        newpoly = Vector{Vertex}()
        N = length(polygon)

        # Discard outside vertices and add up to two wall intersections from consecutive vertices.
        for i in 1:N
            wall = intersect_wall_between(polygon[(i-2+N) % N + 1], polygon[i])
            if wall != nothing
                push!(newpoly, wall)
                wall = intersect_wall_between(wall, polygon[i])
                if wall != nothing
                    push!(newpoly, wall)
                end
            end
            if is_inside(polygon[i])
                push!(newpoly, polygon[i])
            end
        end

        # Add all bounding box corners inside the original polygon.
        foreach((Vertex(0.,0.),Vertex(0.,1.),Vertex(1.,0.),Vertex(1.,1.))) do corner
            if is_in_polygon(corner, polygon)
                push!(newpoly, corner)
            end
        end

        # Sort polygon vertices anticlockwise.
        Polygons[v] = sort(newpoly, by = (p -> atan(p.y-v.y, p.x-v.x)))
    end

    return Polygons
end

function areas_2(V::Dict{Vertex, Vector{Vertex}})::Dict{Int, Float64}
    Areas = Dict{Int, Float64}()
    foreach(keys(V)) do v
        p = v.player isa Nothing ? 0 : v.player
        Areas[p] = (haskey(Areas, p) ? Areas[p] : 0.) + polygon_area(V[v])
    end
    return Areas
end
