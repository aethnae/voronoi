include("datastructure.jl")

#using Colors, Makie, Gtk4, GtkObservables, Graphics
using Gtk4, Graphics

# Spielfeld.
canvas = GtkCanvas()
frame = GtkFrame(canvas, "Spielfeld")
window = GtkWindow(frame, "Canvas")

# Wird automatisch ausgeführt, um das Spielfeld neu zu zeichnen.
@guarded draw(canvas) do widget
    graphics_context = getgc(canvas)
    h = height(canvas)
    w = width(canvas)
    # Paint red rectangle
    rectangle(graphics_context, 0, 0, w, h/2)
    set_source_rgb(graphics_context, 1, 0, 0)
    fill(graphics_context)
    # Paint blue rectangle
    rectangle(graphics_context, 0, 3h/4, w, h/4)
    set_source_rgb(graphics_context, 0, 0, 1)
    fill(graphics_context)
end

#================================== MAUS AKTION ==================================================#

# Observer für Mausaktionen.
gesture = GtkGestureClick()
push!(canvas, gesture)

"""
  on_pressed(controller, n_press, x, y)
Handle mouse button events for adding points, edges, and areas in the drawing area.
The rules are as follows:
- Left click (button 1, no modifiers): Add a point at the clicked position. If a point already exists at that position, it will be deleted.
- Right click (button 3): Start drawing an area. If the first point is clicked again, the area is completed and added to the graph.
- Middle click (button 2): Delete all points, edges, and areas in the graph.
- Left click with Shift modifier (button 1, modifiers == 1): Start drawing a new edge. The first click sets the start point, and the second click sets the end point of the edge.
"""
function on_pressed(controller::GtkGestureClick, n_press, x, y)
  #state = controller.get_current_event_state()

  println("[$(n_press) ($(x), $(y))")
  w=widget(controller)
  graphics_context = getgc(w)
  set_source_rgb(graphics_context, 0, 1, 0)
  arc(graphics_context, x, y, 5, 0, 2pi)
  stroke(graphics_context)
  reveal(w)

  #=
  if btn.button == 1 && btn.modifiers == 0
    add_point!(btn)
  elseif btn.button == 3 
    add_area(btn)
  elseif btn.button == 2
    delete_all!()
  elseif btn.button == 1 && btn.modifiers == 1
    add_edge!(btn)
  end
  =#
end
signal_connect(on_pressed, gesture, "pressed")

"""
  main()
Create a main window with a drawing area and controls.

# Example
```julia
main() 
```
"""
function main()
  # Observer für zeichnen
  graph = Observable(Graph2D())  
  new_edge =  Observable(Point2D[]) 
  drawing_area = Observable(Point2D[]) 

  # Spielfeld
  canvas = GtkCanvas(w=400, h=400)
  frame = GtkFrame(canvas, "Spielfeld")
  gesture = GtkGestureClick()
  signal_connect(on_pressed, gesture, "pressed")  
  push!(canvas,gesture)

  # GUI-Elemente
  label = GtkLabel("Spieleranzahl")
  scale = GtkScale(:h, 2:10)
  button = GtkButton(:mnemonic, "Neues Spiel")

  # GUI-Layout (Spalte, Zeile)
  grid = GtkGrid()
  grid[1, 1] = label
  grid[1, 2] = scale
  grid[2, 1:2] = button
  grid[1:2, 3] = frame
  grid.column_homogeneous = true 
  grid.column_spacing = 15

  # Fenster
  window = GtkWindow(grid, "Voronoi")

  return nothing
end

#============================ MOUSE BUTTON LOGIC VON ÜBUNG =======================================#

#=



MouseButton = Observable{Makie.MouseButtonEvent}


"""
  add_area(btn::MouseButton)
Add an area to the graph by clicking points. If the first point is clicked again, the area is completed and added to the graph.
"""
function add_area(btn::MouseButton)
  p1 = which_point(graph[], btn.position.x, btn.position.y)
  isnothing(p1) && return nothing
  if length(drawing_area[])>0 && p1 == drawing_area.val[1]
    push!(graph[], area(drawing_area[]))
    drawing_area[] = Point2D[]  
    notify(graph) 
  else
    push!(drawing_area[], p1)  
    notify(drawing_area)  
  end
end

"""
  add_point!(btn::MouseButton)
Add a point to the graph at the clicked position. If a point already exists at that position, it will be deleted.
"""
function add_point!(btn::MouseButton)
  point_clicked = which_point(graph[], btn.position.x, btn.position.y)
  
  if !isnothing(point_clicked)
    delete_point(graph[], point_clicked)  # delete the point if it exists
  else
    p = point2d(btn.position.x, btn.position.y)  
    push!(graph[], p)
  end
  notify(graph)  # notify observers of the change
end

"""
  add_edge!(btn::MouseButton)
Add an edge to the graph by clicking two points. The first click sets the start point, and the second click sets the end point of the edge.
"""
function add_edge!(btn::MouseButton)
  if isempty(new_edge[])
    p1 = which_point(graph[], btn.position.x, btn.position.y)
    isnothing(p1) && return nothing

    new_edge.val = [p1]  # start a new line with the first poin
    notify(new_edge)  
    return nothing
  end

  p1 = first(new_edge[])
  p2 = which_point(graph[], btn.position.x, btn.position.y)
  isnothing(p2) && return nothing 

  push!(graph[],(p1, p2))  # add the edge to the graph
  new_edge.val = []  
  notify(graph)  # notify observers
  return nothing 

end

"""
  delete_all!(graph::Graph2D)
Delete all points, edges, and areas in the graph.
"""
function delete_all!()
  delete_all(graph[])  # clear the graph
  notify(graph) 
end

=#

#====================================== DRAWING LOGIC VON ÜBUNG ========================================#
#=
"""
  Drawing Logic
This function is called whenever the drawing area needs to be redrawn.
It draws the background and the graph, including points, edges, and areas.
"""
redraw = draw(c, graph) do cnvs, grs
    #fill!(cnvs, colorant"white")   # Einfacher Hintergrund
    draw_bg(c)

    set_coordinates(cnvs, BoundingBox(0, 100, 0, 100))  # set coordinates to 0..1 along each axis
    ctx = Gtk4.getgc(cnvs)   # gets the "graphics context" object (see Cairo/Gtk)

    draw_graph(ctx,grs)
end

"""
  draw_graph(ctx::CairoContext, g::Graph2D)
Draw the graph on the given Cairo context. It draws areas in green, edges in blue, and points in black.
"""

function draw_graph(ctx, g::Graph2D)
    for a in g.areas
        #set_source_rgb(ctx, 0, 1, 0)  
        draw_area(ctx, a, colorant"green")
    end


    for e in g.edges
        #set_source_rgb(ctx, 0, 0, 1)  
        draw_edge(ctx, e, colorant"blue")  
    end
    for p in g.points
      draw_point(ctx, p.x, p.y, p.r, colorant"black")  
    end

end

"""
  draw_edge(ctx::CairoContext, e::Tuple{Point2D, Point2D}, color)
Draw an edge on the given Cairo context. The edge is defined by a tuple of two points.
"""
function draw_edge(ctx, e::Tuple{Point2D, Point2D}, color)
    start_point, end_point = e
    set_source(ctx, color)
    move_to(ctx, start_point.x, start_point.y)
    line_to(ctx, end_point.x, end_point.y)
    stroke(ctx)
end

"""
  draw_area(ctx::CairoContext, a::Area2D, color)
Draw an area on the given Cairo context. The area is defined by a list of points.
It displays the value of the area in the center of the area.
"""
function draw_area(ctx, a::Area2D, color)
  isempty(a.v) && return
  set_source(ctx, color)
  points = a.v
  p1 = first(points)
  move_to(ctx, p1.x, p1.y)
  for nxt in points[2:end]
    line_to(ctx, nxt.x, nxt.y)
  end
  line_to(ctx, p1.x, p1.y)  # close the area
  fill(ctx)  
  
  center = center_of_area(a)
  x, y = center.x, center.y
  # draw text 
  select_font_face(ctx, "Sans", Cairo.FONT_SLANT_NORMAL,
  Cairo.FONT_WEIGHT_NORMAL);
  set_font_size(ctx, 2.0);
  set_source_rgb(ctx, 0.8,0.8,0.8)
  move_to(ctx, x, y)
  val = floor(compute!(a), digits=2)
  show_text(ctx, string(val))
  fill(ctx)
end

"""
  draw_point(ctx::CairoContext, x, y, r, color)
Draw a point on the given Cairo context. The point is drawn as a circle with the given radius and color.
"""
function draw_point(ctx, x, y, r, color)
    set_source(ctx, color)
    arc(ctx, x, y, r, 0, 2 * π)  # draw a circle at (x,y) with radius r
    fill(ctx)
end

"""
  draw_bg(c::Canvas)
Draw the background of the canvas.
"""
function draw_bg(c)
    ctx = getgc(c)
    h = height(c)
    w = width(c)
    # Paint red rectangle
    rectangle(ctx, 0, 0, w, h)
    set_source_rgb(ctx, 0.8, 0.8, 0.8)
    fill(ctx)
    # Paint blue rectangle
    rectangle(ctx, 0, 3h/4, w, h/4)
    fill(ctx)
end
=#