export main

"""
  main()
Create a main window with a drawing area and controls.

# Example
```julia
main() 
```
"""
function main()
	win = GtkWindow("Voronoi")
	g = GtkGrid()
	header = GtkLabel("Spielfeld")
	b = GtkCheckButton("Wichtige Checkbox")

	c = canvas(UserUnit)
	d = GtkScale(:h, 2:10)    
	frame = GtkFrame(c)

	g[1, 1] = header
	g[2, 1] = b
	g[1:2, 2] = d
	g[1:2, 3] = frame
	g.column_homogeneous = true 
	g.column_spacing = 15  

	push!(win, g)

	# Wir erstellen Observer-Objekte, die auf Änderungen reagieren

	graph = Observable(Delaunay())

	return c
end
