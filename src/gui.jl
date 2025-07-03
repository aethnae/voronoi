export main

using GtkObservables.Gtk4
using Gtk4

"""
  main()

Prompt for number of turns, then create the main window.
"""
function main()
    # --- Prompt dialog for number of turns ---
    dialog = GtkDialog("Voronoi: Anzahl Züge pro Spieler", (), 0, nothing)
    content = dialog.child
    entry = GtkEntry()
    set_gtk_property!(entry, :text, "3")
    push!(content, GtkLabel("Bitte Anzahl Züge pro Spieler eingeben:"))
    push!(content, entry)
    ok_button = GtkButton("OK")
    push!(content, ok_button)
    t = 3  # default

    done = Ref(false)
    signal_connect(ok_button, "clicked") do _
        t_val = tryparse(Int, get_gtk_property(entry, :text, String))
        if t_val !== nothing && t_val > 0
            t = t_val
        end
        done[] = true
        destroy(dialog)
    end

    show(dialog)
    while !done[]
        sleep(0.05)
    end

    # --- Create main window ---
    win = GtkWindow("Voronoi")
    g = GtkGrid()
    c = canvas(UserUnit)
    frame = GtkFrame(c)
    g[1:2, 1] = frame
    g.column_homogeneous = true
    g.column_spacing = 15
    push!(win, g)

    state = Observable(new_game(t))

    function redraw(canvas)
        ctx = getgc(canvas)
        set_source_rgb(ctx, 1, 1, 1)
        paint(ctx)

        # Draw player points
        for (player, points) in state[].player_points
            color = player == 1 ? RGB(0.2,0.4,1.0) : RGB(1.0,0.4,0.2)
            set_source_rgb(ctx, color.r, color.g, color.b)
            for v in points
                arc(ctx, v.x, v.y, 8, 0, 2pi)
                fill(ctx)
            end
        end

        # Draw Voronoi diagram (if enough points)
        try
            V, _ = voronoi(state[].delaunay)
            set_source_rgb(ctx, 0, 0, 0)
            for corners in values(V)
                verts = [(v.x, v.y) for v in corners]
                if !isempty(verts)
                    move_to(ctx, verts[1]...)
                    for v in verts[2:end]
                        line_to(ctx, v...)
                    end
                    close_path(ctx)
                    stroke(ctx)
                end
            end
        catch
        end
    end

    on(state) do _
        redraw(c)
    end

    # --- GTK4 mouse click handling ---
    gesture = GtkGestureClick() # TODO: Ist noch verbuggt, probiert es mal zu zeichnen, dann seht ihr den Error.
    push!(c.widget, gesture)
    signal_connect(gesture, "pressed") do gesture, n_press, x, y
        if place_point!(state[], x, y)
            notify(state)
        end
        return true
    end

    signal_connect(c.widget, "realize") do widget
        redraw(c)
    end

    show(win)
    return win
end
