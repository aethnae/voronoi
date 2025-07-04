export main

"""
  main()

Prompt for number of turns, then create the main window.
"""
function main()
    """
        Prompt dialog for number of turns
    """
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

    """
        Window and canvas setup

    `canvas(UserUnit)` lets GtkObservables automatically convert pixel position to [0,1] x [0,1].
    """
    win = GtkWindow("Voronoi")
    g = GtkGrid()
    c = canvas(UserUnit, 400,400)
    frame = GtkFrame(c)
    push!(win,frame)

    state = Observable(new_game(t))

    """
        Canvas drawing

    `redraw` is a Gtk4 inbuilt function which is called whenever the window resizes.
    Additionally, GtkObservables allows us to hook additional observables which should
    also trigger a function call.
    """
    redraw = draw(c, state) do canvas, st
        ctx = getgc(canvas)
        set_coordinates(ctx, BoundingBox(0,1,0,1))
        set_source_rgb(ctx, 1, 1, 1)
        paint(ctx)

        # Draw Voronoi diagram (if enough points)
        try
            V, _ = voronoi(st.delaunay)
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

        # Draw player points
        for (player, points) in st.player_points
            color = player == 1 ? RGB(0.2,0.4,1.0) : RGB(1.0,0.4,0.2)
            set_source_rgb(ctx, color.r, color.g, color.b)
            for v in points
                arc(ctx, v.x, v.y, 0.01, 0, 2pi)
                fill(ctx)
            end
        end
    end

    """
        Mouse click handling
    """
    on(c.mouse.buttonpress) do btn
        if place_point!(state[], btn.position.x.val, btn.position.y.val)
            notify(state)
        end
    end

    return win
end
