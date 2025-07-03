export main

using GtkObservables.Gtk4

"""
  main()

Prompt for number of turns, then create the main window.
"""
function main()
    # Prompt dialog for number of turns
    dialog = GtkDialog("Voronoi: Anzahl Züge pro Spieler", 0, nothing)
    content = dialog.child
    entry = GtkEntry()
    set_gtk_property!(entry, :text, "3")
    push!(content, GtkLabel("Bitte Anzahl Züge pro Spieler eingeben:"))
    push!(content, entry)
    response = nothing

    signal_connect(dialog, "response") do widget, resp
        response = resp
        destroy(dialog)
    end

    dialog.add_button("OK", 1)
    showall(dialog)
    wait(dialog)

    t = tryparse(Int, get_gtk_property(entry, :text, String))
    if t === nothing || t < 1
        t = 3
    end

    # Create main window
    win = GtkWindow("Voronoi")
    g = GtkGrid()
    c = canvas(UserUnit)
    frame = GtkFrame(c)
    g[1:2, 1] = frame
    g.column_homogeneous = true
    g.column_spacing = 15
    push!(win, g)

    state = Observable(new_game(t))

    showall(win)
    return win
end
