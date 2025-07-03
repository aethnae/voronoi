module Voronoi

using GtkObservables, Colors
using GtkObservables.Gtk4
using GtkObservables.Cairo

export run_tests
export GameState, new_game, place_point!, game_over, compute_areas, winner
export Delaunay, Vertex, insert_point!, voronoi, areas

run_tests() = include("../test/runtests.jl")

include("dcel.jl")
include("triangulation.jl")
include("polygons.jl")
include("gui.jl")
include("controller.jl")

end # module Voronoi
