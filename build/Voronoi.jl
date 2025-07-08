module Voronoi

using GtkObservables, Colors
using GtkObservables.Gtk4
using GtkObservables.Cairo
using LinearAlgebra

export GameState, new_game, place_point!, game_over, compute_areas, winner
export Delaunay, Vertex			# dcel.jl
export insert_point!			# triangulation.jl
export voronoi, areas 			# polygons.jl
export main						# gui.jl
export run_tests				# runtests.jl

run_tests() = include("test/runtests.jl")

include("dcel.jl")
include("triangulation.jl")
include("polygons.jl")
include("gui.jl")
include("controller.jl")

end # module Voronoi
