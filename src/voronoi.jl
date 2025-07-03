module Voronoi

using GtkObservables, Colors
using GtkObservables.Gtk4
using GtkObservables.Cairo

export run_tests

run_tests() = include("../test/runtests.jl")

include("dcel.jl")
include("triangulation.jl")
include("polygons.jl")
include("gui.jl")


end # module Voronoi
