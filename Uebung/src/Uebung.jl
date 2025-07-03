module Uebung

using GtkObservables, Colors
using GtkObservables.Gtk4
using GtkObservables.Cairo

export run_tests


greet() = print("Hello World!")
run_tests() = include("../test/runtests.jl")

include("datastructure.jl")
include("visual.jl")

end # module Uebung
