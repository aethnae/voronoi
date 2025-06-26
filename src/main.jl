include("logic.jl")

using Random, .Logic, .Logic.DCEL, Debugger

Random.seed!(42);

D = Delaunay(Set{Triangle}());

p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000);

D = insert_point!(p, D)

p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)