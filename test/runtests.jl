using Test, Random

function randomVertex(b::Int)
	return round(Vertex(rand(Float64),rand(Float64)))
end

include("triangulation.jl")
include("polygons.jl")