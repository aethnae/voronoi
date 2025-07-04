using Test, Random

function randomVertex(b::Int)
	return Vertex(round(rand(Float64)*2^b)/(2^b), round(rand(Float64)*2^b)/(2^b))
end

include("triangulation.jl")
include("polygons.jl")