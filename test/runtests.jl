using Test, Random

function randomVertex(d::Int)
	return Vertex(round(rand(Float64)*10^d)/(10^d), round(rand(Float64)*10^d)/(10^d))
end

include("triangulation.jl")
include("polygons.jl")