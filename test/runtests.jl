using Test, Random

Tests = 128

"""
	randomVertex(; rounded=true)::Vertex

Builds a `Vertex` with randomly generated x and y coordinates.
If `rounded==true`, rounds to 10 digit decimal precision.
"""
function randomVertex(; rounded=true)::Vertex
	V = Vertex(rand(Float64),rand(Float64))
	return rounded ? round(V) : V
end

include("triangulation.jl")
include("polygons.jl")