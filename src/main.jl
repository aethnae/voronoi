include("logic.jl")

using Random, .Logic, .Logic.DCEL, Debugger

function checkConnected(D)::Bool
	for T in D.triangles
		e1, e2, e3 = T.edge, T.edge.next, T.edge.next.next
		@assert e2.prev === e1
		@assert e3.prev === e2
		@assert e3.next === e1
	end
	return true
end



Random.seed!(42);

D = Delaunay();

p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000);

D = insert_point!(p, D)

p = Vertex(round(rand(Float64)*1000)/1000, round(rand(Float64)*1000)/1000)