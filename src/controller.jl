
#==================== Game State ====================#

mutable struct GameState
    t::Int                                      # Number of points per player
    current_player::Int                         # 1 or 2
    turns_left::Dict{Int,Int}                   # player => turns left
    delaunay::Delaunay                          # Current Delaunay triangulation
    player_points::Dict{Int, Vector{Vertex}}    # player => list of their points
end

function new_game(t::Int)
    GameState(
        t,
        1,
        Dict(1 => t, 2 => t),
        Delaunay(),
        Dict(1 => Vector{Vertex}(), 2 => Vector{Vertex}())
    )
end

#==================== Game Logic ====================#

"""
    place_point!(state::GameState, x::Float64, y::Float64)

Attempts to place a point for the current player at (x, y).
Returns true if successful, false otherwise.
"""
function place_point!(state::GameState, x::Float64, y::Float64)
    if state.turns_left[state.current_player] == 0
        return false
    end
    v = Vertex(x, y, state.current_player)
    insert_point!(v, state.delaunay)
    push!(state.player_points[state.current_player], v)
    state.turns_left[state.current_player] -= 1
    # Switch player if both have turns left
    if state.turns_left[3 - state.current_player] > 0
        state.current_player = 3 - state.current_player
    end
    return true
end

"""
    game_over(state::GameState)::Bool

Returns true if both players have placed all their points.
"""
game_over(state::GameState) = all(v == 0 for v in values(state.turns_left))

"""
    compute_areas(state::GameState)

Returns a Dict mapping player numbers to their Voronoi area.
"""
function compute_areas(state::GameState)
    V, _ = voronoi(state.delaunay)
    return areas(V)
end

"""
    winner(state::GameState)

Returns the player number with the largest area, or nothing if tied.
"""
function winner(state::GameState)
    area_dict = compute_areas(state)
    if isempty(area_dict)
        return nothing
    end
    max_area = maximum(values(area_dict))
    winners = [p for (p, a) in area_dict if a == max_area]
    return length(winners) == 1 ? winners[1] : nothing
end
