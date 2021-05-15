include("../Elements/elements.jl")
"""
A `Grid` is a set of nodes that discretize a domain

### Fields
* `interior` - set of grid points on the interior of the domain
* `boundary` - set of grid points on the boundary of the domain
"""
struct Grid
    interior::Vector{Node}
    boundary::Vector{Node}
end
function Grid(h::Float64)
    grid = RectangleGrid(0:h:1, 0:h:1)
    interior = []; boundary = []
    for point in grid
        if (point[1] != 0.) && (point[1] != 1.) && (point[2] != 0.) && (point[2] != 1.)
            push!(interior, NTuple{2, Float64}(point))
        else
            push!(boundary, NTuple{2, Float64}(point))
        end
    end
    return Grid(interior, boundary)
end
