"""
A `Grid` is a set of nodes that discretize a domain

### Fields
* `interior` - set of grid points on the interior of the domain
* `boundary` - set of grid points on the boundary of the domain
"""
struct Grid
    points::Vector{Node}
end
function Grid(h::Float64)
    grid = RectangleGrid(0:h:1, 0:h:1)
    points = [Node(point) for point in grid]
    return Grid(points)
end
