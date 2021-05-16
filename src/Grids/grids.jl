"""
A `Grid` is a set of nodes that discretize a domain

### Fields
* `interior` - set of grid points on the interior of the domain
* `boundary` - set of grid points on the boundary of the domain
"""
struct Grid
    h::Float64
    points::Vector{Node}
    interior_points::Vector{Node}
end
function Grid(h::Float64)
    grid = RectangleGrid(0:h:1, 0:h:1)
    points = [Node(point) for point in grid]
    interior_points = [Node(point) for point in grid if (point[1] != 0) && (point[1] != 1) && (point[2] != 0) && (point[2] != 1)]
    return Grid(h, points, interior_points)
end

"""
`generate_elements` - generates a set of elements that tile a given (rectangular) grid using a step size h, using triangles with positive slope

### Fields
* `grid` - rectangular grid object
* `h` - step size defining grid
"""
function generate_elements(grid::Grid)
    M = Int(1/grid.h)
    elements = []
    for (index, node) in enumerate(grid.points)
        #= An upper element is the upper triangle in a grid square, lower elements are the lower triangles in the grid square,
            where the defining node of a square is the bottom left one =#
        if (index > M * (M + 1)) || (mod(index, M + 1) == 0) # check if node is a valid bottom left endpoint of a square; if not, skip node
            continue
        end
        up_node = grid.points[index + M + 1]
        right_node = grid.points[index + 1]
        diagonal_node = grid.points[index + M + 2]

        push!(elements, Element([node, up_node, diagonal_node])) # Push the upper element for the square
        push!(elements, Element([node, right_node, diagonal_node])) # Push the lower element for the square
    end
    return elements
end
