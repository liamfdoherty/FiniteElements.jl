using FiniteElements

# Define a grid
h = 1/3
grid = Grid(h)

# Define the set of elements on the grid
elements = generate_elements(grid)

# Compute the basis functions at the interior nodes
basis_functions = []
for point in grid.interior_points
    # compute the basis function at that point
    support = Vector{Element}()
    for element in elements
        if point in element.vertices
            push!(support, element)
        end
    end
    push!(basis_functions, BasisFunction(Set(support)))
end
