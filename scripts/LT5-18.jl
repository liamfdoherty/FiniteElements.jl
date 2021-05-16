using FiniteElements, LinearAlgebra

function finite_element_2d(f::Function, h::Float64)
    # Define a grid
    M = Int(1/h)
    grid = Grid(h)

    # Define the set of elements on the grid
    elements = generate_elements(grid)

    # Compute the basis functions at the interior nodes
    basis_functions = []
    for point in grid.interior_points
        support = Vector{Element}()
        for element in elements
            if point in element.vertices
                push!(support, element)
            end
        end
        push!(basis_functions, BasisFunction(Set(support)))
    end

    # Assemble stiffness matrix
    A = zeros((M - 1)^2, (M - 1)^2)
    for (i, Φᵢ) in enumerate(basis_functions), (j, Φⱼ) in enumerate(basis_functions)
        A[i, j] = compute_bilinear_form(Φᵢ, Φⱼ)
    end

    # Assemble load vector
    b = zeros((M - 1)^2)
    for (index, basis_function) in enumerate(basis_functions)
        b[index] = compute_linear_functional(f, basis_function)
    end

    # Solve system
    U = A\b

    return U
end

# Define the problem and true solution
function f(x, y)
    return sin(π*x)sin(π*y) + sin(π*x)sin(2*π*y)
end

function TrueSolution(x, y)
   return 1/(2*π^2) * sin(π*x)sin(π*y) +  1/(5*π^2) * sin(π*x)sin(2*π*y)
end


# Obtain solution
h_vals = [1/10, 1/20, 1/50]
for h in h_vals
    solution = finite_element_2d(f, h)
    plt = plot(h:h:1-h, h:h:1-h, solution, st=:surface, title = "Numerical Solution; h = $(h)")
    display(plt)
end
