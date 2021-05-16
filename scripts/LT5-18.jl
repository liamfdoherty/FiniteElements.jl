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
    for i in 1:length(basis_functions), j in 1:length(basis_functions)
        Φᵢ = basis_functions[i]; Φⱼ = basis_functions[j]
        support_i = collect(keys(Φᵢ.interpolants)); support_j = collect(keys(Φⱼ.interpolants))
        overlapping_support = intersect(support_i, support_j)
        if length(overlapping_support) == 0
            continue
        else
            # Compute the sum of the bilinear form over each element in the common support
            sum = 0
            for element in overlapping_support
                sum += dot(Φᵢ.interpolants[element][1:2], Φⱼ.interpolants[element][1:2]) * get_area(element)
            end
            A[i, j] = sum
        end
    end

    # Assemble load vector
    b = zeros((M - 1)^2)


    return A
end
