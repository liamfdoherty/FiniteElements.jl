"""
`compute_bilinear_form` - computes entries of stiffness matrix

### Fields:
* `Φ₁` - first basis function
* `Φ₂` - second basis function
"""
function compute_bilinear_form(Φ₁::BasisFunction, Φ₂::BasisFunction)
    support1 = collect(keys(Φ₁.interpolants)); support2 = collect(keys(Φ₂.interpolants))
    overlapping_support = intersect(support1, support2)
    if length(overlapping_support) == 0
        return 0.
    end
    sum = 0.
    for element in overlapping_support
        area = get_area(element)
        sum += dot(Φ₁.interpolants[element][1:2], Φ₂.interpolants[element][1:2]) * area
    end
    return sum
end

"""
`compute_linear_functional` - computes entries of load vector

### Fields:
* `f` - right hand side function
* `Φ` - basis function
"""
function compute_linear_functional(f::Function, Φ::BasisFunction)
    elements = collect(keys(Φ.interpolants))
    sum = 0.
    for element in elements
        barycenter = get_barycenter(element)
        area = get_area(element)
        coefficients = Φ.interpolants[element]
        sum += area * f(barycenter[1], barycenter[2]) * (coefficients[1]*barycenter[1] + coefficients[2]*barycenter[2] + coefficients[3])
    end
    return sum
end

"""
`compute_error` - cmopute the L² norm of the error

### Fields
* `TrueSolution` - the true solution of the PDE
* `NumericalSolution` - the numerical solution of the PDE using finite elements
"""
function compute_error(TrueSolution::Function, NumericalSolution::Vector{Float64}, h::Float64)
    grid = Grid(h)
    elements = generate_elements(grid)

    error = 0.

    for element in elements
        vertices = element.vertices
        barycenter = get_barycenter(element)
        area = get_area(element)

        value_at_barycenter = 0.
        for vertex in vertices
            index = findall(x->x==vertex, grid.points)[1]
            value_at_vertex = NumericalSolution[index]
            value_at_barycenter += value_at_vertex
        end
        value_at_barycenter = value_at_barycenter/3

        true_barycenter_value = TrueSolution(barycenter[1], barycenter[2])
        error += (true_barycenter_value - value_at_barycenter)^2 * area
    end
    error = sqrt(error)
    return error
end
