"""
`Node`s are spatial points, defined by a tuple of coordinates.

### Fields:
* `x` - tuple of coordinates for node
"""
struct Node
    x::NTuple{2, Float64}
end
function Node(x::Vector{Float64})
    return Node(NTuple{2, Float64}(x))
end

"""
`get_coordinates` - get the coordinates of an input node

### Fields
`n` - input node
"""
get_coordinates(n::Node) = n.x

"""
An `Element` is a triangle in the triangulation of the domain

### Fields:
* `vertices` - set of vertices that define an element
"""
struct Element
    vertices::Set{Node}
end
function Element(nodes::Vector{Node})
    return Element(Set(nodes))
end

"""
`get_barycenter` - compute the barycenter of an element (for the barycentric quadrature rule)

### Fields:
* `e` - input element
"""
function get_barycenter(e::Element)
    barycenter = zeros(2)
    for vertex in e.vertices
        barycenter = barycenter .+ get_coordinates(vertex)
    end
    barycenter = barycenter ./ 3
    return barycenter
end

"""
`get_area` - compute the area of an element

### Fields:
* `e` - element to compute the area of
"""
function get_area(e::Element)
    A, B, C = e.vertices
    A_coords = get_coordinates(A)
    B_coords = get_coordinates(B)
    C_coords = get_coordinates(C)
    return abs(A_coords[1]*(B_coords[2] - C_coords[2]) + B_coords[1]*(C_coords[2] - A_coords[2]) + C_coords[1]*(A_coords[2] - B_coords[2]))/2
end

"""
A `BasisFunction` is a (linear) function defined over a set of elements

### Fields:
* `interpolants` - set of element-interpolant pairs where (linear) interpolant is defined by coefficients over a given element;
                    the union of these compactly supported planes defines the basis function at an interior node
"""
struct BasisFunction
    interpolants::Dict{Element, Vector{Float64}}
end
function BasisFunction(support::Set{Element})
    # Check to make sure that the correct number of elements are present (not general, but a check for this problem)
    if length(support) != 6
        throw(ErrorException("Invalid support!"))
    end

    # Find all nodes that define the region of support for the basis function, and find the common node
    all_nodes = Set([])
    for element in support
        union!(all_nodes, element.vertices)
    end
    remaining_nodes = all_nodes
    for element in support
        intersect!(remaining_nodes, element.vertices)
    end
    common_node = remaining_nodes

    # Check to make sure that there is only one common node
    if length(common_node) != 1
        throw(ErrorException("Invalid support!"))
    end

    # Initialize set of interpolants and function values at the nodes
    interpolants = []
    for element in support
        plane_points = []
        apex = []
        for node in element.vertices
            point = [coordinate for coordinate in get_coordinates(node)]
            if node in common_node # Check if the node is the common node; if so, assign it a function value of 1
                push!(point, 1.)
                apex = point
            else # Otherwise, assign it a function value of 0
                push!(point, 0.)
            end
            push!(plane_points, point)
        end

        # Compute the linear interpolant over the element
        vecA = [x - y for (x,y) in zip(plane_points[1], plane_points[2])]
        vecB = [x - y for (x,y) in zip(plane_points[2], plane_points[3])]
        n = cross(vecA, vecB)
        interpolant = Vector{Float64}(undef, 3)
        interpolant[1] = -n[1]/n[3]
        interpolant[2] = -n[2]/n[3]
        interpolant[3] = n[1]*apex[1]/n[3] + n[2]*apex[2]/n[3] + apex[3]
        push!(interpolants, interpolant)
    end

    # Package the interpolants with their corresponding elements
    interpolants = Dict(zip(support, interpolants))
    return BasisFunction(interpolants)
end
