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
get_coordinates(n::Node) = n.x

"""
An `Element` is a triangle in the triangulation of the domain

### Fields:
* `vertices` - set of vertices that define an element
"""
struct Element
    vertices::Set{Node}
end
funtion get_barycenter(e::Element)
    barycenter = zeros(2)
    for vertex in e.vertices
        barycenter = barycenter .+ vertex
    end
    barycenter = barycenter ./ 3
    return barycenter
end

"""
A `BasisFunction` is a (linear) function defined over a set of elements

### Fields:
* `interpolants` - set of element-interpolant pairs where (linear) interpolant is defined by coefficients over a given element;
                    the union of these compactly supported planes defines the basis function at an interior node
"""
struct BasisFunction
    interpolants::Dict{Set{Element}, Vector{Float64}}
end
function BasisFunction(support::Set{Element})
    interpolants = ()
    interpolants = Dict(zip(support, interpolants))
    return BasisFunction(interpolants)
end
