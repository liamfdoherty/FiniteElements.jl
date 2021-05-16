"""
`barycentric_quadrature` - computes the integral of a function over an element using the barycentric quadrature rule

### Fields:
* `f` - function to compute the integral of
* `e` - element that defines the domain of integration
"""
function barycentric_quadrature(f::Function, e::Element)
    area = get_area(element)
    barycenter = get_barycenter(element)
    return area * f(barycenter[1], barycenter[2])
end
