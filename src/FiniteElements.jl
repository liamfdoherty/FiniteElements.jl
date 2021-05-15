module FiniteElements
using StaticArrays
using GridInterpolations

# Grids
include("Grids/grids.jl")

# Elements
include("Elements/elements.jl")

# Quadrature
include("Quadrature/quadrature.jl")

# Export list
include("exports.jl")
end
