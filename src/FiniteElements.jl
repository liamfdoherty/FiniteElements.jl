module FiniteElements
using GridInterpolations
using LinearAlgebra

# Elements
include("Elements/elements.jl")

# Grids
include("Grids/grids.jl")

# Quadrature
include("Quadrature/quadrature.jl")

# Export list
include("exports.jl")
end
