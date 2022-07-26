module OMmodp

# Write your package code here.

using Nemo
using DataStructures

include("subroutines.jl")
include("newtonpolygon.jl")

export AppRoot
export TaylorExp
export PhiExp
export PhiVal
export PhiNewtonPolygon

end