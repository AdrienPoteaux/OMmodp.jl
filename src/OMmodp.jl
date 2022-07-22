module OMmodp

# Write your package code here.

using Nemo

include("subroutines.jl")
include("newtonpolygon.jl")

export AppRoot
export TaylorExp
export PhiExp
export PhiVal

end