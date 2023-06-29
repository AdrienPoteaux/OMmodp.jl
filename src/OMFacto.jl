module OMFacto

# Write your package code here.

using Nemo
using DataStructures

include("subroutines.jl")
include("valuations.jl")

export AppRoot
export TaylorExp
export PhiExp
export PhiVal
export PhiNewtonPolygon
export AllCoeffGivenV # no export in the final version ? (needed for the test to work)
export PhiResidualPol

end