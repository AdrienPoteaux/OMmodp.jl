module OMFacto

# Write your package code here.

using Nemo
using DataStructures

include("subroutines.jl")
include("valuations.jl")
include("irreducibility.jl")
include("hensel.jl")

function (F::fqPolyRepField)(n::fpFieldElem) return F(lift(n)) end # conversion from Fp to Fq. Do not know how to add that from "base ring file" and OMFacto.()() for instance

export AppRoot
export TaylorExp
export PhiExp
export PhiExpMonomials
export PhiMonomialsEval
export PhiVal
export PhiNewtonPolygon
export AllCoeffGivenV # no export in the final version ? (needed for the test to work)
export PhiResidualPol
export FirstApproximants
export Representant
export PhiHensel

end