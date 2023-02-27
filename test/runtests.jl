#using OMmodp # si je comprends bien il ne faut plus l'import car il est deja lu dans le main car def ici ?
using Test
using Nemo
using DataStructures

include("t-subroutines-Kt.jl")
include("t-valuations-Fpt.jl")

@testset "OMmodp.jl" begin
    # Testing subroutines
    @test TestTaylorExp()
    @test TestAppRoot()
    @test TestPhiExp()
    # Testing valuation / Newton polygon
    @test TestPhiVal()
    @test TestPhiNewtonPolygon()
    @test TestAllCoeffGivenV()# to remove in the end ?
end
