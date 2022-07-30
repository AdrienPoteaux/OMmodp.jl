using OMmodp
using Test
using Nemo
using DataStructures

include("test-subroutines.jl")
include("test-valuations.jl")

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
