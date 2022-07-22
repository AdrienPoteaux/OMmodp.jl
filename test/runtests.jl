using OMmodp
using Test
using Nemo

include("test-subroutines.jl")
include("test-newtonpolygon.jl")

@testset "OMmodp.jl" begin
    # Testing subroutines
    @test TestTaylorExp()
    @test TestAppRoot()
    @test TestPhiExp()
    # Testing valuation / Newton polygon
    @test TestPhiVal()
end
