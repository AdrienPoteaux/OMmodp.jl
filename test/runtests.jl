using OMFacto
using Test
using Nemo
using DataStructures

include("t-subroutines-Kt.jl")
include("t-valuations-Kt.jl")
include("t-valuations-padic.jl")
include("t-valuations-Kt1t2.jl")# in progress

@testset "OMmodp.jl" begin
    # Testing subroutines for GF(211)[[t]][x]
    @test TestTaylorExpKt()
    @test TestAppRootKt()
    @test TestPhiExpKt()
    # Testing valuation / Newton polygon for GF(211)[[t]][x]
    @test TestPhiValKt()
    @test TestPhiNewtonPolygonKt()
    @test TestAllCoeffGivenVKt()# to remove in the end ?
    # Testing valuation / Newton polygon for Qp(3)[x]
    @test TestPhiValPadic()
    @test TestPhiNewtonPolygonPadic()
    @test TestAllCoeffGivenVPadic()# to remove in the end ?
end
