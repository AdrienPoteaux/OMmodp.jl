using OMFacto
using Test
using Nemo
using DataStructures

include("t-subroutines-Kt.jl")
include("t-valuations-Kt.jl")
include("t-valuations-padic.jl")
include("t-valuations-Kt1t2.jl")

@testset "OMmodp.jl" begin
    # Testing subroutines for GF(211)[[t]][x]
    @test TestTaylorExpKt()
    @test TestAppRootKt()
    @test TestPhiExpKt()
    # Testing valuation / Newton polygon for GF(211)[[t]][x]
    @test TestPhiValKt()
    @test TestPhiNewtonPolygonKt()
    @test TestAllCoeffGivenVKt()# to remove in the end ?
    @test TestPhiResidualPolKt1t2()
    # Testing valuation / Newton polygon for Qp(3)[x]
    @test TestPhiValPadic()
    @test TestPhiNewtonPolygonPadic()
    @test TestAllCoeffGivenVPadic()# to remove in the end ?
    @test TestPhiResidualPolPadic()
    # Testing valuation / Newton polygon for GF(149)[[t1,t2]][x]
    @test TestPhiValKt1t2()
    @test TestPhiNewtonPolygonKt1t2()
    @test TestAllCoeffGivenVKt1t2()# to remove in the end ?
    @test TestPhiResidualPolKt1t2()
end
