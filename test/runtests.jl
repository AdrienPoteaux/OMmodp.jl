using OMFacto
using Test
using Nemo
using DataStructures

# Functions specific to base rings
include("functions-Kt.jl")
include("functions-Padics.jl")
include("functions-Kt1t2.jl")

# Testing subroutines
include("t-subroutines-Kt.jl")
include("t-valuations-Kt.jl")
include("t-valuations-Padics.jl")
include("t-valuations-Kt1t2.jl")

# Testing main features (irreducibility test, hensel, factorisation algorithm)
include("t-main-Kt1t2.jl")
include("t-main-Kt.jl")

@testset "OMmodp.jl" begin
    # Testing subroutines for GF(211)[[t]][x]
    @test TestTaylorExpKt()
    @test TestAppRootKt()
    @test TestPhiExpKt()
    # Testing valuation / Newton polygon for GF(p)[[t]][x]
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
    @test TestRepresentantKt1t2()
    # Testing irreducibility test over several base rings
    @test TestFirstApproximantsKt1t2()
    @test TestFirstApproximantsKt()
end

## Remark : faudrait tester les Padics avec Q comme corps de base et la valuation p-adique dessus (sans troncation par defaut donc)