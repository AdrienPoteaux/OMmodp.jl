using OMmodp
using Test
using Nemo

include("test-subroutines.jl")

@testset "OMmodp.jl" begin
    # Write your tests here.
    # Testing Taylor expansion
    # Il faut definir le corps de base (F_211)[x][y] ou un truc du genre
    @test TestTaylorExp()
    @test TestAppRoot()
    @test TestPhiExp()
end
