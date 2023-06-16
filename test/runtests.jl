using STOPIT
using Test

@testset "STOPIT.jl" begin

    @testset "Absorber Tests" begin
        include("absorbertests.jl")
    end

    @testset "Desorb Tests" begin
        include("desorbtests.jl")
    end

end