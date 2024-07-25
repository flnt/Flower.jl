using Flower
using Test

@testset "Example tests" begin

    @testset "Math tests" begin
        include("validation_test.jl")
    end

    # @testset "Zero velocity" begin
    #     include("zerovel.jl")
    # end

    # @testset "Constant velocity" begin
    #     include("constantvel.jl")
    # end

    # @testset "Poiseuille" begin
    #     include("Poiseuille.jl")
    # end

end