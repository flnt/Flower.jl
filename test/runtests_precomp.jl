
import Flower
include(joinpath(pkgdir(Flower), "test", "validation_test_precomp.jl"))

# # using Flower
# import Flower
# # using Test
# include("validation_test_precomp.jl")

# @testset "Simple tests" begin
#     include("validation_test_precomp.jl")
# end

# @testset "Example tests" begin

#     # @testset "Zero velocity" begin
#     #     include("zerovel.jl")
#     # end

#     # @testset "Constant velocity" begin
#     #     include("constantvel.jl")
#     # end

#     # @testset "Poiseuille" begin
#     #     include("Poiseuille.jl")
#     # end

# end