using Flower
using Test


@testset "Simple tests with YAML file" begin
    ARGS = String["half_circle_imposed_radius.yml"]
    include("../examples/main_current_folder.jl")
end
@testset "Simple tests" begin
    include("validation_test.jl")
end

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