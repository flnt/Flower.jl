using Flower
using Test


module testyamlfile2 #enables to performa test with ARGS to give an input file
ARGS = String["../examples/half_circle_Poiseuille.yml"]
include("../examples/main_current_folder.jl")
end
@testset "Simple tests with YAML file" begin
    using testyamlfile2
end
module testyamlfile #enables to performa test with ARGS to give an input file
ARGS = String["../examples/half_circle_imposed_radius.yml"]
include("../examples/main_current_folder.jl")
end
@testset "Simple tests with YAML file" begin
    using testyamlfile
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