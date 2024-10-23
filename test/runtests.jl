using Flower
using Test


module testyamlfile4 #enables to perform a test with ARGS to give an input file
ARGS = String["../examples/levelset_Butler_two_LS.yml"]
include("../examples/main_current_folder.jl")
end

module testyamlfile3 #enables to perform a test with ARGS to give an input file
ARGS = String["../examples/levelset_Butler.yml"]
include("../examples/main_current_folder.jl")
end

module testyamlfile2 #enables to perform a test with ARGS to give an input file
ARGS = String["../examples/half_circle_Poiseuille.yml"]
include("../examples/main_current_folder.jl")
end

@testset "Case like Khalighi, bubble at the wall, two levelsets" begin
    using testyamlfile4
end
@testset "Case like Khalighi, bubble at the wall" begin
    using testyamlfile3
end
@testset "Simple tests with YAML file" begin
    using testyamlfile2
end
module testyamlfile #enables to perform a test with ARGS to give an input file
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