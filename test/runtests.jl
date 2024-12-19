using Flower
using Test


# @testset "Interpolation tests" begin
#     include("validation_test.jl")
# end

# Test Poiseuille: 
# Test that with factor 4/3 we find the right pressure gradient
module testyamlfile5 #enables to perform a test with ARGS to give an input file
ARGS = String["../Flower.jl/test/channel_Dirichlet_pressure.yml"]
include("../examples/main_concise.jl")
end


# Manufactured solutions

# Poisson equation from  
@testset "Poisson equation inside circle" begin
    include("poisson.jl")
end

@testset "Poisson equation inside square" begin
    include("poisson_square.jl")
end


module testyamlfile3 #enables to perform a test with ARGS to give an input file
ARGS = String["../examples/levelset_Butler.yml"]
include("../examples/main_current_folder.jl")
end

module testyamlfile2 #enables to perform a test with ARGS to give an input file
ARGS = String["../examples/half_circle_Poiseuille.yml"]
include("../examples/main_current_folder.jl")
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


# Deprecated: 2 LS

# module testyamlfile6 #enables to perform a test with ARGS to give an input file
# ARGS = String["../examples/test_levelset_Butler_two_LS_2.yml"]
# include("../examples/main_current_folder.jl")
# end

# module testyamlfile5 #enables to perform a test with ARGS to give an input file
# ARGS = String["../examples/test_levelset_Butler_two_LS.yml"]
# include("../examples/main_current_folder.jl")
# end

# module testyamlfile4 #enables to perform a test with ARGS to give an input file
# ARGS = String["../examples/levelset_Butler_two_LS.yml"]
# include("../examples/main_current_folder.jl")
# end

# @testset "Test Poisson with 2 LS, bubble (here circle with Butler-Volmer BC) at wall" begin
#     using testyamlfile6
# end

# @testset "Test Poisson with 2 LS, bubble (here circle with Butler-Volmer BC) away from wall" begin
#     using testyamlfile5
# end

# @testset "Case like Khalighi, bubble at the wall, two levelsets" begin
#     using testyamlfile4
# end