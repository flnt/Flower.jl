using Flower
using Test


module testyamlfile
ARGS = String["../examples/half_circle_imposed_radius.yml"]
include("../examples/main_current_folder.jl")
end
@testset "Simple tests with YAML file" begin

    # module testyamlfiles
    # localARGS = String["../examples/half_circle_imposed_radius.yml"]
    # ARGS = String["../examples/half_circle_imposed_radius.yml"]
    # include("../examples/main_current_folder.jl")
    # end
    localARGS = String["../examples/half_circle_imposed_radius.yml"]
    ARGS = String["../examples/half_circle_imposed_radius.yml"]
    # localARGS = String["half_circle_imposed_radius.yml"]
    # ARGS = String["half_circle_imposed_radius.yml"]
    # print(localARGS)
    # include("../examples/main_current_folder.jl")

    # localARGS = String["../examples/half_circle_imposed_radius.yml"]; include("../examples/main_current_folder.jl")

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