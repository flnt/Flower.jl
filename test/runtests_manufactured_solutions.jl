using Flower
using Test

using LsqFit

# Manufactured solutions

@testset "Poisson equation inside square" begin
    include("poisson_square.jl")
end

# Poisson equation from 
@testset "Poisson equation inside circle" begin
    include("poisson.jl")
end




# @testset "Interpolation tests" begin
#     include("validation_test.jl")
# end
