using Flower
using Test

using LsqFit

# Manufactured solutions

# atol =1e-13 in most places
#Constant Dirichlet everywhere
@testset "Poisson equation inside square constant Dirichlet: solve_poisson" begin
    include("poisson_square_constant_solve_poisson.jl")
end

@testset "Poisson equation inside square linear: solve_poisson" begin
    include("poisson_square_constant_solve_poisson_lin.jl")
end


#Linear function
@testset "Poisson equation inside square: solve_poisson" begin
    include("poisson_square_solve_poisson.jl")
end



# Skip

# @testset "Poisson equation inside square constant Dirichlet" begin
#     include("poisson_square_constant.jl")
# end


# @testset "Poisson equation inside square" begin
#     include("poisson_square.jl")
# end

# # Poisson equation from 
# @testset "Poisson equation inside circle" begin
#     include("poisson.jl")
# end




# @testset "Interpolation tests" begin
#     include("validation_test.jl")
# end
