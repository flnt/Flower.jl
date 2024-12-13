using Flower
using Test

import Test: Test, finish
using Test: DefaultTestSet, Broken
using Test: parse_testset_args

"""
Skip a testset

Use `@testset_skip` to replace `@testset` for some tests which should be skipped.

Usage
-----
Replace `@testset` with `@testset "reason"` where `"reason"` is a string saying why the
test should be skipped (which should come before the description string, if that is
present).
"""
macro testset_skip(args...)
    isempty(args) && error("No arguments to @testset_skip")
    length(args) < 2 && error("First argument to @testset_skip giving reason for "
                              * "skipping is required")

    skip_reason = args[1]

    desc, testsettype, options = parse_testset_args(args[2:end-1])

    ex = quote
        # record the reason for the skip in the description, and mark the tests as
        # broken, but don't run tests
        local ts = DefaultTestSet(string($desc, " - ", $skip_reason))
        push!(ts.results, Broken(:skipped, "skipped tests"))
        local ret = finish(ts)
        ret
    end

    return ex
end

# using LsqFit


# Gradient
# @testset "Gradient" begin
#     include("gradient.jl")
# end
# @testset_skip "Check localization of u component of gradient" "Gradient" begin
#     include("gradient.jl")
# end


# Manufactured solutions

# Poisson


# module testyamlfile6 #enables to perform a test with ARGS to give an input file
# ARGS = String["../Flower.jl/test/poisson_circular_interface_wall.yml"]
# include("poisson_circular_interface_wall.jl")
# end

#TODO

# Poisson without interface

module testyamlfile #enables to perform a test with ARGS to give an input file
ARGS = String["../Flower.jl/test/poisson_no_interface.yml"]
include("poisson_no_interface_right_Neumann.jl")
end

module testyamlfile #enables to perform a test with ARGS to give an input file
ARGS = String["../Flower.jl/test/poisson_no_interface.yml"]
include("poisson_no_interface.jl")
end

module testyamlfile2 #enables to perform a test with ARGS to give an input file
ARGS = String["../Flower.jl/test/poisson_circular_interface_Dirichlet.yml"]
include("poisson_circular_interface_Dirichlet.jl")
end


module testyamlfile5 #enables to perform a test with ARGS to give an input file
ARGS = String["../Flower.jl/test/poisson_circular_interface_Neumann.yml"]
include("poisson_circular_interface_Neumann.jl")
end

# module testyamlfile4 #enables to perform a test with ARGS to give an input file
# ARGS = String["../Flower.jl/test/poisson_circular_interface_wall.yml"]
# include("poisson_circular_interface_wall.jl")
# end

# module testyamlfile3 #enables to perform a test with ARGS to give an input file
# ARGS = String["../Flower.jl/test/poisson_circular_interface.yml"]
# include("poisson_square_circle_solve_poisson_cos_cos.jl")
# end




# module testyamlfile #enables to perform a test with ARGS to give an input file
# ARGS = String["../Flower.jl/test/poisson_square_solve_poisson.yml"]
# include("poisson_square_solve_poisson.jl")
# end


# @testset "Poisson equation inside square with interface: solve_poisson coscos" begin
#     using testyamlfile3
# end
# @testset "Poisson equation inside square: solve_poisson coscos" begin
#     using testyamlfile2
# end

# @testset "Poisson equation inside square: solve_poisson" begin
#     using testyamlfile
# end

# #Convergence study
# @testset "Poisson equation inside square: solve_poisson" begin
#     include("poisson_square_solve_poisson.jl")
# end

# #Constant Dirichlet everywhere
# @testset "Poisson equation inside square constant Dirichlet: solve_poisson" begin
#     include("poisson_square_constant_solve_poisson.jl")
# end

# #Linear function
# @testset "Poisson equation inside square linear: solve_poisson" begin
#     include("poisson_square_solve_poisson_lin.jl")
# end



# Heat equation
#TODO 



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
