__precompile__()

module Flower

using Reexport

import Base.<
import Base.isnan
import Base.*
import Base.+
import Base.-
import Base.abs

@reexport using Printf
@reexport using Statistics
@reexport using Base.Threads
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using BenchmarkTools
@reexport using Parameters
@reexport using Optim
@reexport using StaticArrays
@reexport using OffsetArrays
@reexport using Roots
@reexport using SpecialFunctions
@reexport using Test
@reexport using Makie
@reexport using CairoMakie
@reexport using Gnuplot
@reexport using IterativeSolvers
@reexport using PETSc
@reexport using LsqFit
@reexport using Polynomials

include("types.jl")
include("init.jl")
include("common.jl")
include("levelset.jl")
include("cutcell.jl")
include("heat.jl")
include("run.jl")
include("optimize.jl")
include("viz.jl")

export_all()

end
