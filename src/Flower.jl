# __precompile__()

module Flower

using Reexport

import Base.<
import Base.isnan
import Base.*
import Base.+
import Base.-
import Base.abs
import Base.OneTo
import Base.ones
import Base.zeros
import Base.reshape

# import Base.union

# @reexport using Base.Union

@reexport using Printf
@reexport using Base.Threads
@reexport using LinearAlgebra
@reexport using SparseArrays

@reexport using Parameters
@reexport using StaticArrays
@reexport using OffsetArrays
@reexport using Roots
@reexport using SpecialFunctions
@reexport using GeometryBasics
@reexport using GeoInterface

# @reexport using LibGEOS
@reexport import LibGEOS #avoids warnings like:
# "WARNING: both MUMPS and Roots export "solve!"; uses of it in module Flower must be qualified"


# @reexport using HDF5

@reexport import YAML
@reexport using Test
# @reexport using MPI
@reexport using PropertyDicts
# @reexport using PrettyTables
# @reexport using Polynomials
@reexport using Peaks
# @reexport using Profile

# @reexport using ProfileCanvas

# @reexport using BenchmarkTools

@reexport using IterativeSolvers #for gmres

@reexport using Statistics #for mean
@reexport using MPI

# @reexport using MUMPS
@reexport import MUMPS


# #Long version
# @reexport using Printf
# # @reexport using Statistics
# @reexport using Base.Threads
# @reexport using LinearAlgebra
# @reexport using SparseArrays
# @reexport using BenchmarkTools
# @reexport using Parameters
# @reexport using StaticArrays
# @reexport using OffsetArrays
# @reexport using Roots
# @reexport using SpecialFunctions
# # @reexport using Test 
# # @reexport using CairoMakie
# # @reexport using LaTeXStrings
# # @reexport using Gnuplot
# @reexport using IterativeSolvers
# # @reexport using LsqFit
# # @reexport using Polynomials
# # @reexport using JLD
# # @reexport using Peaks
# @reexport using GeometryBasics
# @reexport using GeoInterface
# @reexport using LibGEOS



# @reexport using PackageCompiler

###################################################################################################
# For plotting with python
###################################################################################################
# @reexport using PyCall
# @reexport using LaTeXStrings
# @reexport using PyPlot

# const anim = PyNULL()
# # const plt = PyNULL()
# const mpl_colors = PyNULL()
# const mpl_tickers = PyNULL()
# const pd = PyNULL()

# # from matplotlib.colors import BoundaryNorm
# # from matplotlib.ticker import MaxNLocator

# function __init__()
#     copy!(anim, pyimport_conda("matplotlib.animation", "matplotlib"))
#     # copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib"))
#     copy!(mpl_colors, pyimport_conda("matplotlib.colors", "matplotlib"))
#     copy!(mpl_tickers, pyimport_conda("matplotlib.ticker", "matplotlib"))
#     copy!(pd, pyimport_conda("pandas", "pandas"))
# end
###################################################################################################

include("types.jl")
include("types_PDI.jl")
include("init.jl")
include("common.jl")
include("levelset.jl")
include("cutcell.jl")
include("operators.jl")
include("operators_coupled.jl")
include("heat.jl")
include("heat_coupled.jl")
include("navier_stokes_coupled.jl")
include("run.jl")
include("common_run.jl")
include("contact_line.jl")
include("optimize.jl")
# include("viz.jl")
include("tools.jl")
include("electrolysis_init.jl")
include("electrolysis_utils.jl")
include("electrolysis.jl")
include("electrolysis_2.jl")
include("electrolysis_operators.jl")
include("electrolysis_tests.jl")
include("electrolysis_print.jl")


include("run_profile.jl") #for profiling


###################################################################################################
# For plotting with python
###################################################################################################
# include("electrolysis_plot.jl")
###################################################################################################


# include("Plotpython.jl")
# using .Plotpython

# include("electrolysis_viz.jl")


export_all()
#this is a second test
end
