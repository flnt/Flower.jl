__precompile__()

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

@reexport using Printf
@reexport using Statistics
@reexport using Base.Threads
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using BenchmarkTools
@reexport using Parameters
@reexport using StaticArrays
@reexport using OffsetArrays
@reexport using Roots
@reexport using SpecialFunctions
@reexport using Test
# @reexport using Makie
# @reexport using CairoMakie
# @reexport using LaTeXStrings
# @reexport using Gnuplot
@reexport using IterativeSolvers
@reexport using LsqFit
@reexport using Polynomials
@reexport using JLD
@reexport using Peaks
@reexport using GeometryBasics
@reexport using GeoInterface
@reexport using LibGEOS
@reexport using PackageCompiler


@reexport using PyCall
@reexport using LaTeXStrings
@reexport using PyPlot

const anim = PyNULL()
# const plt = PyNULL()


function __init__()
    copy!(anim, pyimport_conda("matplotlib.animation", "matplotlib"))
    # copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib"))
end

include("types.jl")
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
include("electrolysis.jl")

include("electrolysis_plot.jl")


# include("Plotpython.jl")
# using .Plotpython

# include("electrolysis_viz.jl")


export_all()
#this is a second test
end
