abstract type NumericalParameters end

abstract type MutatingFields end

abstract type AbstractOptimizer end

@with_kw struct Numerical <: NumericalParameters
    n::Int64 = 128
    L0::Float64 = 2.0
    CFL::Float64 = 0.5
    TEND::Float64 = 0.0
    Δ::Float64 = L0/(n-1)
    shift::Float64 = 0.0
    shifted::Float64 = shift*Δ
    τ::Float64 = CFL*(Δ)^2
    max_iterations::Int = TEND÷τ
    current_i::Int = 1
    reinit_every::Int = 1
    nb_reinit::Int = n÷8
    ϵ::Float64 = 0.00
    NB::Int64 = nb_reinit÷2
    Y::Array{Float64, 2} = [-L0 / 2 + i * Δ for i = 0:n-1, j = 0:n-1]
    X::Array{Float64, 2} = transpose(Y)
    H::Array{Float64, 1} = [-L0 / 2 + i * Δ for i = 0:n-1]
    B::SMatrix{3, 3, Float64, 9} = @SMatrix [0.5 -1.0 0.5; -0.5 -0.0 0.5; 0.0 1.0 0.0]
    BT::SMatrix{3, 3, Float64, 9} = @SMatrix [0.5 -0.5 0.0; -1.0 -0.0 1.0; 0.5 0.5 0.0]
    T_inf::Float64 = 0.0
    θd::Float64 = 0.0
    ϵ_κ::Float64 = 0.0
    ϵ_V::Float64 = 0.0
    case::String = "notmycase"
    cases::String = "Planar, Sphere, Crystal, Mullins"
    A::Float64 = 0.05
    N::Int64 = 2
    R::Float64 = 0.5
    m::Int64 = 4
    θ₀::Float64 = pi/4
    aniso::Bool = false
end

@with_kw mutable struct Indices{T <: Int64} <: NumericalParameters
    inside::CartesianIndices{2, Tuple{OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}, OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}}}
    periodicL::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    periodicR::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    periodicB::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    periodicT::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_left::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_bottom::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_right::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_top::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
end

struct Point{T <: Number}
   x::T
   y::T
end

mutable struct Gradient{T <: Float64}
    flag::Bool
    angle::T
    mid_point::Point{T}
    point1::Point{T}
    point2::Point{T}
    d1::T
    d2::T
    pos::Point{T}
end

mutable struct TempArrays{T <: Float64} <: MutatingFields
    SCUT::Array{T,1}
    LCUT::Array{T,1}
    SOL::Array{T,3}
    LIQ::Array{T,3}
    sol_projection::Array{Gradient{T}, 2}
    liq_projection::Array{Gradient{T}, 2}
    AS::SparseMatrixCSC{T, Int64}
    AL::SparseMatrixCSC{T, Int64}
    BS::SparseMatrixCSC{T, Int64}
    BL::SparseMatrixCSC{T, Int64}
    LSA::SparseMatrixCSC{T, Int64}
    LSB::SparseMatrixCSC{T, Int64}
end

mutable struct Forward{T <: Float64} <: MutatingFields
    iso::Array{T, 2}
    u::Array{T,2}
    TS::Array{T,2}
    TL::Array{T,2}
    Tall::Array{T,2}
    V::Array{T,2}
    κ::Array{T, 2}
    usave::Array{T,3}
    TSsave::Array{T,3}
    TLsave::Array{T,3}
    Tsave::Array{T,3}
    Vsave::Array{T,3}
    κsave::Array{T,3}
    lengthsave::Array{T,1}
end

mutable struct Desired{T <: Float64} <: AbstractOptimizer
    x_desired::Vector{T}
    u::Matrix{T}
    usave::Array{T, 3}
    TL::Matrix{T}
    TS::Matrix{T}
end

mutable struct Optim_parameters{T <: Float64, N <: Int64} <: AbstractOptimizer
    nprobes::N
    ind::Vector{N}
    bc_indices::Vector{CartesianIndex{2}}
    γ::Vector{T}
    p::Vector{Vector{T}}
    TLsave::Vector{Matrix{T}}
    TSsave::Vector{Matrix{T}}
    usave::Vector{Array{T, 3}}
end

mutable struct my_Adjoint{T <: Float64} <: MutatingFields
    iso::Array{T, 2}
    u::Array{T,2}
    TS::Array{T,2}
    TL::Array{T,2}
    V::Array{T,2}
    κ::Array{T, 2}
end

@with_kw mutable struct Boundary{L, T, N} <: NumericalParameters
   ind::L = ([CartesianIndex(0,0)], [CartesianIndex(0,0)])
   f::Function = neumann
   val::N = 0.0
   λ::T = 0.0
end

@with_kw mutable struct Boundaries <: NumericalParameters
    left::Boundary = Boundary()
    right::Boundary = Boundary()
    bottom::Boundary = Boundary()
    top::Boundary = Boundary()
end
