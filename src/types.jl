abstract type NumericalParameters end

abstract type MutatingFields end

abstract type AbstractOptimizer end

@with_kw struct Numerical <: NumericalParameters
    CFL::Float64 = 0.5
    Re::Float64 = 1.0
    TEND::Float64 = 0.0
    x::Vector{Float64} = [-0.5 - 1/127 / 2 + i * 1/127 for i = 0:128]
    y::Vector{Float64} = [-0.5 - 1/127 / 2 + i * 1/127 for i = 0:128]
    L0::Float64 = max(x[end]-x[1], y[end]-y[1])
    Δ::Float64 = min(diff(x)...)
    shift::Float64 = 0.0
    shifted::Float64 = shift*Δ
    τ::Float64 = min(CFL*Δ^2*Re, CFL*Δ)
    max_iterations::Int = TEND÷τ
    current_i::Int = 1
    save_every::Int = 1
    reinit_every::Int = 1
    nb_reinit::Int = length(x)÷8
    ϵ::Float64 = 0.00
    NB::Int64 = nb_reinit÷2
    T_inf::Float64 = 0.0
    u_inf::Float64 = 1.0
    v_inf::Float64 = 0.0
    θd::Float64 = 0.0
    ϵ_κ::Float64 = 0.0
    ϵ_V::Float64 = 0.0
    case::String = "notmycase"
    cases::String = "Planar, Sphere, Cylinder, Crystal, Mullins, Nothing, Airfoil"
    A::Float64 = 0.05
    N::Int64 = 2
    R::Float64 = 0.5
    m::Int64 = 4
    θ₀::Float64 = pi/4
    x_airfoil::Array{Float64} = [0.0]
    y_airfoil::Array{Float64} = [0.0]
    aniso::Bool = false
end

@with_kw mutable struct Indices{T <: Int64} <: NumericalParameters
    all_indices::Array{CartesianIndex{2},2}
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

struct Line{T <: Number}
    p1::Point{T}
    p2::Point{T}
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

mutable struct GeometricInfo{T} <: MutatingFields
    cap::Array{T,3}
    dcap::Array{T,3}
    projection::Array{Gradient{T},2}
    centroid::Array{Point{T},2}
    emptied::Array{Bool,2}
end

abstract type Grid <: MutatingFields end

struct GridCC <: Grid end
struct GridFCx <: Grid end
struct GridFCy <: Grid end

mutable struct Mesh{G,T,N} <: Grid where {G<:Grid}
    x_nodes::Vector{T}
    y_nodes::Vector{T}
    x::Array{T,2}
    y::Array{T,2}
    nx::N
    ny::N
    dx::Array{T,2}
    dy::Array{T,2}
    ind::Indices{N}
    u::Array{T,2}
    iso::Array{T,2}
    faces::Array{T,3}
    geoS::GeometricInfo{T}
    geoL::GeometricInfo{T}
    mid_point::Array{Point{T},2}
    cut_points::Array{Vector{Point{T}},2}
    α::Array{T,2}
    κ::Array{T,2}
    V::Array{T,2}
    LSA::SparseMatrixCSC{T,N}
    LSB::SparseMatrixCSC{T,N}
end

mutable struct Operators{T <: Float64} <: MutatingFields
    CUTT::Array{T,1}
    CUTp::Array{T,1}
    CUTu::Array{T,1}
    CUTv::Array{T,1}
    CUTDx::Array{T,1}
    CUTDy::Array{T,1}
    CUTCT::Array{T,1}
    CUTGxT::Array{T,1}
    CUTGyT::Array{T,1}
    CUTGxp::Array{T,1}
    CUTGyp::Array{T,1}
    CUTGxϕ::Array{T,1}
    CUTGyϕ::Array{T,1}
    CUTCu::Array{T,1}
    CUTCv::Array{T,1}
    LT::SparseMatrixCSC{T,Int64}
    Lp::SparseMatrixCSC{T,Int64}
    Lu::SparseMatrixCSC{T,Int64}
    Lv::SparseMatrixCSC{T,Int64}
    A::SparseMatrixCSC{T,Int64}
    B::SparseMatrixCSC{T,Int64}
    Gxp::SparseMatrixCSC{T,Int64}
    Gyp::SparseMatrixCSC{T,Int64}
    Gxϕ::SparseMatrixCSC{T,Int64}
    Gyϕ::SparseMatrixCSC{T,Int64}
    Dxu::SparseMatrixCSC{T,Int64}
    Dyv::SparseMatrixCSC{T,Int64}
    Ap::SparseMatrixCSC{T,Int64}
    Au::SparseMatrixCSC{T,Int64}
    Av::SparseMatrixCSC{T,Int64}
    CT::SparseMatrixCSC{T,Int64}
    GxT::SparseMatrixCSC{T,Int64}
    GyT::SparseMatrixCSC{T,Int64}
    ftcGxT::SparseMatrixCSC{T,Int64}
    ftcGyT::SparseMatrixCSC{T,Int64}
    Cu::SparseMatrixCSC{T,Int64}
    Cv::SparseMatrixCSC{T,Int64}
    E11::SparseMatrixCSC{T,Int64}
    E12_x::SparseMatrixCSC{T,Int64}
    E12_y::SparseMatrixCSC{T,Int64}
    E22::SparseMatrixCSC{T,Int64}
    utp::SparseMatrixCSC{T,Int64}
    vtp::SparseMatrixCSC{T,Int64}
end

mutable struct Phase{T <: Float64} <: MutatingFields
    T::Array{T,2}
    p::Array{T,2}
    ϕ::Array{T,2}
    Gxm1::Array{T,1}
    Gym1::Array{T,1}
    u::Array{T,2}
    v::Array{T,2}
    ucorr::Array{T,2}
    vcorr::Array{T,2}
    DT::Array{T,2}
    Du::Array{T,2}
    Dv::Array{T,2}
    tmp::Array{T,2}
end

mutable struct Forward{T <: Float64} <: MutatingFields
    Tall::Array{T,2}
    usave::Array{T,3}
    uusave::Array{T,3}
    uvsave::Array{T,3}
    TSsave::Array{T,3}
    TLsave::Array{T,3}
    Tsave::Array{T,3}
    psave::Array{T,3}
    ϕsave::Array{T,3}
    Uxsave::Array{T,3}
    Uysave::Array{T,3}
    Uxcorrsave::Array{T,3}
    Uycorrsave::Array{T,3}
    Vsave::Array{T,3}
    κsave::Array{T,3}
    lengthsave::Array{T,1}
    Cd::Array{T,1}
    Cl::Array{T,1}
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
    DTS::Array{T,2}
    DTL::Array{T,2}
    V::Array{T,2}
    κ::Array{T, 2}
end

abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition end
const dir = Dirichlet()
struct Neumann <: BoundaryCondition end
const neu = Neumann()
struct Periodic <: BoundaryCondition end
const per = Periodic()

@with_kw mutable struct Boundary{BC, L, T, N} <: NumericalParameters
   t::BC = neu
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
