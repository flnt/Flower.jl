abstract type NumericalParameters end

abstract type MutatingFields end

abstract type AbstractOptimizer end

@with_kw struct Numerical <: NumericalParameters
    n::Int64 = 128
    L0::Float64 = 2.0
    CFL::Float64 = 0.5
    Re::Float64 = 1.0
    TEND::Float64 = 0.0
    Δ::Float64 = L0/(n-1)
    shift::Float64 = 0.0
    shifted::Float64 = shift*Δ
    τ::Float64 = min(CFL*Δ^2*Re, CFL*Δ)
    max_iterations::Int = TEND÷τ
    current_i::Int = 1
    save_every::Int = 1
    reinit_every::Int = 1
    nb_reinit::Int = n÷8
    ϵ::Float64 = 0.00
    NB::Int64 = nb_reinit÷2
    Y::Array{Float64, 2} = [-L0 / 2 + i * Δ for i = 0:n-1, j = 0:n-1]
    X::Array{Float64, 2} = transpose(Y)
    Yu::Array{Float64, 2} = [Y[i,1] for i = 1:n, j = 1:n+1]
    Xu::Array{Float64, 2} = transpose([-L0 / 2 - Δ / 2 + i * Δ for i = 0:n, j = 1:n])
    Yv::Array{Float64, 2} = [-L0 / 2 - Δ / 2 + i * Δ for i = 0:n, j = 1:n]
    Xv::Array{Float64, 2} = transpose([X[1,i] for i = 1:n, j = 1:n+1])
    H::Array{Float64, 1} = [-L0 / 2 + i * Δ for i = 0:n-1]
    B::SMatrix{3, 3, Float64, 9} = @SMatrix [0.5 -1.0 0.5; -0.5 -0.0 0.5; 0.0 1.0 0.0]
    BT::SMatrix{3, 3, Float64, 9} = @SMatrix [0.5 -0.5 0.0; -1.0 -0.0 1.0; 0.5 0.5 0.0]
    T_inf::Float64 = 0.0
    u_inf::Float64 = 1.0
    v_inf::Float64 = 0.0
    θd::Float64 = 0.0
    ϵ_κ::Float64 = 0.0
    ϵ_V::Float64 = 0.0
    case::String = "notmycase"
    cases::String = "Planar, Sphere, Cylinder, Crystal, Mullins, Nothing, Airfoil, Square"
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

mutable struct TempArrays{T <: Float64} <: MutatingFields
    SCUTT::Array{T,1}
    LCUTT::Array{T,1}
    SCUTp::Array{T,1}
    LCUTp::Array{T,1}
    SCUTu::Array{T,1}
    LCUTu::Array{T,1}
    SCUTv::Array{T,1}
    LCUTv::Array{T,1}
    SCUTDx::Array{T,1}
    SCUTDy::Array{T,1}
    LCUTDx::Array{T,1}
    LCUTDy::Array{T,1}
    SCUTCT::Array{T,1}
    LCUTCT::Array{T,1}
    SCUTGxT::Array{T,1}
    LCUTGxT::Array{T,1}
    SCUTGyT::Array{T,1}
    LCUTGyT::Array{T,1}
    SCUTCu::Array{T,1}
    LCUTCu::Array{T,1}
    SCUTCv::Array{T,1}
    LCUTCv::Array{T,1}
    SOL::Array{T,3}
    LIQ::Array{T,3}
    SOLu::Array{T,3}
    LIQu::Array{T,3}
    SOLv::Array{T,3}
    LIQv::Array{T,3}
    sol_projection::Array{Gradient{T}, 2}
    liq_projection::Array{Gradient{T}, 2}
    sol_projectionu::Array{Gradient{T}, 2}
    liq_projectionu::Array{Gradient{T}, 2}
    sol_projectionv::Array{Gradient{T}, 2}
    liq_projectionv::Array{Gradient{T}, 2}
    sol_centroid::Array{Point{T},2}
    liq_centroid::Array{Point{T},2}
    mid_point::Array{Point{T},2}
    sol_centroidu::Array{Point{T},2}
    liq_centroidu::Array{Point{T},2}
    mid_pointu::Array{Point{T},2}
    sol_centroidv::Array{Point{T},2}
    liq_centroidv::Array{Point{T},2}
    mid_pointv::Array{Point{T},2}
    α::Array{T,2}
    αu::Array{T,2}
    αv::Array{T,2}
    LTS::SparseMatrixCSC{T, Int64}
    LTL::SparseMatrixCSC{T, Int64}
    LpS::SparseMatrixCSC{T, Int64}
    LpL::SparseMatrixCSC{T, Int64}
    LuS::SparseMatrixCSC{T, Int64}
    LuL::SparseMatrixCSC{T, Int64}
    LvS::SparseMatrixCSC{T, Int64}
    LvL::SparseMatrixCSC{T, Int64}
    AS::SparseMatrixCSC{T, Int64}
    AL::SparseMatrixCSC{T, Int64}
    BS::SparseMatrixCSC{T, Int64}
    BL::SparseMatrixCSC{T, Int64}
    LSA::SparseMatrixCSC{T, Int64}
    LSB::SparseMatrixCSC{T, Int64}
    GxpS::SparseMatrixCSC{T, Int64}
    GxpL::SparseMatrixCSC{T, Int64}
    GypS::SparseMatrixCSC{T, Int64}
    GypL::SparseMatrixCSC{T, Int64}
    DxuS::SparseMatrixCSC{T, Int64}
    DxuL::SparseMatrixCSC{T, Int64}
    DyvS::SparseMatrixCSC{T, Int64}
    DyvL::SparseMatrixCSC{T, Int64}
    ApS::SparseMatrixCSC{T, Int64}
    ApL::SparseMatrixCSC{T, Int64}
    AuS::SparseMatrixCSC{T, Int64}
    AuL::SparseMatrixCSC{T, Int64}
    AvS::SparseMatrixCSC{T, Int64}
    AvL::SparseMatrixCSC{T, Int64}
    CTS::SparseMatrixCSC{T, Int64}
    CTL::SparseMatrixCSC{T, Int64}
    GxTS::SparseMatrixCSC{T, Int64}
    GxTL::SparseMatrixCSC{T, Int64}
    GyTS::SparseMatrixCSC{T, Int64}
    GyTL::SparseMatrixCSC{T, Int64}
    ftcGxTS::SparseMatrixCSC{T, Int64}
    ftcGxTL::SparseMatrixCSC{T, Int64}
    ftcGyTS::SparseMatrixCSC{T, Int64}
    ftcGyTL::SparseMatrixCSC{T, Int64}
    CuS::SparseMatrixCSC{T, Int64}
    CuL::SparseMatrixCSC{T, Int64}
    CvS::SparseMatrixCSC{T, Int64}
    CvL::SparseMatrixCSC{T, Int64}
end

mutable struct Forward{T <: Float64} <: MutatingFields
    iso::Array{T, 2}
    isou::Array{T, 2}
    isov::Array{T, 2}
    u::Array{T,2}
    uu::Array{T,2}
    uv::Array{T,2}
    TS::Array{T,2}
    TL::Array{T,2}
    pS::Array{T,2}
    pL::Array{T,2}
    ϕS::Array{T,2}
    ϕL::Array{T,2}
    Gxm1S::Array{T,1}
    Gym1S::Array{T,1}
    Gxm1L::Array{T,1}
    Gym1L::Array{T,1}
    uS::Array{T,2}
    uL::Array{T,2}
    vS::Array{T,2}
    vL::Array{T,2}
    Tall::Array{T,2}
    DTS::Array{T,2}
    DTL::Array{T,2}
    V::Array{T,2}
    Vu::Array{T,2}
    Vv::Array{T,2}
    κ::Array{T, 2}
    κu::Array{T, 2}
    κv::Array{T, 2}
    usave::Array{T,3}
    uusave::Array{T,3}
    uvsave::Array{T,3}
    TSsave::Array{T,3}
    TLsave::Array{T,3}
    Tsave::Array{T,3}
    psave::Array{T,3}
    Uxsave::Array{T,3}
    Uysave::Array{T,3}
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

abstract type Grid end

struct GridCC <: Grid end
const gcc = GridCC()

struct GridFCx <: Grid end
const gfcx = GridFCx()

struct GridFCy <: Grid end
const gfcy = GridFCy()
