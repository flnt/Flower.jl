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
    Δ::Float64 = min(diff(x)..., diff(y)...)
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
    σ::Float64 = 1.0
    case::String = "notmycase"
    cases::String = "Planar, Sphere, Cylinder, Ellipse, Crystal, Mullins, Nothing, Airfoil, Jet, Drop"
    A::Float64 = 0.05
    N::Int64 = 2
    R::Float64 = 0.5
    m::Int64 = 4
    θ₀::Float64 = pi/4
    g::Float64 = 0.0
    β::Float64 = 0.0
    x_airfoil::Array{Float64} = [0.0]
    y_airfoil::Array{Float64} = [0.0]
    aniso::Bool = false
    subdomains::Int64 = 2
    overlaps::Int64 = 1
    tolu::Float64 = 1.0e-6
    tolp::Float64 = 1.0e-6
    tolt::Float64 = 1.0e-6
    # contact angle models parameters
    Ca::Float64 = 0.00 # Capillary number
    εCA::Float64 = 0.00 # width
    λCA::Float64 = 0.00 # slip lenght
    θe::Float64 = 0.00 # prescribed contact angle
end

@with_kw mutable struct Indices{T <: Int64} <: NumericalParameters
    all_indices::Array{CartesianIndex{2},2}
    inside::CartesianIndices{2, Tuple{OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}, OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}}}
    periodic_x::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    periodic_y::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_left::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_bottom::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_right::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    b_top::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    MIXED::Vector{CartesianIndex{2}}
    MIXED_ext::Vector{CartesianIndex{2}}
    LIQUID::Vector{CartesianIndex{2}}
    LIQUID_ext::Vector{CartesianIndex{2}}
    SOLID::Vector{CartesianIndex{2}}
    SOLID_ext::Vector{CartesianIndex{2}}
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
    fresh::Array{Bool,2}
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
    Young::Array{T,2} # maybe transform in sparse ?
    κ::Array{T,2}
    V::Array{T,2}
    LSA::SparseMatrixCSC{T,N}
    LSB::SparseMatrixCSC{T,N}
    domdec::DomainDecomposition{TS,3,D,A,O,W} where {TS,D,A,O,W}
    pou::DomainDecomposedVector{T,3,DomainDecomposition{TS,3,D,A,O,W},A2} where {TS,D,A,O,W,A2}
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

mutable struct OperatorsCoupled{T <: Float64} <: MutatingFields
    AxT::SparseMatrixCSC{T,Int64}
    AyT::SparseMatrixCSC{T,Int64}
    Bx::SparseMatrixCSC{T,Int64}
    By::SparseMatrixCSC{T,Int64}
    BxT::SparseMatrixCSC{T,Int64}
    ByT::SparseMatrixCSC{T,Int64}
    Hx::SparseMatrixCSC{T,Int64}
    Hy::SparseMatrixCSC{T,Int64}
    HxT::SparseMatrixCSC{T,Int64}
    HyT::SparseMatrixCSC{T,Int64}
    tmp_x::SparseMatrixCSC{T,Int64}
    tmp_y::SparseMatrixCSC{T,Int64}
    M::Diagonal{T,Vector{T}}
    iMx::Diagonal{T,Vector{T}}
    iMy::Diagonal{T,Vector{T}}
    χ::Diagonal{T,Vector{T}}
    Rx::SparseMatrixCSC{T,Int64}
    Ry::SparseMatrixCSC{T,Int64}
    Gx::SparseMatrixCSC{T,Int64}
    Gy::SparseMatrixCSC{T,Int64}
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
    Dϕ::Array{T,2}
    Du::Array{T,2}
    Dv::Array{T,2}
    TD::Array{T,1}
    pD::Array{T,1}
    ϕD::Array{T,1}
    uD::Array{T,1}
    vD::Array{T,1}
    ucorrD::Array{T,1}
    vcorrD::Array{T,1}
end

mutable struct Forward{T <: Float64} <: MutatingFields
    T::Array{T,3}
    u::Array{T,3}
    ux::Array{T,3}
    uy::Array{T,3}
    V::Array{T,3}
    κ::Array{T,3}
    length::Array{T,1}
    t::Array{T,1}
    Cd::Array{T,1}
    Cl::Array{T,1}
end

mutable struct ForwardPhase{T <: Float64} <: MutatingFields
    T::Array{T,3}
    p::Array{T,3}
    ϕ::Array{T,3}
    u::Array{T,3}
    v::Array{T,3}
    ucorr::Array{T,3}
    vcorr::Array{T,3}
    TD::Array{T,2}
    pD::Array{T,2}
    ucorrD::Array{T,2}
    vcorrD::Array{T,2}
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

mutable struct adjoint_phase{T <: Float64} <: MutatingFields
    TD::Array{T,2}
    pD::Array{T,2}
    u::Array{T,2}
    v::Array{T,2}
    ucorrD::Array{T,2}
    vcorrD::Array{T,2}
end

mutable struct adjoint_fields{T <: Float64} <: MutatingFields
    u::Array{T,2}
    phS::adjoint_phase{T}
    phL::adjoint_phase{T}
end

mutable struct adjoint_derivatives{T <: Float64} <: MutatingFields
    RheatS_ls::SparseMatrixCSC{T,Int64} # Derivative of solid phase heat eq. wrt. levelset
    RheatL_ls::SparseMatrixCSC{T,Int64} # Derivative of liquid phase heat eq. wrt. levelset
    RlsS_ls::SparseMatrixCSC{T,Int64} # Derivative of Stefan levelset advection eq. wrt. levelset
    RlsS_TS::SparseMatrixCSC{T,Int64} # Derivative of Stefan levelset advection eq. wrt. solid temperature 
    RlsS_TL::SparseMatrixCSC{T,Int64} # Derivative of Stefan levelset advection eq. wrt. liquid temperature
    RucorrS_ls0::SparseMatrixCSC{T,Int64} # Derivative of solid phase prediction step in x wrt. levelset at prev. it.
    RucorrS_ls1::SparseMatrixCSC{T,Int64} # Derivative of solid phase prediction step in x wrt. levelset
    RucorrL_ls0::SparseMatrixCSC{T,Int64} # Derivative of liquid phase prediction step in x wrt. levelset at prev. it.
    RucorrL_ls1::SparseMatrixCSC{T,Int64} # Derivative of liquid phase prediction step in x wrt. levelset
    RvcorrS_ls0::SparseMatrixCSC{T,Int64} # Derivative of solid phase prediction step in y wrt. levelset at prev. it.
    RvcorrS_ls1::SparseMatrixCSC{T,Int64} # Derivative of solid phase prediction step in y wrt. levelset
    RvcorrL_ls0::SparseMatrixCSC{T,Int64} # Derivative of liquid phase prediction step in y wrt. levelset at prev. it.
    RvcorrL_ls1::SparseMatrixCSC{T,Int64} # Derivative of liquid phase prediction step in y wrt. levelset
    RpS_ls::SparseMatrixCSC{T,Int64} # Derivative of solid phase Poisson eq. wrt. levelset
    RpL_ls::SparseMatrixCSC{T,Int64} # Derivative of liquid phase Poisson eq. wrt. levelset
    RuS_ls::SparseMatrixCSC{T,Int64} # Derivative of solid phase projection step in x wrt. levelset
    RuL_ls::SparseMatrixCSC{T,Int64} # Derivative of liquid phase projection step in x wrt. levelset
    RvS_ls::SparseMatrixCSC{T,Int64} # Derivative of solid phase projection step in y wrt. levelset
    RvL_ls::SparseMatrixCSC{T,Int64} # Derivative of liquid phase projection step in y wrt. levelset
    RlsFS_ls::SparseMatrixCSC{T,Int64} # Derivative of Free Surface levelset advection eq. wrt. levelset
    RlsFS_ucorrS::SparseMatrixCSC{T,Int64} # Derivative of Free Surface levelset advection eq. wrt. solid horizontal velocity 
    RlsFS_ucorrL::SparseMatrixCSC{T,Int64} # Derivative of Free Surface levelset advection eq. wrt. liquid horizontal velocity 
    RlsFS_vcorrS::SparseMatrixCSC{T,Int64} # Derivative of Free Surface levelset advection eq. wrt. solid vertical velocity
    RlsFS_vcorrL::SparseMatrixCSC{T,Int64} # Derivative of Free Surface levelset advection eq. wrt. liquid vertical velocity
end

abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition end
const dir = Dirichlet()
struct Neumann <: BoundaryCondition end
const neu = Neumann()
struct Robin <: BoundaryCondition end
const rob = Robin()
struct Periodic <: BoundaryCondition end
const per = Periodic()
struct Navier <: BoundaryCondition end
const nav = Navier()

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
