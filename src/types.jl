abstract type NumericalParameters end

abstract type MutatingFields end

abstract type AbstractOptimizer end

@with_kw struct Numerical{T <: Real, D <: Integer} <: NumericalParameters
    CFL::T = 0.5
    Re::T = 1.0
    TEND::T = 0.0
    x::Union{Vector{T},LinRange{T,D}} = [-0.5 - 1/127 / 2 + i * 1/127 for i = 0:128]
    y::Union{Vector{T},LinRange{T,D}} = [-0.5 - 1/127 / 2 + i * 1/127 for i = 0:128]
    L0::T = max(x[end]-x[1], y[end]-y[1])
    Δ::T = min(diff(x)..., diff(y)...)
    shift::T = 0.0
    shifted::T = shift*Δ
    τ::T = min(CFL*Δ^2*Re, CFL*Δ)
    max_iterations::D = TEND÷τ
    current_i::D = 1
    save_every::D = 1
    reinit_every::D = 1
    nb_reinit::D = length(x)÷8
    δreinit::T = 0.01
    ϵ::T = 0.00
    NB::D = nb_reinit÷2
    T_inf::T = 0.0
    u_inf::T = 1.0
    v_inf::T = 0.0
    uD::T = 0.0
    vD::T = 0.0
    θd::T = 0.0
    ϵ_κ::T = 0.0
    ϵ_V::T = 0.0
    σ::T = 0.0
    case::String = "notmycase"
    cases::String = "Planar, Sphere, Cylinder, Ellipse, Crystal, Mullins, Nothing, Airfoil, Jet, Drop"
    A::T = 0.05
    N::D = 2
    R::T = 0.5
    m::D = 4
    θ₀::T = pi/4
    g::T = 0.0
    β::T = 0.0
    n_ext_cl::D = 5
    x_airfoil::Array{T} = [0.0]
    y_airfoil::Array{T} = [0.0]
    aniso::Bool = false
    nLS::D = 1
    _nLS::D = nLS == 1 ? 1 : nLS + 1
    subdomains::D = 2
    overlaps::D = 1
end

@with_kw struct Indices{T <: Integer} <: NumericalParameters
    all_indices::Array{CartesianIndex{2},2}
    inside::CartesianIndices{2, Tuple{OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}, OffsetArrays.IdOffsetRange{T, Base.OneTo{T}}}}
    periodic_x::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
    periodic_y::Tuple{Vector{CartesianIndex{2}}, Vector{CartesianIndex{2}}}
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

mutable struct Gradient{T <: Real}
    flag::Bool
    angle::T
    mid_point::Point{T}
    point1::Point{T}
    point2::Point{T}
    d1::T
    d2::T
    pos::Point{T}
end

struct GeometricInfo{T <: Real} <: MutatingFields
    cap::Array{T,3}
    dcap::Array{T,3}
    projection::Array{Gradient{T},2}
    centroid::Array{Point{T},2}
    vertices::Array{Vector{Point{T}},2}
    emptied::Array{Bool,2}
    fresh::Array{Bool,2}
end

mutable struct Levelset{T,N}
    u::Matrix{T}
    iso::Matrix{T}
    faces::Array{T,3}
    geoS::GeometricInfo{T}
    geoL::GeometricInfo{T}
    mid_point::Matrix{Point{T}}
    cut_points::Matrix{Vector{Point{T}}}
    α::Matrix{T}
    κ::Matrix{T}
    A::SparseMatrixCSC{T,N}
    B::SparseMatrixCSC{T,N}
    MIXED::Vector{CartesianIndex{2}}
    LIQUID::Vector{CartesianIndex{2}}
    SOLID::Vector{CartesianIndex{2}}
    cl::Vector{CartesianIndex{2}}
end

abstract type Grid <: MutatingFields end

struct GridCC <: Grid end
struct GridFCx <: Grid end
struct GridFCy <: Grid end

struct Mesh{G,T,N} <: Grid where {G<:Grid}
    x_nodes::Vector{T}
    y_nodes::Vector{T}
    x::Array{T,2}
    y::Array{T,2}
    nx::N
    ny::N
    dx::Array{T,2}
    dy::Array{T,2}
    LS::Vector{Levelset}
    ind::Indices{N}
    V::Array{T,2}
    domdec::DomainDecomposition{TS,3,D,A,O,W} where {TS,D,A,O,W}
    pou::DomainDecomposedVector{T,3,DomainDecomposition{TS,3,D,A,O,W},A2} where {TS,D,A,O,W,A2}
end

struct OperatorsConvection{T <: Real, D <: Integer} <: MutatingFields
    CUTCT::Array{T,1}
    CUTCu::Array{T,1}
    CUTCv::Array{T,1}
    CT::SparseMatrixCSC{T,D}
    Cu::SparseMatrixCSC{T,D}
    Cv::SparseMatrixCSC{T,D}
    E11::SparseMatrixCSC{T,D}
    E12_x::SparseMatrixCSC{T,D}
    E12_y::SparseMatrixCSC{T,D}
    E22::SparseMatrixCSC{T,D}
end

struct Operators{T <: Real, D <: Integer} <: MutatingFields
    AxT::SparseMatrixCSC{T,D}
    AyT::SparseMatrixCSC{T,D}
    Bx::SparseMatrixCSC{T,D}
    By::SparseMatrixCSC{T,D}
    BxT::SparseMatrixCSC{T,D}
    ByT::SparseMatrixCSC{T,D}
    Hx::Vector{SparseMatrixCSC{T,D}}
    Hy::Vector{SparseMatrixCSC{T,D}}
    HxT::Vector{SparseMatrixCSC{T,D}}
    HyT::Vector{SparseMatrixCSC{T,D}}
    tmp_x::SparseMatrixCSC{T,D}
    tmp_y::SparseMatrixCSC{T,D}
    M::Diagonal{T,Vector{T}}
    iMx::Diagonal{T,Vector{T}}
    iMy::Diagonal{T,Vector{T}}
    χ::Vector{Diagonal{T,Vector{T}}}
    Rx::SparseMatrixCSC{T,D}
    Ry::SparseMatrixCSC{T,D}
    Gx::Vector{SparseMatrixCSC{T,D}}
    Gy::Vector{SparseMatrixCSC{T,D}}
    Hx_b::SparseMatrixCSC{T,D}
    Hy_b::SparseMatrixCSC{T,D}
    HxT_b::SparseMatrixCSC{T,D}
    HyT_b::SparseMatrixCSC{T,D}
    iMx_b::SparseMatrixCSC{T,D}
    iMy_b::SparseMatrixCSC{T,D}
    iMx_bd::Diagonal{T,Vector{T}}
    iMy_bd::Diagonal{T,Vector{T}}
    Gx_b::SparseMatrixCSC{T,D}
    Gy_b::SparseMatrixCSC{T,D}
    χ_b::Diagonal{T,Vector{T}}
end

struct DiscreteOperators{T <: Real, D <: Integer}
    opS::OperatorsConvection{T,D}
    opL::OperatorsConvection{T,D}
    opC_TS::Operators{T,D}
    opC_TL::Operators{T,D}
    opC_pS::Operators{T,D}
    opC_pL::Operators{T,D}
    opC_uS::Operators{T,D}
    opC_uL::Operators{T,D}
    opC_vS::Operators{T,D}
    opC_vL::Operators{T,D}
end

struct Phase{T <: Real} <: MutatingFields
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

struct Forward{T <: Real} <: MutatingFields
    T::Array{T,3}
    u::Array{T,4}
    ux::Array{T,4}
    uy::Array{T,4}
    V::Array{T,3}
    κ::Array{T,4}
    length::Array{T,1}
    t::Array{T,1}
    Cd::Array{T,1}
    Cl::Array{T,1}
end

struct ForwardPhase{T <: Real} <: MutatingFields
    T::Array{T,3}
    p::Array{T,3}
    ϕ::Array{T,3}
    u::Array{T,3}
    v::Array{T,3}
    TD::Array{T,2}
    pD::Array{T,2}
    ucorrD::Array{T,2}
    vcorrD::Array{T,2}
    Vratio::Vector{T}
end

mutable struct Desired{T <: Real} <: AbstractOptimizer
    x_desired::Vector{T}
    u::Matrix{T}
    usave::Array{T, 3}
    TL::Matrix{T}
    TS::Matrix{T}
end

mutable struct Optim_parameters{T <: Real, N <: Integer} <: AbstractOptimizer
    nprobes::N
    ind::Vector{N}
    bc_indices::Vector{CartesianIndex{2}}
    γ::Vector{T}
    p::Vector{Vector{T}}
    TLsave::Vector{Matrix{T}}
    TSsave::Vector{Matrix{T}}
    usave::Vector{Array{T, 3}}
end

mutable struct my_Adjoint{T <: Real} <: MutatingFields
    iso::Array{T, 2}
    u::Array{T,2}
    TS::Array{T,2}
    TL::Array{T,2}
    DTS::Array{T,2}
    DTL::Array{T,2}
    V::Array{T,2}
    κ::Array{T, 2}
end

mutable struct adjoint_phase{T <: Real} <: MutatingFields
    TD::Array{T,2}
    ϕD::Array{T,2}
    pD::Array{T,2}
    uD::Array{T,2}
    vD::Array{T,2}
    ucorrD::Array{T,2}
    vcorrD::Array{T,2}
end

mutable struct adjoint_fields{T <: Real} <: MutatingFields
    u::Array{T,2}
    phS::adjoint_phase{T}
    phL::adjoint_phase{T}
end

mutable struct adjoint_derivatives{T <: Real, D <: Integer} <: MutatingFields
    RheatS_ls::SparseMatrixCSC{T,D} # Derivative of solid phase heat eq. wrt. levelset
    RheatL_ls::SparseMatrixCSC{T,D} # Derivative of liquid phase heat eq. wrt. levelset
    RlsS_ls::SparseMatrixCSC{T,D} # Derivative of Stefan levelset advection eq. wrt. levelset
    RlsS_TS::SparseMatrixCSC{T,D} # Derivative of Stefan levelset advection eq. wrt. solid temperature 
    RlsS_TL::SparseMatrixCSC{T,D} # Derivative of Stefan levelset advection eq. wrt. liquid temperature
    RucorrS_ls0::SparseMatrixCSC{T,D} # Derivative of solid phase prediction step in x wrt. levelset at prev. it.
    RucorrS_ls1::SparseMatrixCSC{T,D} # Derivative of solid phase prediction step in x wrt. levelset
    RucorrL_ls0::SparseMatrixCSC{T,D} # Derivative of liquid phase prediction step in x wrt. levelset at prev. it.
    RucorrL_ls1::SparseMatrixCSC{T,D} # Derivative of liquid phase prediction step in x wrt. levelset
    RvcorrS_ls0::SparseMatrixCSC{T,D} # Derivative of solid phase prediction step in y wrt. levelset at prev. it.
    RvcorrS_ls1::SparseMatrixCSC{T,D} # Derivative of solid phase prediction step in y wrt. levelset
    RvcorrL_ls0::SparseMatrixCSC{T,D} # Derivative of liquid phase prediction step in y wrt. levelset at prev. it.
    RvcorrL_ls1::SparseMatrixCSC{T,D} # Derivative of liquid phase prediction step in y wrt. levelset
    RpS_ls::SparseMatrixCSC{T,D} # Derivative of solid phase Poisson eq. wrt. levelset
    RpL_ls::SparseMatrixCSC{T,D} # Derivative of liquid phase Poisson eq. wrt. levelset
    RuS_ls::SparseMatrixCSC{T,D} # Derivative of solid phase projection step in x wrt. levelset
    RuL_ls::SparseMatrixCSC{T,D} # Derivative of liquid phase projection step in x wrt. levelset
    RvS_ls::SparseMatrixCSC{T,D} # Derivative of solid phase projection step in y wrt. levelset
    RvL_ls::SparseMatrixCSC{T,D} # Derivative of liquid phase projection step in y wrt. levelset
    RlsFS_ls::SparseMatrixCSC{T,D} # Derivative of Free Surface levelset advection eq. wrt. levelset
    RlsFS_ucorrS::SparseMatrixCSC{T,D} # Derivative of Free Surface levelset advection eq. wrt. solid horizontal velocity 
    RlsFS_ucorrL::SparseMatrixCSC{T,D} # Derivative of Free Surface levelset advection eq. wrt. liquid horizontal velocity 
    RlsFS_vcorrS::SparseMatrixCSC{T,D} # Derivative of Free Surface levelset advection eq. wrt. solid vertical velocity
    RlsFS_vcorrL::SparseMatrixCSC{T,D} # Derivative of Free Surface levelset advection eq. wrt. liquid vertical velocity
end

abstract type TemporalIntegration end
struct ForwardEuler <: TemporalIntegration end
struct CrankNicolson <: TemporalIntegration end

abstract type LevelsetDiscretization end
struct WENO5 <: LevelsetDiscretization end
struct ENO2 <: LevelsetDiscretization end

abstract type BoundaryCondition end

@with_kw mutable struct Neumann{T,N} <: BoundaryCondition
   val::N = 0.0
   λ::T = 0.0
end

@with_kw mutable struct Neumann_cl{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
    θe::T = π / 2
end

@with_kw mutable struct Neumann_inh{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

@with_kw mutable struct Dirichlet{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

@with_kw mutable struct Periodic{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

@with_kw mutable struct Robin{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

"""
    Navier{T,N} <: BoundaryCondition

Navier boundary condition applied at the contact lines. Elsewhere, Dirichlet is
applied. The boundary condition is given by:

``u - \\lambda \\dfrac{\\partial u}{\\partial n} = U``
"""
@with_kw mutable struct Navier{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

"""
    Navier_cl{T,N} <: BoundaryCondition

Navier boundary condition applied at the contact lines. Elsewhere, Dirichlet is
applied. The boundary condition is given by:

``u - \\lambda \\dfrac{\\partial u}{\\partial n} = U``
"""
@with_kw mutable struct Navier_cl{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

"""
    GNBC{T,N} <: BoundaryCondition

Generalized Navier boundary condition. It is only applied at the contact lines. Elsewhere,
Dirichlet is applied. The boundary condition is given by:

``u - \\lambda \\dfrac{\\partial u}{\\partial n} = U f ( x / \\epsilon ) \\sigma / \\mu ( \\cos \\theta _ d - \\cos \\theta _ e)``
"""
@with_kw mutable struct GNBC{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
    ϵ::T = 1.0
    σ::T = 1.0
    μ::T = 1.0
    θe::T = π / 2
end

@with_kw mutable struct Boundaries <: NumericalParameters
    left::BoundaryCondition = Neumann()
    right::BoundaryCondition = Neumann()
    bottom::BoundaryCondition = Neumann()
    top::BoundaryCondition = Neumann()
end

@with_kw mutable struct DummyBC{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end
@with_kw mutable struct Stefan{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

@with_kw mutable struct FreeSurface{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
end

@with_kw mutable struct Wall{T,N} <: BoundaryCondition
    val::N = 0.0
    λ::T = 0.0
    θe::T = π / 2
end
