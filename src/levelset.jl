function diamond(a, p, n)
    U = @SVector[0.25(a[p]+a[p-n]+a[p-n-1]+a[p-1]),
    0.25(a[p]+a[p-1]+a[p+n-1]+a[p+n]),
    0.25(a[p]+a[p+n]+a[p+n+1]+a[p+1]),
    0.25(a[p]+a[p+1]+a[p-n+1]+a[p-n])]
    return U
end

function diamond(a, II, nx, ny, per_x, per_y)
    U = @SVector[0.25(a[II]+a[δx⁻(II, nx, per_x)]+a[δx⁻(δy⁻(II, ny, per_y), nx, per_x)]+a[δy⁻(II, ny, per_y)]),
                0.25(a[II]+a[δy⁻(II, ny, per_y)]+a[δy⁻(δx⁺(II, nx, per_x), ny, per_y)]+a[δx⁺(II, nx, per_x)]),
                0.25(a[II]+a[δx⁺(II, nx, per_x)]+a[δx⁺(δy⁺(II, ny, per_y), nx, per_x)]+a[δy⁺(II, ny, per_y)]),
                0.25(a[II]+a[δy⁺(II, ny, per_y)]+a[δy⁺(δx⁻(II, nx, per_x), ny, per_y)]+a[δx⁻(II, nx, per_x)])]
    return U
end

function quadratic_recons(a, U::SArray{Tuple{4},Float64,1,4},
    p::Int64, n::Int64)
    D = @SVector [0.25*(a[p]+a[p-n]+U[1]+U[4]),
    0.25*(a[p]+a[p-1]+U[1]+U[2]),
    0.25*(a[p]+a[p+n]+U[2]+U[3]),
    0.25*(a[p]+a[p+1]+U[3]+U[4])]
    return D
end

function quadratic_recons(a, U::SArray{Tuple{4},Float64,1,4},
    II::CartesianIndex{2}, nx::Int64, ny::Int64, per_x::Bool, per_y::Bool)
    D = @SVector [0.25*(a[II]+a[δx⁻(II, nx, per_x)]+U[1]+U[4]),
                0.25*(a[II]+a[δy⁻(II, ny, per_y)]+U[1]+U[2]),
                0.25*(a[II]+a[δx⁺(II, nx, per_x)]+U[2]+U[3]),
                0.25*(a[II]+a[δy⁺(II, ny, per_y)]+U[3]+U[4])]
    return D
end

function quadratic_recons(D::SArray{Tuple{4},Float64,1,4})
    D_ = @SVector [0.25*(D[1]+D[2]+D[3]+D[4])]
    return D_
end

function normal_grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    D::SArray{Tuple{4},Float64,1,4}, D_::Float64,
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[g[p]*(D_-D[1])/(sqrt(((a[p]-a[p-n])/h)^2 + ((U[1] - U[4])/h)^2)),
    g[p]*(D_-D[2])/(sqrt(((a[p]-a[p-1])/h)^2 + ((U[1] - U[2])/h)^2)),
    g[p]*(D_-D[3])/(sqrt(((a[p]-a[p+n])/h)^2 + ((U[2] - U[3])/h)^2)),
    g[p]*(D_-D[4])/(sqrt(((a[p]-a[p+1])/h)^2 + ((U[3] - U[4])/h)^2))]
    return F
end

function normal_grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[g[p]*(a[p]-a[p-n])/(sqrt(((a[p]-a[p-n])/h)^2 + ((U[1] - U[4])/h)^2)),
    g[p]*(a[p]-a[p-1])/(sqrt(((a[p]-a[p-1])/h)^2 + ((U[1] - U[2])/h)^2)),
    g[p]*(a[p]-a[p+n])/(sqrt(((a[p]-a[p+n])/h)^2 + ((U[2] - U[3])/h)^2)),
    g[p]*(a[p]-a[p+1])/(sqrt(((a[p]-a[p+1])/h)^2 + ((U[3] - U[4])/h)^2))]
    return F
end

function normal_grad(a,
    U::SArray{Tuple{4},Float64,1,4},
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[(a[p]-a[p-n])/(sqrt(((a[p]-a[p-n])/h)^2 + ((U[1] - U[4])/h)^2)),
    (a[p]-a[p-1])/(sqrt(((a[p]-a[p-1])/h)^2 + ((U[1] - U[2])/h)^2)),
    (a[p]-a[p+n])/(sqrt(((a[p]-a[p+n])/h)^2 + ((U[2] - U[3])/h)^2)),
    (a[p]-a[p+1])/(sqrt(((a[p]-a[p+1])/h)^2 + ((U[3] - U[4])/h)^2))]
    return F
end

function grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    D::SArray{Tuple{4},Float64,1,4}, D_::Float64,
    p::Int64, n::Int64)
    F = SA_F64[g[p]*(D_-D[1]),
    g[p]*(D_-D[2]),
    g[p]*(D_-D[3]),
    g[p]*(D_-D[4])]
    return F
end

function grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    p::Int64, n::Int64)
    F = SA_F64[g[p]*(a[p]-a[p-n]),
    g[p]*(a[p]-a[p-1]),
    g[p]*(a[p]-a[p+n]),
    g[p]*(a[p]-a[p+1])]
    return F
end

function grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    II::CartesianIndex{2}, nx::Int64, ny::Int64, per_x::Bool, per_y::Bool)
    F = SA_F64[g[II]*(a[II]-a[δx⁻(II, nx, per_x)]),
            g[II]*(a[II]-a[δy⁻(II, ny, per_y)]),
            g[II]*(a[II]-a[δx⁺(II, nx, per_x)]),
            g[II]*(a[II]-a[δy⁺(II, ny, per_y)])]
    return F
end

function grad_IIOE(a, gx, gy,
    U::SArray{Tuple{4},Float64,1,4},
    D::SArray{Tuple{4},Float64,1,4}, D_::Float64,
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[gx[p]*(D_-D[1]),
    gy[p]*(D_-D[2]),
    gx[p]*(D_-D[3]),
    gy[p]*(D_-D[4])]
    return F
end

function grad_IIOE(a, gx, gy,
    U::SArray{Tuple{4},Float64,1,4},
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[gx[p]*(a[p]-a[p-n]),
    gy[p]*(a[p]-a[p-1]),
    gx[p]*(a[p]-a[p+n]),
    gy[p]*(a[p]-a[p+1])]
    return F
end

function advection(a, gx, gy,
    U::SArray{Tuple{4},Float64,1,4},
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[gx[p],
    gy[p],
    -gx[p],
    -gy[p]] .* h
    return F
end

function advection(gx, gy,
    hx::Matrix{Float64}, hy::Matrix{Float64},
    i::CartesianIndex)
    F = SA_F64[gx[i] * hy[i],
               gy[i] * hx[i],
               -gx[δx⁺(i)] * hy[i],
               -gy[δy⁺(i)] * hx[i]]
    return F
end

function θout(a_out, u, umax, umin, II, τ, nx, ny, mp, per_x, per_y)
    θ_out = zeros(4)

    a = (δx⁻(II, nx, per_x),
         δy⁻(II, ny, per_y), 
         δx⁺(II, nx, per_x),
         δy⁺(II, ny, per_y))
    @inbounds for (i,j) in zip(1:4,a)
        n_out = - sum(sign.(a_out))
        cond = a_out[i] * (u[j] - u[II])
        if abs(cond) < 1e-12
            θ_out[i] = 0.5
        elseif cond >= 1e-12
            θ = mp * (umax - u[II]) / (τ * n_out * cond)
            θ_out[i] = min(0.5, θ)
        else
            θ = mp * (umin - u[II]) / (τ * n_out * cond)
            θ_out[i] = min(0.5, θ)
        end
    end

    return θ_out
end

function θin(θ_out, nx, ny, per_x, per_y, II)
    θ_in = ones(4)

    a = (δx⁻(II, nx, per_x),
         δy⁻(II, ny, per_y), 
         δx⁺(II, nx, per_x),
         δy⁺(II, ny, per_y))
    b = (3,4,1,2)
    @inbounds for (i,j,k) in zip(1:4,a,b)
        θ_in[i] -= θ_out[j,k]
    end
    
    return θ_in
end

function inflow_outflow(F::SArray{Tuple{4},Float64,1,4})
    a_in = @SVector[max(0,F[1]),
    max(0,F[2]),
    max(0,F[3]),
    max(0,F[4])]
    a_ou = @SVector[min(0,F[1]),
    min(0,F[2]),
    min(0,F[3]),
    min(0,F[4])]
    return a_in, a_ou
end

function sumloc(a_in::SArray{Tuple{4},Float64,1,4},
    a_ou::SArray{Tuple{4},Float64,1,4})
    S = @SVector[sum(a_in),
    sum(a_ou)]
    return S
end

function sumloc(a_in::SizedVector{4,Float64,Vector{Float64}},
    a_ou::SizedVector{4,Float64,Vector{Float64}})
    S = @SVector[sum(a_in),
    sum(a_ou)]
    return S
end

"""
    IIOE_normal!(grid, A, B, u, V, CFL, periodic_x, periodic_y)

Advection of the levelset in the normal direction  (Mikula et al. 2014).
"""
function IIOE_normal!(grid, A, B, u, V, CFL, periodic_x, periodic_y)
    @unpack nx, ny, ind = grid
    @unpack inside, all_indices = ind

    if !periodic_x && !periodic_y
        indices = inside
    elseif periodic_x && periodic_y
        indices = all_indices
    elseif periodic_x
        indices = @view all_indices[2:end-1,:]
    else
        indices = @view all_indices[:,2:end-1]
    end

    @inbounds @threads for II in indices
        p = lexicographic(II, ny)
        U = diamond(u, II, nx, ny, periodic_x, periodic_y)
        F = grad(u, V, U, II, nx, ny, periodic_x, periodic_y)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        if S[1] < -S[2]
            D = quadratic_recons(u, U, II, nx, ny, periodic_x, periodic_y)
            D_ = quadratic_recons(D)
            F = 2*grad(u, V, U, D, D_[1], p, ny)
            a_in, a_ou = inflow_outflow(F)
            S = sumloc(a_in, a_ou)
        end
        A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, S, A, B, II, nx, ny, CFL, periodic_x, periodic_y)
    end
end


"""
    IIOE!(grid, grid_u, grid_v, A, B, θ_out, τ, periodic_x, periodic_y)

Advection of the levelset using the basic IIOE scheme (Mikula et al. 2014).
"""
function IIOE!(grid, grid_u, grid_v, A, B, θ_out, τ, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack inside, all_indices = ind

    if !periodic_x && !periodic_y
        indices = inside
    elseif periodic_x && periodic_y
        indices = all_indices
    elseif periodic_x
        indices = @view all_indices[2:end-1,:]
    else
        indices = @view all_indices[:,2:end-1]
    end

    θ_in = 0.5 .* ones(4)

    @inbounds @threads for II in indices
        p = lexicographic(II, ny)
        F = advection(grid_u.V, grid_v.V, dx, dy, II)
        a_in, a_ou = inflow_outflow(F)
        θ_out[II,:] .= 0.5
        S = sumloc(a_in .* θ_in, a_ou .* θ_out[II,:])
        mp = dx[II] * dy[II]
        A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, θ_in, θ_out[II,:], S, A, B, II, nx, ny, mp, τ, periodic_x, periodic_y)
    end
end

"""
    S2IIOE!(grid, grid_u, grid_v, A, B, utmp, u, θ_out, τ, periodic_x, periodic_y)

Advection of the levelset using the S2IIOE scheme (Mikula et al. 2014).
"""
function S2IIOE!(grid, grid_u, grid_v, A, B, utmp, u, θ_out, τ, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack inside, all_indices = ind

    if !periodic_x && !periodic_y
        indices = inside
    elseif periodic_x && periodic_y
        indices = all_indices
    elseif periodic_x
        indices = @view all_indices[2:end-1,:]
    else
        indices = @view all_indices[:,2:end-1]
    end

    θ_in = zeros(4)

    @inbounds @threads for II in indices
        umax = max(u[II], u[δx⁻(II, nx, periodic_x)], u[δy⁻(II, ny, periodic_y)], u[δx⁺(II, nx, periodic_x)], u[δy⁺(II, ny, periodic_y)])
        umin = min(u[II], u[δx⁻(II, nx, periodic_x)], u[δy⁻(II, ny, periodic_y)], u[δx⁺(II, nx, periodic_x)], u[δy⁺(II, ny, periodic_y)])
        if utmp[II] > umax || utmp[II] < umin
            F = advection(grid_u.V, grid_v.V, dx, dy, II)
            _, a_ou = inflow_outflow(F)
            mp = dx[II] * dy[II]
            θ_out[II,:] .= θout(a_ou, u, umax, umin, II, τ, nx, ny, mp, periodic_x, periodic_y)
        end
    end
    @inbounds @threads for II in indices
        p = lexicographic(II, ny)
        umax = max(u[II], u[δx⁻(II, nx, periodic_x)], u[δy⁻(II, ny, periodic_y)], u[δx⁺(II, nx, periodic_x)], u[δy⁺(II, ny, periodic_y)])
        umin = min(u[II], u[δx⁻(II, nx, periodic_x)], u[δy⁻(II, ny, periodic_y)], u[δx⁺(II, nx, periodic_x)], u[δy⁺(II, ny, periodic_y)])
        if utmp[II] > umax || utmp[II] < umin
            F = advection(grid_u.V, grid_v.V, dx, dy, II)
            a_in, a_ou = inflow_outflow(F)
            θ_in .= θin(θ_out, nx, ny, periodic_x, periodic_y, II)
            S = sumloc(a_in .* θ_in, a_ou .* θ_out[II,:])
            mp = dx[II] * dy[II]
            A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, θ_in, θ_out[II,:], S, A, B, II, nx, ny, mp, τ, periodic_x, periodic_y)
        end
    end
end

@inline function fill_matrices2!(a_in::SArray{Tuple{4},Float64,1,4}, a_ou::SArray{Tuple{4},Float64,1,4},
    S, A, B, II, nx, ny, CFL, per_x, per_y)
    p = lexicographic(II, ny)
    a = (lexicographic(δx⁻(II, nx, per_x), ny),
         lexicographic(δy⁻(II, ny, per_y), ny), 
         lexicographic(δx⁺(II, nx, per_x), ny),
         lexicographic(δy⁺(II, ny, per_y), ny))
    @inbounds for (i,j) in zip(1:4,a)
        @inbounds A[p, j] = -CFL*a_in[i]
        @inbounds B[p, j] = CFL*a_ou[i]
    end
    return 2 + CFL*S[1], 2 - CFL*S[2]
end

@inline function fill_matrices2!(a_in::SArray{Tuple{4},Float64,1,4}, a_ou::SArray{Tuple{4},Float64,1,4},
    θ_in, θ_out, S, A, B, II, nx, ny, mp, τ, per_x, per_y)
    p = lexicographic(II, ny)
    a = (lexicographic(δx⁻(II, nx, per_x), ny),
         lexicographic(δy⁻(II, ny, per_y), ny), 
         lexicographic(δx⁺(II, nx, per_x), ny),
         lexicographic(δy⁺(II, ny, per_y), ny))
    @inbounds for (i,j) in zip(1:4,a)
        @inbounds A[p, j] = -τ * θ_in[i] * a_in[i] / mp
        @inbounds B[p, j] = τ * θ_out[i] * a_ou[i] / mp
    end
    return 1 + τ*S[1]/mp, 1 - τ*S[2]/mp
end

@inline central_differences(u, II, dx, dy, nx, ny, per_x, per_y) =
    @SVector[minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_l(u, II, dx, dy, nx, ny, per_x, per_y, dxx) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_b(u, II, dx, dy, nx, ny, per_x, per_y, dyy) =
    @SVector[minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_r(u, II, dx, dy, nx, ny, per_x, per_y, dxx) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.)),
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_t(u, II, dx, dy, nx, ny, per_x, per_y, dyy) =
    @SVector[minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(Dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), Dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_bl(u, II, dx, dy, nx, ny, per_x, per_y, dxx, dyy) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_br(u, II, dx, dy, nx, ny, per_x, per_y, dxx, dyy) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁺(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁺(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_tr(u, II, dx, dy, nx, ny, per_x, per_y, dxx, dyy) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁻(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁻(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]
@inline central_differences_tl(u, II, dx, dy, nx, ny, per_x, per_y, dxx, dyy) =
    @SVector[minmod(dxx(u, II, dx, nx, per_x), ifelse(in_bounds(δx⁺(II, nx, per_x)[2], nx, per_x), dxx(u, δx⁺(II, nx, per_x), dx, nx, per_x), 0.)),
            minmod(dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))
]

@inline finite_difference_eno(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[2],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[3],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[4]
]
@inline finite_difference_eno_l(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[2],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[3]
]
@inline finite_difference_eno_b(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[2],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[3]
]
@inline finite_difference_eno_r(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[-2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[1],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[2],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[3]
]
@inline finite_difference_eno_t(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[2],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[3]
]
@inline finite_difference_eno_bl(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[2]
]
@inline finite_difference_eno_br(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[-2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[1],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[2]
]
@inline finite_difference_eno_tr(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[-2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[2]
]
@inline finite_difference_eno_tl(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[2]
]

@inline Godunov(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + max(⁻(a[3])^2, ⁺(a[4])^2)),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + max(⁺(a[3])^2, ⁻(a[4])^2)))
@inline Godunov_l(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁻(a[1])^2 + max(⁻(a[2])^2, ⁺(a[3])^2)),
    sqrt(⁺(a[1])^2 + max(⁺(a[2])^2, ⁻(a[3])^2)))
@inline Godunov_b(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + ⁻(a[3])^2),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + ⁺(a[3])^2))
@inline Godunov_r(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁺(a[1])^2 + max(⁻(a[2])^2, ⁺(a[3])^2)),
    sqrt(⁻(a[1])^2 + max(⁺(a[2])^2, ⁻(a[3])^2)))
@inline Godunov_t(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + ⁺(a[3])^2),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + ⁻(a[3])^2))
@inline Godunov_bl(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁻(a[1])^2 + ⁻(a[2])^2),
    sqrt(⁺(a[1])^2 + ⁺(a[2])^2))
@inline Godunov_br(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁺(a[1])^2 + ⁻(a[2])^2),
    sqrt(⁻(a[1])^2 + ⁺(a[2])^2))
@inline Godunov_tr(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁺(a[1])^2 + ⁺(a[2])^2),
    sqrt(⁻(a[1])^2 + ⁻(a[2])^2))
@inline Godunov_tl(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(⁻(a[1])^2 + ⁺(a[2])^2),
    sqrt(⁺(a[1])^2 + ⁻(a[2])^2))

@inline root_extraction(u, uxx, D, II, JJ, h, eps) =
    ifelse(abs(uxx) > eps,
    h*(0.5 + ((u[II] - u[JJ] - mysign(u[II] - u[JJ]) * sqrt(D))/uxx)),
    h*(u[II]/(u[II]-u[JJ])))

"""
    Φweno(a)

Intermediate function for the 5th-order WENO as defined in (Jiang and Pen, 2000).
"""
Φweno(a) = (
    1.0 / 3.0 * ω0(a) * (a[1] - 2.0 * a[2] + a[3]) + 
    1.0 / 6.0 * (ω2(a) - 0.5) * (a[2] - 2.0 * a[3] + a[4])
)

"""
    Φeno(a)

Intermediate function for the 3rd-order ENO as defined in (Jiang and Pen, 2000).
"""
function Φeno(a)
    u1 = - 1.0 / 12.0 * (a[2] - 2.0 * a[3] + a[4])

    if abs(a[2]) < abs(a[3]) && abs(a[1] - a[2]) < abs(a[2] - a[3])
        return u1 + 1.0 / 3.0 * (a[1] - 2.0 * a[2] + a[3])
    elseif abs(a[2]) > abs(a[3]) && abs(a[2] - a[3]) > abs(a[3] - a[4])
        return -u1
    else
        return u1
    end
end

"""
    finite_difference_weno5(u, II, nx, ny, dx, dy, per_x, per_y)

5th order WENO scheme as defined in (Jiang and Pen, 2000).

Replace Φweno by Φeno to use the 3rd-order ENO scheme.
"""
function finite_difference_weno5(grid, u, II, nx, ny, dx, dy, per_x, per_y)
    a = Vector{Float64}(undef, 4)

    II_1 = δx⁻(II, nx, per_x)
    II_2 = δx⁻(II_1, nx, per_x)
    II_3 = δx⁻(II_2, nx, per_x)
    II1 = δx⁺(II, nx, per_x)
    II2 = δx⁺(II1, nx, per_x)
    II3 = δx⁺(II2, nx, per_x)

    if II[2] > 4
        d_3 = ∇x⁺(u, II_3, nx, dx, per_x)
        d_2 = ∇x⁺(u, II_2, nx, dx, per_x)
        d_1 = ∇x⁺(u, II_1, nx, dx, per_x)
        if within_bounds(II1, grid)
            d = ∇x⁺(u, II, nx, dx, per_x)
        else
            d = d_1
        end
        if within_bounds(II2, grid)
            d1 = ∇x⁺(u, II1, nx, dx, per_x)
        else
            d1 = d
        end
        if within_bounds(II3, grid)
            d2 = ∇x⁺(u, II2, nx, dx, per_x)
        else
            d2 = d1
        end
    else
        d2 = ∇x⁺(u, II2, nx, dx, per_x)
        d1 = ∇x⁺(u, II1, nx, dx, per_x)
        d = ∇x⁺(u, II, nx, dx, per_x)
        if within_bounds(II_1, grid)
            d_1 = ∇x⁺(u, II_1, nx, dx, per_x)
        else
            d_1 = d
        end
        if within_bounds(II_2, grid)
            d_2 = ∇x⁺(u, II_2, nx, dx, per_x)
        else
            d_2 = d_1
        end
        if within_bounds(II_3, grid)
            d_3 = ∇x⁺(u, II_3, nx, dx, per_x)
        else
            d_3 = d_2
        end
    end

    cnst = 1.0 / 12.0 * (-d_2 + 7.0 * d_1 + 7.0 * d - d1)

    @inbounds a[1] = d2 - d1
    @inbounds a[2] = d1 - d
    @inbounds a[3] = d - d_1
    @inbounds a[4] = d_1 - d_2
    fdx⁺ = cnst + Φweno(a)

    @inbounds a[1] = d_2 - d_3
    @inbounds a[2] = d_1 - d_2
    @inbounds a[3] = d - d_1
    @inbounds a[4] = d1 - d
    fdx⁻ = cnst - Φweno(a)

    II_1 = δy⁻(II, ny, per_y)
    II_2 = δy⁻(II_1, ny, per_y)
    II_3 = δy⁻(II_2, ny, per_y)
    II1 = δy⁺(II, ny, per_y)
    II2 = δy⁺(II1, ny, per_y)
    II3 = δy⁺(II2, ny, per_y)

    if II[1] > 4
        d_3 = ∇y⁺(u, II_3, ny, dy, per_y)
        d_2 = ∇y⁺(u, II_2, ny, dy, per_y)
        d_1 = ∇y⁺(u, II_1, ny, dy, per_y)
        if within_bounds(II1, grid)
            d = ∇y⁺(u, II, ny, dy, per_y)
        else
            d = d_1
        end
        if within_bounds(II2, grid)
            d1 = ∇y⁺(u, II1, ny, dy, per_y)
        else
            d1 = d
        end
        if within_bounds(II3, grid)
            d2 = ∇y⁺(u, II2, ny, dy, per_y)
        else
            d2 = d1
        end
    else
        d2 = ∇y⁺(u, II2, ny, dy, per_y)
        d1 = ∇y⁺(u, II1, ny, dy, per_y)
        d = ∇y⁺(u, II, ny, dy, per_y)
        if within_bounds(II_1, grid)
            d_1 = ∇y⁺(u, II_1, ny, dy, per_y)
        else
            d_1 = d
        end
        if within_bounds(II_2, grid)
            d_2 = ∇y⁺(u, II_2, ny, dy, per_y)
        else
            d_2 = d_1
        end
        if within_bounds(II_3, grid)
            d_3 = ∇y⁺(u, II_3, ny, dy, per_y)
        else
            d_3 = d_2
        end
    end

    cnst = 1.0 / 12.0 * (-d_2 + 7.0 * d_1 + 7.0 * d - d1)

    @inbounds a[1] = d2 - d1
    @inbounds a[2] = d1 - d
    @inbounds a[3] = d - d_1
    @inbounds a[4] = d_1 - d_2
    fdy⁺ = cnst + Φweno(a)

    @inbounds a[1] = d_2 - d_3
    @inbounds a[2] = d_1 - d_2
    @inbounds a[3] = d - d_1
    @inbounds a[4] = d1 - d
    fdy⁻ = cnst - Φweno(a)

    return @SVector[fdx⁺, fdx⁻, fdy⁺, fdy⁻]
end

s0(a) = 13.0 * (a[1] - a[2])^2 + 3.0 * (a[1] - 3.0 * a[2])^2
s1(a) = 13.0 * (a[2] - a[3])^2 + 3.0 * (a[2] + a[3])^2
s2(a) = 13.0 * (a[3] - a[4])^2 + 3.0 * (3.0 * a[3] - a[4])^2

α0(a) = 1.0 / (1e-6 + s0(a))^2
α1(a) = 6.0 / (1e-6 + s1(a))^2
α2(a) = 3.0 / (1e-6 + s2(a))^2

ω0(a) = α0(a) / (α0(a) + α1(a) + α2(a))
ω2(a) = α2(a) / (α0(a) + α1(a) + α2(a))

"""
    δ0(u, II, ϵ, nx, ny, per_x, per_y)

δ function for (Russo and Smereka, 2000) subcell fix.

Not used in the code as (Hartmann et al., 2010) works better.

TODO: only works with periodic BCs for now!!
"""
function δ0(u, II, ϵ, nx, ny, per_x, per_y)
    return maximum([
        abs(u[δx⁺(II,nx,per_x)] - u[II]), abs(u[δx⁻(II,nx,per_x)] - u[II]),
        abs(u[δy⁺(II,ny,per_y)] - u[II]), abs(u[δy⁻(II,ny,per_y)] - u[II]),
        sqrt((u[δx⁺(II,nx,per_x)] - u[δx⁻(II,nx,per_x)])^2 + 
             (u[δy⁺(II,ny,per_y)] - u[δy⁻(II,ny,per_y)])^2)/ 2.0, ϵ]
    )
end

"""
    reinit_min(scheme, grid, u, u0, periodic_x, periodic_y)

Run one reinitilization iteration using the subcell fix of (Min, 2010).

`scheme` can be either the 5th-order WENO scheme `weno5` (Jiang and Pen, 2000) or
the 2nd-order ENO `eno2`.
"""
function reinit_min(scheme, grid, u, u0, indices, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    local cfl = 0.45
    u1 = copy(u0)

    @inbounds @threads for II in indices
        sign_u0 = sign(u0[II])
        if (II in inside || 
            ((II in b_left[1][2:end-1] || II in b_right[1][2:end-1]) && periodic_x) || 
            ((II in b_bottom[1][2:end-1] || II in b_top[1][2:end-1]) && periodic_y) ||
            ((II == b_left[1][1] || II == b_left[1][end] || II == b_right[1][1] || II == b_right[1][end]) && periodic_x && periodic_y)
            )
            II_0 = II
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface
            a, b, c = (1, 2, 3, 4), (-1, 1, -1, 1), (2, 2, 1, 1)
            f, d, e = (Dxx, Dxx, Dyy, Dyy), (nx, nx, ny, ny), (periodic_x, periodic_x, periodic_y, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y), δy⁻(II, ny, periodic_y))
            shift = central_differences(u, II, dx, dy, nx, ny, periodic_x, periodic_y)
            eno = finite_difference_eno(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
            end
            god = Godunov
        elseif II == b_bottom[1][1]
            II_0 = δy⁺(δx⁺(II))
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_bl
            a, b, c = (1, 2), (-1, -1), (2, 1)
            f, d, e = (Dxx_l, Dyy_b), (nx, ny), (periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            shift = central_differences_bl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_b)
            eno = finite_difference_eno_bl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_bl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II == b_bottom[1][end]
            II_0 = δy⁺(δx⁻(II))
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_br
            a, b, c = (1, 2), (1, -1), (2, 1)
            f, d, e = (Dxx_r, Dyy_b), (nx, ny), (periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            neighbour = (δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            shift = central_differences_br(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_b)
            eno = finite_difference_eno_br(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_br
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II == b_top[1][1]
            II_0 = δy⁻(δx⁺(II))
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_tl
            a, b, c = (1, 2), (-1, 1), (2, 1)
            f, d, e = (Dxx_l, Dyy_t), (nx, ny), (periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δy⁻(II, ny, periodic_y))
            shift = central_differences_tl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_t)
            eno = finite_difference_eno_tl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_tl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II == b_top[1][end]
            II_0 = δy⁻(δx⁻(II))
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_tr
            a, b, c = (1, 2), (1, 1), (2, 1)
            f, d, e = (Dxx_r, Dyy_t), (nx, ny), (periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y))
            shift = central_differences_tr(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_t)
            eno = finite_difference_eno_tr(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_tr
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II in b_left[1]
            II_0 = δx⁺(II)
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_l
            a, b, c = (1, 2, 3), (-1, -1, 1), (2, 1, 1)
            f, d, e = (Dxx_l, Dyy, Dyy), (nx, ny, ny), (periodic_x, periodic_y, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]),
                 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y), δy⁻(II, ny, periodic_y))
            shift = central_differences_l(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l)
            eno = finite_difference_eno_l(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_l
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II in b_bottom[1]
            II_0 = δy⁺(II)
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_b
            a, b, c = (1, 2, 3), (-1, 1, -1), (2, 2, 1)
            f, d, e = (Dxx, Dxx, Dyy_b), (nx, nx, ny), (periodic_x, periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            shift = central_differences_b(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_b)
            eno = finite_difference_eno_b(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_b
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II in b_right[1]
            II_0 = δx⁻(II)
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_r
            a, b, c = (1, 2, 3), (1, -1, 1), (2, 1, 1)
            f, d, e = (Dxx_r, Dyy, Dyy), (nx, ny, ny), (periodic_x, periodic_y, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y), δy⁻(II, ny, periodic_y))
            shift = central_differences_r(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r)
            eno = finite_difference_eno_r(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_r
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        elseif II in b_top[1]
            II_0 = δy⁻(II)
            st = static_stencil(u, II_0, nx, ny, periodic_x, periodic_y)
            near_interface = is_near_interface_t
            a, b, c = (1, 2, 3), (-1, 1, 1), (2, 2, 1)
            f, d, e = (Dxx, Dxx, Dyy_t), (nx, nx, ny), (periodic_x, periodic_x, periodic_y)
            g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            neighbour = (δx⁺(II, nx, periodic_x), δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y))
            shift = central_differences_t(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_t)
            eno = finite_difference_eno_t(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            if is_eno(scheme)
                diffs = copy(eno)
                god = Godunov_t
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
        end

        if near_interface(u0, II, nx, ny, periodic_x, periodic_y)
            eno_interface = convert(Vector{Float64}, eno)
            h_ = 1e30
            for (JJ, i, j, k) in zip(neighbour, a, c, g)
                if u0[II] * u0[JJ] < 0
                    uxx = minmod(
                        f[i](u, II, d[i], e[i]), 
                        ifelse(
                            in_bounds(
                                JJ[j], 
                                (j == 1 ? ny : nx), 
                                (j == 1 ? periodic_y : periodic_x)
                            ), 
                            f[i](u, JJ, d[i], e[i]), 
                            0.0
                        )
                    )
                    D = (uxx / 2.0 - u0[II] - u0[JJ])^2 - 4.0 * u0[II] * u0[JJ]
                    Δx = root_extraction(u0, uxx, D, II, JJ, k, 1e-10)
                    if Δx < h_ h_ = Δx end
                    eno_interface[i] = b[i] * (u[II] / Δx + (Δx / 2.0) * shift[i])
                end
            end
            gdv = god(sign_u0, eno_interface)
            u1[II] = u[II] - cfl * h_ * sign_u0 * (gdv - 1.0)
        else
            gdv = god(sign_u0, diffs)
            u1[II] = u[II] - cfl * h_ * sign_u0 * (gdv - 1.0)
        end
    end

    return u1
end

"""
    reinit_rs(scheme, grid, u, u0, periodic_x, periodic_y)

Run one reinitilization iteration using the subcell fix of (Russo and Smereka, 2000).

`scheme` can be either the 5th-order WENO scheme `weno5` (Jiang and Pen, 2000) or
the 2nd-order ENO `eno2`.
"""
function reinit_rs(scheme, grid, u, u0, indices, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    local cfl = 0.45
    u1 = copy(u0)

    @inbounds @threads for II in indices
        sign_u0 = sign(u0[II])
        if (II in inside || 
            ((II in b_left[1][2:end-1] || II in b_right[1][2:end-1]) && periodic_x) || 
            ((II in b_bottom[1][2:end-1] || II in b_top[1][2:end-1]) && periodic_y) ||
            ((II == b_left[1][1] || II == b_left[1][end] || II == b_right[1][1] || II == b_right[1][end]) && periodic_x && periodic_y)
            )
            if is_eno(scheme)
                shift = central_differences(u, II, dx, dy, nx, ny, periodic_x, periodic_y)
                diffs = finite_difference_eno(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            god = Godunov
            near_interface = is_near_interface
        elseif II == b_bottom[1][1]
            if is_eno(scheme)
                shift = central_differences_bl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_b)
                diffs = finite_difference_eno_bl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_bl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_bl
        elseif II == b_bottom[1][end]
            if is_eno(scheme)
                shift = central_differences_br(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_b)
                diffs = finite_difference_eno_br(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_br
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_br
        elseif II == b_top[1][1]
            if is_eno(scheme)
                shift = central_differences_tl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_t)
                diffs = finite_difference_eno_tl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_tl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_tl
        elseif II == b_top[1][end]
            if is_eno(scheme)
                shift = central_differences_tr(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_t)
                diffs = finite_difference_eno_tr(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_tr
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_tr
        elseif II in b_left[1]
            if is_eno(scheme)
                shift = central_differences_l(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l)
                diffs = finite_difference_eno_l(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_l
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_l
        elseif II in b_bottom[1]
            if is_eno(scheme)
                shift = central_differences_b(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_b)
                diffs = finite_difference_eno_b(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_b
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_b
        elseif II in b_right[1]
            if is_eno(scheme)
                shift = central_differences_r(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r)
                diffs = finite_difference_eno_r(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_r
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_r
        elseif II in b_top[1]
            if is_eno(scheme)
                shift = central_differences_t(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_t)
                diffs = finite_difference_eno_t(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_t
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_t
        end

        if near_interface(u0, II, nx, ny, periodic_x, periodic_y)
            d = h_ * u0[II] / δ0(u0, II, h_, nx, ny, periodic_x, periodic_y)
            u1[II] = u[II] - cfl * (sign_u0 * abs(u[II]) - d)
        else
            gdv = god(sign_u0, diffs)
            u1[II] = u[II] - cfl * h_ * sign_u0 * (gdv - 1.0)
        end
    end

    return u1
end

"""
    reinit_hartmann(scheme, grid, u, u0, periodic_x, periodic_y)

Run one reinitilization iteration following (Hartmann et al., 2010).

`scheme` can be either the 5th-order WENO scheme `weno5` (Jiang and Pen, 2000) or
the 2nd-order ENO `eno2`.
"""
function reinit_hartmann(scheme, grid, u, u0, indices, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    local cfl = 0.45
    u1 = copy(u0)
    utmp = copy(u0)

    @inbounds @threads for II in indices
        sign_u0 = sign(u0[II])
        if (II in inside || 
            ((II in b_left[1][2:end-1] || II in b_right[1][2:end-1]) && periodic_x) || 
            ((II in b_bottom[1][2:end-1] || II in b_top[1][2:end-1]) && periodic_y) ||
            ((II == b_left[1][1] || II == b_left[1][end] || II == b_right[1][1] || II == b_right[1][end]) && periodic_x && periodic_y)
            )
            if is_eno(scheme)
                shift = central_differences(u, II, dx, dy, nx, ny, periodic_x, periodic_y)
                diffs = finite_difference_eno(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
            end
            idx = (
                δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y), 
                δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y)
            )
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            god = Godunov
            near_interface = is_near_interface
        elseif II == b_bottom[1][1]
            if is_eno(scheme)
                shift = central_differences_bl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_b)
                diffs = finite_difference_eno_bl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_bl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = [δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y)]
            h_ = min(0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]))
            near_interface = is_near_interface_bl
        elseif II == b_bottom[1][end]
            if is_eno(scheme)
                shift = central_differences_br(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_b)
                diffs = finite_difference_eno_br(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_br
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx =(δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            h_ = min(0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]))
            near_interface = is_near_interface_br
        elseif II == b_top[1][end]
            if is_eno(scheme)
                shift = central_differences_tr(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r, Dyy_t)
                diffs = finite_difference_eno_tr(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_tr
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y))
            h_ = min(0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            near_interface = is_near_interface_tr
        elseif II == b_top[1][1]
            if is_eno(scheme)
                shift = central_differences_tl(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l, Dyy_t)
                diffs = finite_difference_eno_tl(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_tl
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δy⁻(II, ny, periodic_y), δx⁺(II, nx, periodic_x))
            h_ = min(0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            near_interface = is_near_interface_tl
        elseif II in b_left[1]
            if is_eno(scheme)
                shift = central_differences_l(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_l)
                diffs = finite_difference_eno_l(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_l
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δy⁻(II, ny, periodic_y), δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_l
        elseif II in b_bottom[1]
            if is_eno(scheme)
                shift = central_differences_b(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_b)
                diffs = finite_difference_eno_b(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_b
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δx⁻(II, nx, periodic_x), δx⁺(II, nx, periodic_x), δy⁺(II, ny, periodic_y))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_b
        elseif II in b_right[1]
            if is_eno(scheme)
                shift = central_differences_r(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dxx_r)
                diffs = finite_difference_eno_r(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_r
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y), δy⁺(II, ny, periodic_y))
            h_ = min(
                0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_r
        elseif II in b_top[1]
            if is_eno(scheme)
                shift = central_differences_t(u, II, dx, dy, nx, ny, periodic_x, periodic_y, Dyy_t)
                diffs = finite_difference_eno_t(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)
                god = Godunov_t
            else
                diffs = finite_difference_weno5(grid, u, II, nx, ny, dx, dy, periodic_x, periodic_y)
                god = Godunov
            end
            idx = (δx⁻(II, nx, periodic_x), δy⁻(II, ny, periodic_y), δx⁺(II, nx, periodic_x))
            h_ = min(
                0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)])
            )
            near_interface = is_near_interface_t
        end
        
        ssign_u0 = mysign(u0[II], h_)

        if near_interface(u0, II, nx, ny, periodic_x, periodic_y)
            neighbours = Vector{CartesianIndex}(undef, 4)
            m = 0
            for JJ in idx
                if u0[JJ] * u0[II] < 0
                    m += 1
                    neighbours[m] = JJ
                end
            end

            gdv = god(sign_u0, diffs)
            utmp[II] = u[II] - cfl * h_ * ssign_u0 * (gdv - 1.0)

            fix = true
            for i in 1:m
                if u[neighbours[i]] * u[II] > 0
                    fix = false
                end
            end

            if fix
                @inbounds r = u0[II] / sum(u0[neighbours[1:m]])
                @inbounds F = 1.0 / h_ * (r * sum(utmp[neighbours[1:m]]) - utmp[II])
                u1[II] = utmp[II] + cfl * h_ * 0.5 * F
            else
                u1[II] = utmp[II]
            end
        else
            gdv = god(sign_u0, diffs)
            u1[II] = u[II] - cfl * h_ * ssign_u0 * (gdv - 1.0)
        end
    end

    return u1
end

"""
    FE_reinit!(scheme, grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)

Reinitializes a levelset `u` using a Forward Euler integration scheme.

`nb_reinit` iterations are performed to reach the stationary state. The actual work is done
in the `reinit_hartmann` function.
"""
function FE_reinit!(scheme, grid, ind, u, nb_reinit, periodic_x, periodic_y, BC)
    @unpack nx, ny, dx, dy = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    u0 = copy(u)

    indices = vcat(
        vec(inside), 
        !is_neumann_cl(BC.left) ? b_left[1][2:end-1] : [], 
        !is_neumann_cl(BC.bottom) ? b_bottom[1] : [],
        !is_neumann_cl(BC.right) ? b_right[1][2:end-1] : [],
        !is_neumann_cl(BC.top) ? b_top[1] : []
    )

    for nb = 1:nb_reinit
        u .= reinit_hartmann(scheme, grid, u, u0, indices, periodic_x, periodic_y)
    end

    return nothing
end

"""
    RK2_reinit!(scheme, grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)

Reinitializes a levelset using a 2nd-order Runge-Kutta integration scheme.

`nb_reinit` iterations are performed to reach the stationary state. The actual work is done
in the `reinit_hartmann` function.
"""
function RK2_reinit!(scheme, grid, ind, u, nb_reinit, periodic_x, periodic_y, BC)
    @unpack nx, ny, dx, dy = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind

    u0 = copy(u)
    tmp1 = copy(u)
    tmp2 = copy(u)

    indices = vcat(
        vec(inside), 
        !is_neumann_cl(BC.left) ? b_left[1][2:end-1] : [], 
        !is_neumann_cl(BC.bottom) ? b_bottom[1] : [],
        !is_neumann_cl(BC.right) ? b_right[1][2:end-1] : [],
        !is_neumann_cl(BC.top) ? b_top[1] : []
    )

    for nb in 1:nb_reinit
        tmp1 .= reinit_hartmann(scheme, grid, u, u0, indices, periodic_x, periodic_y)
        tmp2 .= reinit_hartmann(scheme, grid, tmp1, u0, indices, periodic_x, periodic_y)
        u .= 0.5 .* (u .+ tmp2)
    end

    return nothing
end


"""
    rg(grid, periodic_x, periodic_y)

Compute the deviation of the levelset from a distance function following (Luddens et al., 2015).

Computes rg(∇ϕ) = (|∇ϕ| - 1) _ {L ^ 1}.
"""
function rg(grid, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind, u = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind

    tmp = zeros(grid)

    @inbounds @threads for II in all_indices
        if II in inside || ((II in b_left[1] || II in b_right[1]) && periodic_x) || ((II in b_bottom[1] || II in b_top[1]) && periodic_y)
            hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
            hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
            gx = c∇x(u, II, hx, nx, periodic_x)
            gy = c∇y(u, II, hy, ny, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II == b_bottom[1][1]
            gx = ∇x⁺(u, II, nx, dx, periodic_x)
            gy = ∇y⁺(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II == b_bottom[1][end]
            gx = ∇x⁻(u, II, nx, dx, periodic_x)
            gy = ∇y⁺(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II == b_top[1][end]
            gx = ∇x⁻(u, II, nx, dx, periodic_x)
            gy = ∇y⁻(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II == b_top[1][1]
            gx = ∇x⁺(u, II, nx, dx, periodic_x)
            gy = ∇y⁻(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II in b_left[1]
            hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
            gx = ∇x⁺(u, II, nx, dx, periodic_x)
            gy = c∇y(u, II, hy, ny, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II in b_bottom[1]
            hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
            gx = c∇x(u, II, hx, nx, periodic_x)
            gy = ∇y⁺(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II in b_right[1]
            hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
            gx = ∇x⁻(u, II, nx, dx, periodic_x)
            gy = c∇y(u, II, hy, ny, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        elseif II in b_top[1]
            hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
            gx = c∇x(u, II, hx, nx, periodic_x)
            gy = ∇y⁻(u, II, ny, dy, periodic_y)
            tmp[II] = sqrt(gx^2 + gy^2) - 1.0
        end
    end

    return sum(abs.(tmp))
end

function velocity_extension!(grid, V, inside, NB)
    @unpack dx, dy, u = grid

    local cfl = 0.45
    local Vt = similar(V)
    for j = 1:NB
        Vt .= V
        @inbounds @threads for II in inside
            sx = mysign(u[II], dx[II])
            sy = mysign(u[II], dy[II])
            nx = mysign(c∇x(u, II), c∇y(u, II))
            ny = mysign(c∇y(u, II), c∇x(u, II))
            V[II] = Vt[II] - cfl*(⁺(sx*nx)*(-∇x⁻(Vt, II)) +
            ⁻(sx*nx)*(∇x⁺(Vt, II)) +
            ⁺(sy*ny)*(-∇y⁻(Vt, II)) +
            ⁻(sy*ny)*(∇y⁺(Vt, II)))
        end
    end
end

function velocity_extension!(grid, u, V, indices, NB, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid

    local cfl = 0.45
    local Vt = similar(V)

    τ = cfl * max(dx..., dy...)

    for j = 1:NB
        Vt .= V
        @inbounds @threads for II in indices
            cfl_x = τ / dx[II]
            cfl_y = τ / dy[II]
            sx = mysign(u[II], dx[II])
            sy = mysign(u[II], dy[II])
            hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
            hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
            nnx = mysign(c∇x(u, II, hx, nx, periodic_x), c∇y(u, II, hy, ny, periodic_y))
            nny = mysign(c∇y(u, II, hy, ny, periodic_y), c∇x(u, II, hx, nx, periodic_x))
            V[II] = Vt[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(Vt, II, nx, periodic_x)) +
                                      ⁻(sx*nnx)*(∇x⁺(Vt, II, nx, periodic_x))) -
                             cfl_y * (⁺(sy*nny)*(-∇y⁻(Vt, II, ny, periodic_y)) +
                                      ⁻(sy*nny)*(∇y⁺(Vt, II, ny, periodic_y)))
        end
    end
end

function field_extension!(grid, u, f, indices, NB, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, ind = grid
    @unpack all_indices, b_left, b_bottom, b_right, b_top = ind

    local cfl = 0.45
    local ft = similar(f)

    τ = cfl * max(dx..., dy...)

    if periodic_x && periodic_y
        _indices = indices
    elseif !periodic_x && periodic_y
        _indices = all_indices[:,2:end-1]
    elseif periodic_x && !periodic_y
        _indices = all_indices[2:end-1,:]
    else
        _indices = all_indices[2:end-1,2:end-1]
    end

    for j = 1:NB
        ft .= f

        if !periodic_x
            @inbounds @threads for II in b_left[1][2:end-1]
                cfl_x = τ / dx[II]
                cfl_y = τ / dy[II]
                sx = mysign(u[II], dx[II])
                sy = mysign(u[II], dy[II])
                II_0 = δx⁺(II, nx, periodic_x)
                hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
                nnx = mysign(∇x⁺(u, II, nx, dx, periodic_x), c∇y(u, II, hy, ny, periodic_y))
                nny = mysign(c∇y(u, II, hy, ny, periodic_y), ∇x⁺(u, II, nx, dx, periodic_x))
                f[II] = ft[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(ft, II_0, nx, periodic_x)) +
                                        ⁻(sx*nnx)*(∇x⁺(ft, II, nx, periodic_x))) -
                                cfl_y * (⁺(sy*nny)*(-∇y⁻(ft, II, ny, periodic_y)) +
                                        ⁻(sy*nny)*(∇y⁺(ft, II, ny, periodic_y)))
            end
            @inbounds @threads for II in b_right[1][2:end-1]
                cfl_x = τ / dx[II]
                cfl_y = τ / dy[II]
                sx = mysign(u[II], dx[II])
                sy = mysign(u[II], dy[II])
                II_0 = δx⁻(II, nx, periodic_x)
                hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
                nnx = mysign(∇x⁻(u, II, nx, dx, periodic_x), c∇y(u, II, hy, ny, periodic_y))
                nny = mysign(c∇y(u, II, hy, ny, periodic_y), ∇x⁻(u, II, nx, dx, periodic_x))
                f[II] = ft[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(ft, II, nx, periodic_x)) +
                                        ⁻(sx*nnx)*(∇x⁺(ft, II_0, nx, periodic_x))) -
                                cfl_y * (⁺(sy*nny)*(-∇y⁻(ft, II, ny, periodic_y)) +
                                        ⁻(sy*nny)*(∇y⁺(ft, II, ny, periodic_y)))
            end
        end
        if !periodic_y
            @inbounds @threads for II in b_bottom[1][2:end-1]
                cfl_x = τ / dx[II]
                cfl_y = τ / dy[II]
                sx = mysign(u[II], dx[II])
                sy = mysign(u[II], dy[II])
                II_0 = δy⁺(II, ny, periodic_y)
                hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
                nnx = mysign(c∇x(u, II, hx, nx, periodic_x), ∇y⁺(u, II, ny, dy, periodic_y))
                nny = mysign(∇y⁺(u, II, ny, dy, periodic_y), c∇x(u, II, hx, nx, periodic_x))
                f[II] = ft[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(ft, II, nx, periodic_x)) +
                                        ⁻(sx*nnx)*(∇x⁺(ft, II, nx, periodic_x))) -
                                cfl_y * (⁺(sy*nny)*(-∇y⁻(ft, II_0, ny, periodic_y)) +
                                        ⁻(sy*nny)*(∇y⁺(ft, II, ny, periodic_y)))
            end
            @inbounds @threads for II in b_top[1][2:end-1]
                cfl_x = τ / dx[II]
                cfl_y = τ / dy[II]
                sx = mysign(u[II], dx[II])
                sy = mysign(u[II], dy[II])
                II_0 = δy⁻(II, ny, periodic_y)
                hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
                nnx = mysign(c∇x(u, II, hx, nx, periodic_x), ∇y⁻(u, II, ny, dy, periodic_y))
                nny = mysign(∇y⁻(u, II, ny, dy, periodic_y), c∇x(u, II, hx, nx, periodic_x))
                f[II] = ft[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(ft, II, nx, periodic_x)) +
                                        ⁻(sx*nnx)*(∇x⁺(ft, II, nx, periodic_x))) -
                                cfl_y * (⁺(sy*nny)*(-∇y⁻(ft, II, ny, periodic_y)) +
                                        ⁻(sy*nny)*(∇y⁺(ft, II_0, ny, periodic_y)))
            end
        end

        @inbounds @threads for II in _indices
            cfl_x = τ / dx[II]
            cfl_y = τ / dy[II]
            sx = mysign(u[II], dx[II])
            sy = mysign(u[II], dy[II])
            hx = dx[II] + dx[δx⁺(II, nx, periodic_x)] / 2.0 + dx[δx⁻(II, nx, periodic_x)] / 2.0
            hy = dy[II] + dy[δy⁺(II, ny, periodic_y)] / 2.0 + dy[δy⁻(II, ny, periodic_y)] / 2.0
            nnx = mysign(c∇x(u, II, hx, nx, periodic_x), c∇y(u, II, hy, ny, periodic_y))
            nny = mysign(c∇y(u, II, hy, ny, periodic_y), c∇x(u, II, hx, nx, periodic_x))
            f[II] = ft[II] - cfl_x * (⁺(sx*nnx)*(-∇x⁻(ft, II, nx, periodic_x)) +
                                    ⁻(sx*nnx)*(∇x⁺(ft, II, nx, periodic_x))) -
                            cfl_y * (⁺(sy*nny)*(-∇y⁻(ft, II, ny, periodic_y)) +
                                    ⁻(sy*nny)*(∇y⁺(ft, II, ny, periodic_y)))
        end
    end
end

function velocity_extension2!(V, u, inside, MIXED, n, h, NB, B, BT, pos)
    local cfl = 0.45
    local Vt = similar(V)
    local tmp1 = Vector{Float64}(undef, 0)
    local tmp2 = Vector{Float64}(undef, 0)
    for k = 1:1
        for j = 1:NB
            Vt .= V
            @inbounds @threads for II in inside
                s = mysign(u[II], h)
                nx = mysign(c∇x(u, II), c∇y(u, II))
                ny = mysign(c∇y(u, II), c∇x(u, II))
                V[II] = Vt[II] - cfl*(⁺(s*nx)*(-∇x⁻(Vt, II)) +
                ⁻(s*nx)*(∇x⁺(Vt, II)) +
                ⁺(s*ny)*(-∇y⁻(Vt, II)) +
                ⁻(s*ny)*(∇y⁺(Vt, II)))
            end
        end
        for II in MIXED
            st = static_stencil(V, II)
            itp = B * st * BT
            a = biquadratic(itp, pos[II].mid_point.x, pos[II].mid_point.y)
            b = abs(V[II]-a)
            c = sqrt(pos[II].mid_point.x^2 + pos[II].mid_point.y^2)
            push!(tmp1, V[II])
            push!(tmp2, V[II] + h*b/(1-c))
        end
    end
    return tmp1, tmp2
end

@inline faces_scalar(itp, II_0, II, x, y, dx, dy) = 
    @SVector [biquadratic(itp, (x[II] - x[II_0] - dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
        (y[II] - y[II_0])/(y[δy⁺(II_0)] - y[δy⁻(II_0)]))
    biquadratic(itp, (x[II] - x[II_0])/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
        (y[II] - y[II_0] - dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)]))
    biquadratic(itp, (x[II] - x[II_0] + dx/2)/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
        (y[II] - y[II_0])/(y[δy⁺(II_0)] - y[δy⁻(II_0)]))
    biquadratic(itp, (x[II] - x[II_0])/(x[δx⁺(II_0)] - x[δx⁻(II_0)]),
        (y[II] - y[II_0] + dy/2)/(y[δy⁺(II_0)] - y[δy⁻(II_0)]))]

function aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)
    st = static_stencil(u, II_0)
    B, BT = B_BT(II_0, x, y)
    itp = B * st * BT
    faces = faces_scalar(itp, II_0, II, x, y, dx[II], dy[II])
    @inbounds u_faces[II,:] .= faces

    return nothing
end

function interpolate_scalar!(grid, grid_u, grid_v, u, uu, uv)
    @unpack x, y, nx, ny, dx, dy, ind = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    u_faces = zeros(ny, nx, 4)

    @inbounds @threads for II in inside
        aux_interpolate_scalar!(II, II, u, x, y, dx, dy, u_faces)
    end

    @inbounds @threads for II in b_left[1][2:end-1]
        aux_interpolate_scalar!(δx⁺(II), II, u, x, y, dx, dy, u_faces)
    end

    II = b_left[1][1]
    II_0 = δy⁺(δx⁺(II))
    aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)

    II = b_left[1][end]
    II_0 = δy⁻(δx⁺(II))
    aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)

    @inbounds @threads for II in b_bottom[1][2:end-1]
        aux_interpolate_scalar!(δy⁺(II), II, u, x, y, dx, dy, u_faces)
    end

    @inbounds @threads for II in b_right[1][2:end-1]
        aux_interpolate_scalar!(δx⁻(II), II, u, x, y, dx, dy, u_faces)
    end

    II = b_right[1][1]
    II_0 = δy⁺(δx⁻(II))
    aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)

    II = b_right[1][end]
    II_0 = δy⁻(δx⁻(II))
    aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)

    @inbounds @threads for II in b_top[1][2:end-1]
        aux_interpolate_scalar!(δy⁻(II), II, u, x, y, dx, dy, u_faces)
    end

    @inbounds uu[:,1] .= @view u_faces[:,1,1]
    @inbounds uu[:,end] .= @view u_faces[:,end,3]

    @inbounds uv[1,:] .= @view u_faces[1,:,2]
    @inbounds uv[end,:] .= @view u_faces[end,:,4]

    @inbounds @threads for i = 1:grid_u.ny
        @inbounds uu[i,2:end-1] .= @views 0.5 * (u_faces[i,1:end-1,3] .+ u_faces[i,2:end,1])
    end
    @inbounds @threads for i = 1:grid_v.nx
        @inbounds uv[2:end-1,i] .= @views 0.5 * (u_faces[1:end-1,i,4] .+ u_faces[2:end,i,2])
    end
    
    return nothing
end

function breakup_f(u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, ϵ_break)
    u_copy = copy(u)
    local count = 0
    for II in NB_indices
        if !((II[2] == 1 || II[2] == nx) && !periodic_x) && II[1] != 1 && II[1] != ny
            if (u_copy[II]*u_copy[δx⁺(II, nx, periodic_x)] < 0) && (u_copy[II]*u_copy[δx⁻(II, nx, periodic_x)] < 0)
                if ((abs(u_copy[δx⁺(II, nx, periodic_x)]) < (dx[δx⁺(II, nx, periodic_x)]/2 - ϵ_break) &&
                     abs(u_copy[δx⁻(II, nx, periodic_x)]) < (dx[δx⁻(II, nx, periodic_x)]/2 - ϵ_break) &&
                     abs(u_copy[II]) > (dx[II]/2)) &&
                    ((abs(u_copy[δx⁺(δy⁺(II), nx, periodic_x)]) < (dx[δx⁺(δy⁺(II), nx, periodic_x)]/2 - ϵ_break) &&
                      abs(u_copy[δx⁻(δy⁺(II), nx, periodic_x)]) < (dx[δx⁻(δy⁺(II), nx, periodic_x)]/2 - ϵ_break) &&
                      abs(u_copy[δy⁺(II)]) > (dx[δy⁺(II)]/2)) ||
                     (abs(u_copy[δx⁺(δy⁻(II), nx, periodic_x)]) < (dx[δx⁺(δy⁻(II), nx, periodic_x)]/2 - ϵ_break) &&
                      abs(u_copy[δx⁻(δy⁻(II), nx, periodic_x)]) < (dx[δx⁻(δy⁻(II), nx, periodic_x)]/2 - ϵ_break) &&
                      abs(u_copy[δy⁻(II)]) > (dx[δy⁻(II)]/2))))
                    u[II] *= -1
                    count += 1
                end
            end
        end
        if !((II[1] == 1 || II[1] == ny) && !periodic_y) && II[2] != 1 && II[2] != nx
            if (u_copy[II]*u_copy[δy⁺(II, ny, periodic_y)] < 0) && (u_copy[II]*u_copy[δy⁻(II, ny, periodic_y)] < 0)
                if ((abs(u_copy[δy⁺(II, ny, periodic_y)]) < (dy[δy⁺(II, ny, periodic_y)]/2 - ϵ_break) &&
                     abs(u_copy[δy⁻(II, ny, periodic_y)]) < (dy[δy⁻(II, ny, periodic_y)]/2 - ϵ_break) &&
                     abs(u_copy[II]) > (dy[II]/2)) &&
                    ((abs(u_copy[δy⁺(δx⁺(II), ny, periodic_y)]) < (dy[δy⁺(δx⁺(II), ny, periodic_y)]/2 - ϵ_break) &&
                      abs(u_copy[δy⁻(δx⁺(II), ny, periodic_y)]) < (dy[δy⁻(δx⁺(II), ny, periodic_y)]/2 - ϵ_break) &&
                      abs(u_copy[δx⁺(II)]) > (dy[δx⁺(II)]/2)) ||
                     (abs(u_copy[δy⁺(δx⁻(II), ny, periodic_y)]) < (dy[δy⁺(δx⁻(II), ny, periodic_y)]/2 - ϵ_break) &&
                      abs(u_copy[δy⁻(δx⁻(II), ny, periodic_y)]) < (dy[δy⁻(δx⁻(II), ny, periodic_y)]/2 - ϵ_break) &&
                      abs(u_copy[δx⁻(II)]) > (dy[δx⁻(II)]/2))))
                    u[II] *= -1
                    count += 1
                end
            end
        end
    end
    return count
end