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

function IIOE(A, B, u, V, inside, CFL, n)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        U = diamond(u, p, n)
        F = grad(u, V, U, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        if S[1] < -S[2]
            D = quadratic_recons(u, U, p, n)
            D_ = quadratic_recons(D)
            F = 2*grad(u, V, U, D, D_[1], p, n)
            a_in, a_ou = inflow_outflow(F)
            S = sumloc(a_in, a_ou)
        end
        A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, S, A, B, p, n, CFL)
    end
end

function IIOE(grid, A, B, u, V, CFL, periodic_x, periodic_y)
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


function level_update_IIOE!(grid, grid_u, grid_v, A, B, θ_out, MIXED, τ, periodic_x, periodic_y)
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

function S2IIOE!(grid, grid_u, grid_v, A, B, utmp, u, θ_out, MIXED, τ, periodic_x, periodic_y)
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

function level_update_koenig!(A, B, u, V, Vx, Vy, inside, CFL, h, n)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        U = diamond(u, p, n)
        F = grad(u, V, U, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        if S[1] < -S[2]
            D = quadratic_recons(u, U, p, n)
            D_ = quadratic_recons(D)
            F = 2*grad(u, V, U, D, D_[1], p, n)
            a_in, a_ou = inflow_outflow(F)
            S = sumloc(a_in, a_ou)
        end
        F2 = advection(u, Vx, Vy, U, h, p, n)
        a_in2, a_ou2 = inflow_outflow(F2)
        a_in_tot = a_in + a_in2
        a_ou_tot = a_ou + a_ou2
        S = sumloc(a_in_tot, a_ou_tot)
        A[p,p], B[p,p] = fill_matrices2!(a_in_tot, a_ou_tot, S, A, B, p, n, CFL)
    end
end

@inline function fill_matrices2!(a_in::SArray{Tuple{4},Float64,1,4}, a_ou::SArray{Tuple{4},Float64,1,4},
    S, A, B, p, n, CFL)
    a = (-n, -1, n, 1)
    @inbounds for (i,j) in zip(1:4,a)
        @inbounds A[p, p + j] = -CFL*a_in[i]
        @inbounds B[p, p + j] = CFL*a_ou[i]
    end
    return 2 + CFL*S[1], 2 - CFL*S[2]
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
            minmod(Dyy(u, II, dy, ny, per_y), ifelse(in_bounds(δy⁻(II, ny, per_y)[1], ny, per_y), Dyy(u, δy⁻(II, ny, per_y), dy, ny, per_y), 0.))]

@inline finite_difference_eno(u, II, a::AbstractArray, dx, dy, nx, ny, per_x, per_y) =
    @SVector[2*∇x⁺(u, II, nx, per_x)/(dx[δx⁺(II, nx, per_x)]+dx[II]) - ((dx[δx⁺(II, nx, per_x)]+dx[II])/4)*a[1],
            -2*∇x⁻(u, II, nx, per_x)/(dx[δx⁻(II, nx, per_x)]+dx[II]) + ((dx[δx⁻(II, nx, per_x)]+dx[II])/4)*a[2],
            2*∇y⁺(u, II, ny, per_y)/(dy[δy⁺(II, ny, per_y)]+dy[II]) - ((dy[δy⁺(II, ny, per_y)]+dy[II])/4)*a[3],
            -2*∇y⁻(u, II, ny, per_y)/(dy[δy⁻(II, ny, per_y)]+dy[II]) + ((dy[δy⁻(II, ny, per_y)]+dy[II])/4)*a[4]]

@inline Godunov(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + max(⁻(a[3])^2, ⁺(a[4])^2)),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + max(⁺(a[3])^2, ⁻(a[4])^2)))

@inline root_extraction(u, uxx, D, II, JJ, h, eps) =
    ifelse(abs(uxx) > eps,
    h*(0.5 + ((u[II] - u[JJ] - mysign(u[II] - u[JJ]) * sqrt(D))/uxx)),
    h*(u[II]/(u[II]-u[JJ])))

function levelset_BC!(u, BC_u, dx, dy, periodic_x, periodic_y)
    if periodic_x
        @sync begin
            @spawn bcs!(u, BC_u.bottom, dy[1,1])
            @spawn bcs!(u, BC_u.top, dy[end,1])
        end
    elseif periodic_y
        @sync begin
            @spawn bcs!(u, BC_u.left, dx[1,1])
            @spawn bcs!(u, BC_u.right, dx[1,end])
        end
    else
        @sync begin
            @spawn bcs!(u, BC_u.left, dx[1,1])
            @spawn bcs!(u, BC_u.right, dx[1,end])
            @spawn bcs!(u, BC_u.bottom, dy[1,1])
            @spawn bcs!(u, BC_u.top, dy[end,1])
        end
    end
end

function FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy = grid
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

    local cfl = 0.45
    levelset_BC!(u, BC_u, dx, dy, periodic_x, periodic_y)
    u0 = copy(u)
    tmp = similar(u)
    a, b, c = (1, 2, 3, 4), (-1, 1, -1, 1), (2, 2, 1, 1)
    f, d, e = (Dxx, Dxx, Dyy, Dyy), (nx, nx, ny, ny), (periodic_x, periodic_x, periodic_y, periodic_y)

    for i = 1:nb_reinit
        @inbounds @threads for II in indices
            h_ = min(0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                     0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
            sign_u0 = sign(u0[II])
            shift = central_differences(u, II, dx, dy, nx, ny, periodic_x, periodic_y)
            eno = finite_difference_eno(u, II, shift, dx, dy, nx, ny, periodic_x, periodic_y)

            if is_near_interface(u0, II, nx, ny, periodic_x, periodic_y)
                eno_interface = convert(Vector{Float64}, eno)
                h_ = 1e30
                g = (0.5*(dx[II]+dx[δx⁺(II, nx, periodic_x)]), 0.5*(dx[II]+dx[δx⁻(II, nx, periodic_x)]),
                     0.5*(dy[II]+dy[δy⁺(II, ny, periodic_y)]), 0.5*(dy[II]+dy[δy⁻(II, ny, periodic_y)]))
                for (JJ, i, j, k) in zip((δx⁺(II, nx, periodic_x), δx⁻(II, nx, periodic_x), δy⁺(II, ny, periodic_y), δy⁻(II, ny, periodic_y)), a, c, g)
                    if u0[II]*u0[JJ] < 0
                        uxx = minmod(f[i](u, II, d[i], e[i]), ifelse(in_bounds(JJ[j], (j == 1 ? ny : nx), (j == 1 ? periodic_y : periodic_x)), f[i](u, JJ, d[i], e[i]), 0.))
                        D = (uxx/2 - u0[II] - u0[JJ])^2 - 4*u0[II]*u0[JJ]
                        Δx = root_extraction(u0, uxx, D, II, JJ, k, 1e-10)
                        if Δx < h_ h_ = Δx end
                        eno_interface[i] = b[i]*(u[II]/Δx + (Δx/2) * shift[i])
                    end
                end
                gdv = Godunov(sign_u0, eno_interface)
            else
                gdv = Godunov(sign_u0, eno)
            end

            tmp[II] = u[II] - cfl * h_ * sign_u0 * (gdv - 1.0)
        end
        u[indices] .= tmp[indices]

        levelset_BC!(u, BC_u, dx, dy, periodic_x, periodic_y)
    end
end

function velocity_extension!(grid, inside, NB)
    @unpack dx, dy, V, u = grid

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

function velocity_extension!(grid, indices, NB, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, V, u, ind = grid

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

function field_extension!(grid, f, indices, NB, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy, u = grid

    local cfl = 0.45
    local ft = similar(f)

    τ = cfl * max(dx..., dy...)

    for j = 1:NB
        ft .= f
        @inbounds @threads for II in indices
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