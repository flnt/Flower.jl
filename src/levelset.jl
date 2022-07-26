function diamond(a, p, n)
    U = @SVector[0.25(a[p]+a[p-n]+a[p-n-1]+a[p-1]),
    0.25(a[p]+a[p-1]+a[p+n-1]+a[p+n]),
    0.25(a[p]+a[p+n]+a[p+n+1]+a[p+1]),
    0.25(a[p]+a[p+1]+a[p-n+1]+a[p-n])]
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


function level_update_IIOE!(A, B, u, Vx, Vy, inside, CFL, h, n)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        U = diamond(u, p, n)
        F = advection(u, Vx, Vy, U, h, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, S, A, B, p, n, CFL)
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

@inline central_differences(u, II, dx, dy, nx, ny) =
    @SVector[minmod(Dxx(u, II, dx), ifelse(in_bounds(δx⁺(II)[2], nx), Dxx(u, δx⁺(II), dx), 0.)),
            minmod(Dxx(u, II, dx), ifelse(in_bounds(δx⁻(II)[2], nx), Dxx(u, δx⁻(II), dx), 0.)),
            minmod(Dyy(u, II, dy), ifelse(in_bounds(δy⁺(II)[1], ny), Dyy(u, δy⁺(II), dy), 0.)),
            minmod(Dyy(u, II, dy), ifelse(in_bounds(δy⁻(II)[1], ny), Dyy(u, δy⁻(II), dy), 0.))]

@inline finite_difference_eno(u, II, a::AbstractArray, dx, dy) =
    @SVector[2*∇x⁺(u, II)/(dx[δx⁺(II)]+dx[II]) - ((dx[δx⁺(II)]+dx[II])/4)*a[1],
            -2*∇x⁻(u, II)/(dx[δx⁻(II)]+dx[II]) + ((dx[δx⁻(II)]+dx[II])/4)*a[2],
            2*∇y⁺(u, II)/(dy[δy⁺(II)]+dy[II]) - ((dy[δy⁺(II)]+dy[II])/4)*a[3],
            -2*∇y⁻(u, II)/(dy[δy⁻(II)]+dy[II]) + ((dy[δy⁻(II)]+dy[II])/4)*a[4]]

@inline Godunov(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + max(⁻(a[3])^2, ⁺(a[4])^2)),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + max(⁺(a[3])^2, ⁻(a[4])^2)))

@inline root_extraction(u, uxx, D, II, JJ, h, eps) =
    ifelse(abs(uxx) > eps,
    h*(0.5 + ((u[II] - u[JJ] - mysign(u[II] - u[JJ]) * sqrt(D))/uxx)),
    h*(u[II]/(u[II]-u[JJ])))

function FE_reinit(grid, ind, u, nb_reinit, BC_u)
    @unpack nx, ny, dx, dy = grid
    @unpack inside = ind

    local cfl = 0.45
    @sync begin
        @spawn bcs!(u, BC_u.left, dx[1,1])
        @spawn bcs!(u, BC_u.right, dx[1,end])
        @spawn bcs!(u, BC_u.bottom, dy[1,1])
        @spawn bcs!(u, BC_u.top, dy[end,1])
    end
    u0 = copy(u)
    tmp = similar(u)
    f, a, b, c = (Dxx, Dxx, Dyy, Dyy), (1, 2, 3, 4), (-1, 1, -1, 1), (2, 2, 1, 1)
    for i = 1:nb_reinit
        @inbounds @threads for II in inside
            h_ = min(0.5*(dx[II]+dx[δx⁺(II)]), 0.5*(dx[II]+dx[δx⁻(II)]),
                     0.5*(dy[II]+dy[δy⁺(II)]), 0.5*(dy[II]+dy[δy⁻(II)]))
            sign_u0 = sign(u0[II])
            shift = central_differences(u, II, dx, dy, nx, ny)
            eno = finite_difference_eno(u, II, shift, dx, dy)

            if is_near_interface(u0, II)
                eno_interface = convert(Vector{Float64}, eno)
                h_ = 1e30
                d = (0.5*(dx[II]+dx[δx⁺(II)]), 0.5*(dx[II]+dx[δx⁻(II)]),
                     0.5*(dy[II]+dy[δy⁺(II)]), 0.5*(dy[II]+dy[δy⁻(II)]))
                for (JJ, i, j, k) in zip((δx⁺(II), δx⁻(II), δy⁺(II), δy⁻(II)), a, c, d)
                    if u0[II]*u0[JJ] < 0
                        uxx = minmod(f[i](u, II), ifelse(in_bounds(JJ[j], (j == 1 ? ny : nx)), f[i](u, JJ), 0.))
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
        u[inside] .= tmp[inside]
        @sync begin
            @spawn bcs!(u, BC_u.left, dx[1,1])
            @spawn bcs!(u, BC_u.right, dx[1,end])
            @spawn bcs!(u, BC_u.bottom, dy[1,1])
            @spawn bcs!(u, BC_u.top, dy[end,1])
        end
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

    @inbounds uu[:,1] .= u_faces[:,1,1]
    @inbounds uu[:,end] .= u_faces[:,end,3]

    @inbounds uv[1,:] .= u_faces[1,:,2]
    @inbounds uv[end,:] .= u_faces[end,:,4]

    @inbounds @threads for i = 1:grid_u.ny
        @inbounds uu[i,2:end-1] .= 0.5 * (u_faces[i,1:end-1,3] .+ u_faces[i,2:end,1])
    end
    @inbounds @threads for i = 1:grid_v.nx
        @inbounds uv[2:end-1,i] .= 0.5 * (u_faces[1:end-1,i,4] .+ u_faces[2:end,i,2])
    end
    
    return nothing
end