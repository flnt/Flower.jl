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
    h::Float64, p::Int64, n::Int64)
    F = SA_F64[g[p]*(D_-D[1]),
    g[p]*(D_-D[2]),
    g[p]*(D_-D[3]),
    g[p]*(D_-D[4])]
    return F
end

function grad(a, g,
    U::SArray{Tuple{4},Float64,1,4},
    h::Float64, p::Int64, n::Int64)
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

function IIOE(A, B, u, V, inside, CFL, h, n)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        U = diamond(u, p, n)
        F = grad(u, V, U, h, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        if S[1] < -S[2]
            D = quadratic_recons(u, U, p, n)
            D_ = quadratic_recons(D)
            F = 2*grad(u, V, U, D, D_[1], h, p, n)
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
        F = grad(u, V, U, h, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        if S[1] < -S[2]
            D = quadratic_recons(u, U, p, n)
            D_ = quadratic_recons(D)
            F = 2*grad(u, V, U, D, D_[1], h, p, n)
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

@inline central_differences(u, II, h, nx, ny) =
    @SVector[minmod(Dxx(u, II, h), ifelse(in_bounds(δx⁺(II)[2], nx), Dxx(u, δx⁺(II), h), 0.)),
            minmod(Dxx(u, II, h), ifelse(in_bounds(δx⁻(II)[2], nx), Dxx(u, δx⁻(II), h), 0.)),
            minmod(Dyy(u, II, h), ifelse(in_bounds(δy⁺(II)[1], ny), Dyy(u, δy⁺(II), h), 0.)),
            minmod(Dyy(u, II, h), ifelse(in_bounds(δy⁻(II)[1], ny), Dyy(u, δy⁻(II), h), 0.))]

@inline finite_difference_eno(u, II, a::AbstractArray, h) =
    @SVector[∇x⁺(u, II)/h - (h/2)*a[1],
            -∇x⁻(u, II)/h + (h/2)*a[2],
            ∇y⁺(u, II)/h - (h/2)*a[3],
            -∇y⁻(u, II)/h + (h/2)*a[4]]

@inline Godunov(s, a::AbstractArray) = ifelse(s >= 0,
    sqrt(max(⁻(a[1])^2, ⁺(a[2])^2) + max(⁻(a[3])^2, ⁺(a[4])^2)),
    sqrt(max(⁺(a[1])^2, ⁻(a[2])^2) + max(⁺(a[3])^2, ⁻(a[4])^2)))

@inline root_extraction(u, uxx, D, II, JJ, h, eps) =
    ifelse(abs(uxx) > eps,
    h*(0.5 + ((u[II] - u[JJ] - mysign(u[II] - u[JJ]) * sqrt(D))/uxx)),
    h*(u[II]/(u[II]-u[JJ])))

function FE_reinit(grid, ind, u, h, nb_reinit, BC_u)
    @unpack nx, ny = grid
    @unpack inside = ind

    local cfl = 0.45
    @sync begin
        @spawn bcs!(u, BC_u.left, h)
        @spawn bcs!(u, BC_u.right, h)
        @spawn bcs!(u, BC_u.bottom, h)
        @spawn bcs!(u, BC_u.top, h)
    end
    u0 = copy(u)
    tmp = similar(u)
    f, a, b, c = (Dxx, Dxx, Dyy, Dyy), (1, 2, 3, 4), (-1, 1, -1, 1), (2, 2, 1, 1)
    for i = 1:nb_reinit
        @inbounds @threads for II in inside
            h_ = h
            sign_u0 = sign(u0[II])
            shift = central_differences(u, II, h, nx, ny)
            eno = finite_difference_eno(u, II, shift, h)

            if is_near_interface(u0, II)
                eno_interface = convert(Vector{Float64}, eno)
                h_ = 1e30
                for (JJ, i, j) in zip((δx⁺(II), δx⁻(II), δy⁺(II), δy⁻(II)), a, c)
                    if u0[II]*u0[JJ] < 0
                        uxx = minmod(f[i](u, II), ifelse(in_bounds(JJ[j], (j == 1 ? ny : nx)), f[i](u, JJ), 0.))
                        D = (uxx/2 - u0[II] - u0[JJ])^2 - 4*u0[II]*u0[JJ]
                        Δx = root_extraction(u0, uxx, D, II, JJ, h, 1e-10)
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
            @spawn bcs!(u, BC_u.left, h)
            @spawn bcs!(u, BC_u.right, h)
            @spawn bcs!(u, BC_u.bottom, h)
            @spawn bcs!(u, BC_u.top, h)
        end
    end
end

function velocity_extension!(grid, inside, NB, Δ)
    @unpack V, u = grid

    local cfl = 0.45
    local Vt = similar(V)
    for j = 1:NB
        Vt .= V
        @inbounds @threads for II in inside
            s = mysign(u[II], Δ)
            nx = mysign(c∇x(u, II), c∇y(u, II))
            ny = mysign(c∇y(u, II), c∇x(u, II))
            V[II] = Vt[II] - cfl*(⁺(s*nx)*(-∇x⁻(Vt, II)) +
            ⁻(s*nx)*(∇x⁺(Vt, II)) +
            ⁺(s*ny)*(-∇y⁻(Vt, II)) +
            ⁻(s*ny)*(∇y⁺(Vt, II)))
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

@inline faces_scalar(itp) = [biquadratic(itp, -0.5, 0.0),
                                biquadratic(itp, 0.0, -0.5),
                                biquadratic(itp, 0.5, 0.0),
                                biquadratic(itp, 0.0, 0.5)]

@inline faces_scalar_left(itp) = [biquadratic(itp, -1.5, 0.0),
                                biquadratic(itp, -1.0, -0.5),
                                biquadratic(itp, -0.5, 0.0),
                                biquadratic(itp, -1.0, 0.5)]

@inline faces_scalar_left_bottom(itp) = [biquadratic(itp, -1.5, -1.0),
                                biquadratic(itp, -1.0, -1.5),
                                biquadratic(itp, -0.5, -1.0),
                                biquadratic(itp, -1.0, -0.5)]

@inline faces_scalar_left_top(itp) = [biquadratic(itp, -1.5, 1.0),
                                biquadratic(itp, -1.0, 0.5),
                                biquadratic(itp, -0.5, .0),
                                biquadratic(itp, -1.0, 1.5)]

@inline faces_scalar_bottom(itp) = [biquadratic(itp, -0.5, -1.0),
                                biquadratic(itp, 0.0, -1.5),
                                biquadratic(itp, 0.5, -1.0),
                                biquadratic(itp, 0.0, -0.5)]

@inline faces_scalar_right(itp) = [biquadratic(itp, 0.5, 0.0),
                                biquadratic(itp, 1.0, -0.5),
                                biquadratic(itp, 1.5, 0.0),
                                biquadratic(itp, 1.0, 0.5)]

@inline faces_scalar_right_bottom(itp) = [biquadratic(itp, 0.5, -1.0),
                                biquadratic(itp, 1.0, -1.5),
                                biquadratic(itp, 1.5, -1.0),
                                biquadratic(itp, 1.0, -0.5)]

@inline faces_scalar_right_top(itp) = [biquadratic(itp, 0.5, 1.0),
                                biquadratic(itp, 1.0, 0.5),
                                biquadratic(itp, 1.5, 1.0),
                                biquadratic(itp, 1.0, 1.5)]

@inline faces_scalar_top(itp) = [biquadratic(itp, -0.5, 1.0),
                                biquadratic(itp, 0.0, 0.5),
                                biquadratic(itp, 0.5, 1.0),
                                biquadratic(itp, 0.0, 1.5)]

function interpolate_scalar!(num, grid, grid_u, grid_v)
    @unpack B, BT = num
    @unpack nx, ny, ind, u = grid
    @unpack inside, b_left, b_bottom, b_right, b_top = ind

    u_faces = zeros(ny, nx, 4)

    @inbounds @threads for II in inside
        st = static_stencil(u, II)
        itp = B * st * BT
        faces = faces_scalar(itp)
        @inbounds u_faces[II,:] .= faces
    end

    @inbounds for II in b_left[1][2:end-1]
        st = static_stencil(u, δx⁺(II))
        itp = B * st * BT
        faces = faces_scalar_left(itp)
        @inbounds u_faces[II,:] .= faces
    end

    II = b_left[1][1]
    st = static_stencil(u, δy⁺(δx⁺(II)))
    itp = B * st * BT
    faces = faces_scalar_left_bottom(itp)
    @inbounds u_faces[II,:] .= faces

    II = b_left[1][end]
    st = static_stencil(u, δy⁻(δx⁺(II)))
    itp = B * st * BT
    faces = faces_scalar_left_top(itp)
    @inbounds u_faces[II,:] .= faces

    @inbounds for II in b_bottom[1][2:end-1]
        st = static_stencil(u, δy⁺(II))
        itp = B * st * BT
        faces = faces_scalar_bottom(itp)
        @inbounds u_faces[II,:] .= faces
    end

    @inbounds for II in b_right[1][2:end-1]
        st = static_stencil(u, δx⁻(II))
        itp = B * st * BT
        faces = faces_scalar_right(itp)
        @inbounds u_faces[II,:] .= faces
    end

    II = b_right[1][1]
    st = static_stencil(u, δy⁺(δx⁻(II)))
    itp = B * st * BT
    faces = faces_scalar_right_bottom(itp)
    @inbounds u_faces[II,:] .= faces

    II = b_right[1][end]
    st = static_stencil(u, δy⁻(δx⁻(II)))
    itp = B * st * BT
    faces = faces_scalar_right_top(itp)
    @inbounds u_faces[II,:] .= faces

    @inbounds for II in b_top[1][2:end-1]
        st = static_stencil(u, δy⁻(II))
        itp = B * st * BT
        faces = faces_scalar_top(itp)
        @inbounds u_faces[II,:] .= faces
    end

    @inbounds grid_u.u[:,1] .= u_faces[:,1,1]
    @inbounds grid_u.u[:,end] .= u_faces[:,end,3]

    @inbounds grid_v.u[1,:] .= u_faces[1,:,2]
    @inbounds grid_v.u[end,:] .= u_faces[end,:,4]

    @inbounds @threads for i = 1:grid_u.ny
        @inbounds grid_u.u[i,2:end-1] .= 0.5 * (u_faces[i,1:end-1,3] .+ u_faces[i,2:end,1])
    end
    @inbounds @threads for i = 1:grid_v.nx
        @inbounds grid_v.u[2:end-1,i] .= 0.5 * (u_faces[1:end-1,i,4] .+ u_faces[2:end,i,2])
    end
    
    return nothing
end