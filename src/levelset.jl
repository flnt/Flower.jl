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
        F = grad_IIOE(u, Vx, Vy, U, h, p, n)
        a_in, a_ou = inflow_outflow(F)
        S = sumloc(a_in, a_ou)
        A[p,p], B[p,p] = fill_matrices2!(a_in, a_ou, S, A, B, p, n, CFL)
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

@inline central_differences(u, II, h, n) =
    @SVector[minmod(Dxx(u, II, h), ifelse(in_bounds(δx⁺(II)[2], n), Dxx(u, δx⁺(II), h), 0.)),
            minmod(Dxx(u, II, h), ifelse(in_bounds(δx⁻(II)[2], n), Dxx(u, δx⁻(II), h), 0.)),
            minmod(Dyy(u, II, h), ifelse(in_bounds(δy⁺(II)[1], n), Dyy(u, δy⁺(II), h), 0.)),
            minmod(Dyy(u, II, h), ifelse(in_bounds(δy⁻(II)[1], n), Dyy(u, δy⁻(II), h), 0.))]

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

function FE_reinit(u, h, n, nb_reinit, BC_u, idx)
    @unpack inside, b_left, b_bottom, b_right, b_top = idx
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
            shift = central_differences(u, II, h, n)
            eno = finite_difference_eno(u, II, shift, h)

            if is_near_interface(u0, II)
                eno_interface = convert(Vector{Float64}, eno)
                h_ = 1e30
                for (JJ, i, j) in zip((δx⁺(II), δx⁻(II), δy⁺(II), δy⁻(II)), a, c)
                    if u0[II]*u0[JJ] < 0
                        uxx = minmod(f[i](u, II), ifelse(in_bounds(JJ[j], n), f[i](u, JJ), 0.))
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
        u .= tmp
        @sync begin
            @spawn bcs!(u, BC_u.left, h)
            @spawn bcs!(u, BC_u.right, h)
            @spawn bcs!(u, BC_u.bottom, h)
            @spawn bcs!(u, BC_u.top, h)
        end
    end
end

function hamiltonian(u, inside, h, n, BC_u)
    ham = similar(u)
    f, a, b, c = (Dxx, Dxx, Dyy, Dyy), (1, 2, 3, 4), (-1, 1, -1, 1), (2, 2, 1, 1)
    @inbounds @threads for II in inside
        h_ = h
        sign_u = sign(u[II])
        shift = central_differences(u, II, h, n)
        eno = finite_difference_eno(u, II, shift, h)
        if is_near_interface(u, II)
            eno_interface = convert(Vector{Float64}, eno)
            h_ = 1e30
            for (JJ, i, j) in zip((δx⁺(II), δx⁻(II), δy⁺(II), δy⁻(II)), a, c)
                if u[II]*u[JJ] < 0
                    uxx = minmod(f[i](u, II), ifelse(in_bounds(JJ[j], n), f[i](u, JJ), 0.))
                    D = (uxx/2 - u[II] - u[JJ])^2 - 4*u[II]*u[JJ]
                    Δx = root_extraction(u, uxx, D, II, JJ, h, 1e-10)
                    if Δx < h_ h_ = Δx end
                    eno_interface[i] = b[i]*(u[II]/Δx + (Δx/2) * shift[i])
                end
            end
            gdv = Godunov(sign_u, eno_interface)
        else
            gdv = Godunov(sign_u, eno)
        end
        ham[II] = gdv
    end
    @sync begin
        @spawn bcs!(u, BC_u.left, h)
        @spawn bcs!(u, BC_u.right, h)
        @spawn bcs!(u, BC_u.bottom, h)
        @spawn bcs!(u, BC_u.top, h)
    end
    return ham
end

function velocity_extension!(V, u, inside, n, h, NB, BC_u)
    local cfl = 0.45
    local Vt = similar(V)
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
            itp = B*st*BT
            a = biquadratic(itp, pos[II].mid_point.x, pos[II].mid_point.y)
            b = abs(V[II]-a)
            c = sqrt(pos[II].mid_point.x^2 + pos[II].mid_point.y^2)
            push!(tmp1, V[II])
            push!(tmp2, V[II] + h*b/(1-c))
        end
    end
    return tmp1, tmp2
end
