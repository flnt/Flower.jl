@inline anisotropy(ϵ, m, θ, θ₀) = ϵ*(1 + 0.4*((8/3)*sin(0.5*m*(θ - θ₀))^4 - 1))

@inline function apply_curvature(bc, D, κ, ϵ_κ, ϵ_V, V, all_indices)
    @inbounds @threads for II in all_indices
        @inbounds bc[II] = D[II] - ϵ_κ*κ[II] - ϵ_V*V[II]
    end
    return nothing
end

@inline function apply_anisotropy(bc, D, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
    @inbounds @threads for II in MIXED
        ϵ_c = anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀)
        ϵ_v = anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀)
        @inbounds bc[II] = D[II] - ϵ_c*κ[II] - ϵ_v*V[II]
    end
    return nothing
end

function crank_nicolson!(L, A, B, C, cap, τ, n, Δ, all_indices)
    @inbounds V = cap[:,:,5] .* Δ^2

    @inbounds A .= -L .* τ
    @inbounds B .= (L .- C)  .* τ
    @inbounds @threads for II in all_indices
        pII = lexicographic(II, n)
        @inbounds B[pII,pII] += V[II] * 2.
        @inbounds A[pII,pII] += V[II] * 2.
    end
    return nothing
end

function set_heat!(HT, DT, θd, bcT, BC_T, κ, ϵ_κ, ϵ_V, V, m, θ₀, τ, aniso, convection,
            LT, CUTT, CAP, A, B, C,
            Du, Dv, Hu, Hv, u, v, BC_u, BC_v, CUTCT,
            Gx, Gy, CUTGx, CUTGy, ftcGx, ftcGy,
            centroid, mid_point, sol_projection,
            all_indices, inside, MIXED, empty, b_left, b_bottom, b_right, b_top, n, Δ)
    # rhs = B*vec(T) .+ τ*CUTT
    # for II in empty
    #     pII = lexicographic(II, n)
    #     rhs[pII] = 0.
    # end
    HT .= 0.
    @inbounds @threads for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HT[II] = distance(mid_point[II], centroid[II]) * Δ
    end

    DT .= θd
    apply_curvature(bcT, DT, κ, ϵ_κ, ϵ_V, V, all_indices)
    if aniso
        apply_anisotropy(bcT, DT, MIXED, κ, ϵ_κ, ϵ_V, V, m, θ₀, sol_projection)
    end
    bcTx, bcTy = set_bc_bnds(dir, bcT, HT, BC_T)

    laplacian!(dir, LT, CUTT, bcTx, bcTy, CAP, n, Δ, BC_T, inside, empty,
                b_left[1], b_bottom[1], b_right[1], b_top[1])
    if convection
        bcU, bcV = set_bc_bnds(dir, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
        scalar_convection!(dir, C, CUTCT, u, v, bcTx, bcTy, bcU, bcV, CAP, n, Δ, BC_T, inside, b_left[1], b_bottom[1], b_right[1], b_top[1])
    end
    crank_nicolson!(LT, A, B, C, CAP, τ, n, Δ, all_indices)

    gradient!(dir, Gx, Gy, CUTGx, CUTGy, bcTx, bcTy, CAP, n, Δ, BC_T, all_indices, b_left[1], b_bottom[1], b_right[1], b_top[1])
    face_to_cell_gradient!(dir, ftcGx, ftcGy, Gx, Gy, CAP, n, all_indices)

    return nothing
end


function Stefan_velocity!(TS, TL, sol_projection, liq_projection, V, MIXED, κ, ϵ_κ, ϵ_V, θd, h, m, θ₀, aniso)
    V .= 0
    @inbounds @threads for II in MIXED
        ϵ_c = ifelse(aniso, anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀), ϵ_κ)
        ϵ_v = ifelse(aniso, anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀), ϵ_V)
        θ_d = (θd - ϵ_c*κ[II] - ϵ_v*V[II])
        dTS = 0.
        dTL = 0.
        if sol_projection[II].flag
            T_1, T_2 = interpolated_temperature(sol_projection[II].angle, sol_projection[II].point1, sol_projection[II].point2, TS, II)
            dTS = normal_gradient(sol_projection[II].d1, sol_projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(sol_projection[II].angle, sol_projection[II].point1, TS, II)
            dTS = normal_gradient(sol_projection[II].d1, T_1, θ_d)
        end
        if liq_projection[II].flag
            T_1, T_2 = interpolated_temperature(liq_projection[II].angle, liq_projection[II].point1, liq_projection[II].point2, TL, II)
            dTL = normal_gradient(liq_projection[II].d1, liq_projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(liq_projection[II].angle, liq_projection[II].point1, TL, II)
            dTL = normal_gradient(liq_projection[II].d1, T_1, θ_d)
        end
        V[II] = (dTL + dTS)/h
    end
    return nothing
end

function Stefan_vel!(TS, TL, V, MIXED, GxTS, GxTL, GyTS, GyTL, ftcGxTS, ftcGxTL, ftcGyTS, ftcGyTL, SCUTGxT, LCUTGxT, SCUTGyT, LCUTGyT, iMxS, iMyS, iMxL, iMyL, n, Δ)
    gS = reshape(sqrt.((ftcGxTS * iMxS * (GxTS * vec(TS) .+ SCUTGxT)).^2 .+
                       (ftcGyTS * iMyS * (GyTS * vec(TS) .+ SCUTGyT)).^2), (n,n))
    gL = reshape(sqrt.((ftcGxTL * iMxL * (GxTL * vec(TL) .+ LCUTGxT)).^2 .+
                       (ftcGyTL * iMyL * (GyTL * vec(TL) .+ LCUTGyT)).^2), (n,n))

    V .= 0.
    @inbounds @threads for II in MIXED
        @inbounds V[II] = -(gS[II] + gL[II])
    end

    return nothing
end

@inline normal_gradient(d1, d2, T_1, T_2, θd) = (1/(d2 - d1)) * ((d2/d1) * (θd - T_1) - (d1/d2) * (θd - T_2))
@inline normal_gradient(d1, T_1, θd) = (θd - T_1)/d1

@inline function interpolated_temperature(α, P1, P2, temp, II)
    T_1 = 0.
    T_2 = 0.
    Ac = @SMatrix [0.5 -1.0 0.5; -0.5 -0.0 0.5; 0.0 1.0 0.0]
    Ap = @SMatrix [0.5 -1.0 0.5; -1.5 2.0 -0.5; 1.0 0.0 0.0]
    Am = @SMatrix [0.5 -1.0 0.5; 0.5 -2.0 1.5; 0.0 0.0 1.0]
    if π/8 < α < 3π/8
        st = static_stencil(temp, δx⁺(δy⁺(II)))
        if α > π/4
            a = @view st[2,1:3]
            b = @view st[3,1:3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.x, P2.x)
        else
            a = @view st[1:3,2]
            b = @view st[1:3,3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.y, P2.y)
        end
    elseif 5π/8 < α < 7π/8
        st = static_stencil(temp, δx⁻(δy⁺(II)))
        if α < 3π/4
            a = @view st[2,1:3]
            b = @view st[3,1:3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.x, P2.x)
        else
            a = @view st[1:3,2]
            b = @view st[1:3,1]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.y, P2.y)
        end
    elseif -3π/8 < α < -π/8
        st = static_stencil(temp, δx⁺(δy⁻(II)))
        if α < -π/4
            a = @view st[2,1:3]
            b = @view st[1,1:3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.x, P2.x)
        else
            a = @view st[1:3,2]
            b = @view st[1:3,3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.y, P2.y)
        end
    elseif -7π/8 < α < -5π/8
        st = static_stencil(temp, δx⁻(δy⁻(II)))
        if α > -3π/4
            a = @view st[2,1:3]
            b = @view st[1,1:3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.x, P2.x)
        else
            a = @view st[1:3,2]
            b = @view st[1:3,1]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.y, P2.y)
        end
    elseif -π/8 <= α <= π/8
        st = static_stencil(temp, δx⁺(II))
        a = @view st[1:3,2]
        b = @view st[1:3,3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.y, P2.y)
    elseif 3π/8 <= α <= 5π/8
        st = static_stencil(temp, δy⁺(II))
        a = @view st[2,1:3]
        b = @view st[3,1:3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.x, P2.x)
    elseif α >= 7π/8 || α <= -7π/8
        st = static_stencil(temp, δx⁻(II))
        a = @view st[1:3,2]
        b = @view st[1:3,1]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.y, P2.y)
    elseif -5π/8 <= α <= -3π/8
        st = static_stencil(temp, δy⁻(II))
        a = @view st[2,1:3]
        b = @view st[1,1:3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.x, P2.x)
    end
    return T_1, T_2
end

@inline function interpolated_temperature(α, P1, temp, II)
    A = @SMatrix [0.5 -1.0 0.5; -0.5 -0.0 0.5; 0.0 1.0 0.0]
    T_1 = 0.
    st = static_stencil(temp, II)
    if π/4 <= α < 3π/4
        a = @view st[3,1:3]
        T_1 = quadratic_interp(A, a, P1.x)
    elseif α >= 3π/4 || α < -3π/4
        a = @view st[1:3,1]
        T_1 = quadratic_interp(A, a, P1.y)
    elseif -π/4 > α >= -3π/4
        a = @view st[1,1:3]
        T_1 = quadratic_interp(A, a, P1.x)
    elseif -π/4 <= α < π/4
        a = @view st[1:3,3]
        T_1 = quadratic_interp(A, a, P1.y)
    end
    return T_1
end

@inline function quadratic_interp(A, a, b, xP1, xP2)
    c_a = A*a
    c_b = A*b
    return c_a[1]*xP1^2 + c_a[2]*xP1 + c_a[3], c_b[1]*xP2^2 + c_b[2]*xP2 + c_b[3]
end

@inline function quadratic_interp(A, a, xP1)
    c_a = A*a
    return c_a[1]*xP1^2 + c_a[2]*xP1 + c_a[3]
end
