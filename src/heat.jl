@inline function diffusion_coefficients(a, i)
    avg = SA_F64[a[i,1]+a[i,3], a[i,2]+a[i,4]] .* 0.5
    V = SA_F64[a[i,5]+a[δx⁻(i),5], a[i,5]+a[δy⁻(i),5], a[i,5]+a[δx⁺(i),5], a[i,5]+a[δy⁺(i),5]] .* 0.5 .+ eps(.01)
    f = SA_F64[a[δx⁻(i),1]+a[δx⁻(i),3], a[δy⁻(i),2]+a[δy⁻(i),4], a[δx⁺(i),3]+a[δx⁺(i),1], a[δy⁺(i),4]+a[δy⁺(i),2]] .* 0.5
    bc = SA_F64[a[δx⁻(i),1]-a[i,3], a[δy⁻(i),2]-a[i,4], a[δx⁺(i),3]-a[i,1], a[δy⁺(i),4]-a[i,2]] .* 0.5
    return avg, V, f, bc
end

@inline function fill_matrices!(avg, V, f, bc, A, B, p, n, CFL)
    s = 0.
    bs = 0.
    a = (-n, -1, n, 1)
    @inbounds for (i,j) in zip(1:4,a)
        @inbounds A[p, p + j] = -(CFL*avg[2-mod(i,2)]*f[i])/V[i]
        @inbounds B[p, p + j] = (CFL*avg[2-mod(i,2)]*f[i])/V[i]
        s += (avg[2-mod(i,2)]^2)/V[i]
        bs -= (avg[2-mod(i,2)]*bc[i])/V[i]
    end
    return 2 + CFL*s, 2 - CFL*s, CFL*bs
end

function crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        Savg, SV, Sf, Sbc = diffusion_coefficients(SOL,II)
        Lavg, LV, Lf, Lbc = diffusion_coefficients(LIQ,II)
        AS[p,p], BS[p,p], SCUT_= fill_matrices!(Savg, SV, Sf, Sbc, AS, BS, p, n, CFL)
        AL[p,p], BL[p,p], LCUT_= fill_matrices!(Lavg, LV, Lf, Lbc, AL, BL, p, n, CFL)
        SCUT[p] = SCUT_*(θd - ϵ_κ*κ[II] - ϵ_V*V[II])
        LCUT[p] = LCUT_*(θd - ϵ_κ*κ[II] - ϵ_V*V[II])
    end
    return nothing
end

function crank_nicolson(SCUT, LCUT, SOL, LIQ, AS, AL, BS, BL, n, inside, CFL, θd, κ, ϵ_κ, ϵ_V, V, m, θ₀, aniso, sol_projection)
    @inbounds @threads for II in inside
        p = lexicographic(II, n)
        Savg, SV, Sf, Sbc = diffusion_coefficients(SOL,II)
        Lavg, LV, Lf, Lbc = diffusion_coefficients(LIQ,II)
        AS[p,p], BS[p,p], SCUT_= fill_matrices!(Savg, SV, Sf, Sbc, AS, BS, p, n, CFL)
        AL[p,p], BL[p,p], LCUT_= fill_matrices!(Lavg, LV, Lf, Lbc, AL, BL, p, n, CFL)
        ϵ_c = anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀)
        ϵ_v = anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀)
        SCUT[p] = SCUT_*(θd - ϵ_c*κ[II] - ϵ_V*V[II])
        LCUT[p] = LCUT_*(θd - ϵ_c*κ[II] - ϵ_v*V[II])
    end
    return nothing
end

@inline anisotropy(ϵ, m, θ, θ₀) = ϵ*(1 + 0.4*((8/3)*sin(0.5*m*(θ - θ₀))^4 - 1))


function Stefan_velocity!(TS, TL, sol_projection, liq_projection, V, MIXED, κ, ϵ_κ, ϵ_V, θd, h, m, θ₀, aniso)
    V .= 0
    @inbounds @threads for II in MIXED
        ϵ_c = ifelse(aniso, anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀), ϵ_κ)
        ϵ_v = ifelse(aniso, anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀), ϵ_κ)
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
