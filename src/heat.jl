@inline anisotropy(ϵ, m, θ, θ₀) = ϵ*(1 + 0.4*((8/3)*sin(0.5*m*(θ - θ₀))^4 - 1))

@inline function apply_curvature(num, grid, bc, D, all_indices)
    @unpack ϵ_κ, ϵ_V = num
    @unpack κ, V = grid

    @inbounds @threads for II in all_indices
        @inbounds bc[II] = D[II] - ϵ_κ*κ[II] - ϵ_V*V[II]
    end
    return nothing
end

@inline function apply_anisotropy(num, grid, bc, D, MIXED, sol_projection)
    @unpack ϵ_κ, ϵ_V, m, θ₀ = num
    @unpack κ, V = grid

    @inbounds @threads for II in MIXED
        ϵ_c = anisotropy(ϵ_κ, m, sol_projection[II].angle, θ₀)
        ϵ_v = anisotropy(ϵ_V, m, sol_projection[II].angle, θ₀)
        @inbounds bc[II] = D[II] - ϵ_c*κ[II] - ϵ_v*V[II]
    end
    return nothing
end

function Stefan_velocity!(num, grid, LS, V, TS, TL, MIXED, periodic_x, periodic_y)
    @unpack θd, ϵ_κ, ϵ_V, m, θ₀, aniso = num
    @unpack geoS, geoL, κ = LS

    V .= 0
    @inbounds @threads for II in MIXED
        ϵ_c = ifelse(aniso, anisotropy(ϵ_κ, m, geoS.projection[II].angle, θ₀), ϵ_κ)
        ϵ_v = ifelse(aniso, anisotropy(ϵ_V, m, geoS.projection[II].angle, θ₀), ϵ_V)
        θ_d = (θd - ϵ_c*κ[II] - ϵ_v*V[II])
        dTS = 0.
        dTL = 0.
        if geoS.projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, geoS.projection[II].angle, geoS.projection[II].point1, geoS.projection[II].point2, TS, II, periodic_x, periodic_y)
            dTS = normal_gradient(geoS.projection[II].d1, geoS.projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(grid, geoS.projection[II].angle, geoS.projection[II].point1, TS, II, periodic_x, periodic_y)
            dTS = normal_gradient(geoS.projection[II].d1, T_1, θ_d)
        end
        if geoL.projection[II].flag
            T_1, T_2 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, geoL.projection[II].point2, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, geoL.projection[II].d2, T_1, T_2, θ_d)
        else
            T_1 = interpolated_temperature(grid, geoL.projection[II].angle, geoL.projection[II].point1, TL, II, periodic_x, periodic_y)
            dTL = normal_gradient(geoL.projection[II].d1, T_1, θ_d)
        end
        V[II] = dTL + dTS
    end
    return nothing
end

@inline normal_gradient(d1, d2, T_1, T_2, θd) = (1/(d2 - d1)) * ((d2/d1) * (θd - T_1) - (d1/d2) * (θd - T_2))
@inline normal_gradient(d1, T_1, θd) = (θd - T_1)/d1

@inline function Acpm(grid, II_0, II)
    @unpack x, y, nx, ny = grid
    if II_0[1] < 3 || II_0[2] < 3
        Ac = B_BT(II, grid)[1]
        Ap = B_BT(II, grid)[1]
        Am = B_BT(II, grid)[1]
    elseif II_0[1] > ny-2 || II_0[2] > nx-2
        Ac = B_BT(II, grid)[1]
        Ap = B_BT(II, grid)[1]
        Am = B_BT(II, grid)[1]
    else
        Ac = B_BT(II_0, x, y)[1]
        Ap = B_BT(II_0, x, y, δx⁺)[1]
        Am = B_BT(II_0, x, y, δx⁻)[1]
    end

    return Ac, Ap, Am
end

@inline function interpolated_temperature(grid, α, P1, P2, temp, II, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy = grid
    T_1 = 0.
    T_2 = 0.
    if π/8 < α < 3π/8
        II_0 = δx⁺(δy⁺(II))
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        if α > π/4
            a = @view st[2,1:3]
            b = @view st[3,1:3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
        else
            a = @view st[1:3,2]
            b = @view st[1:3,3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
        end
    elseif 5π/8 < α < 7π/8
        II_0 = δx⁻(δy⁺(II))
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        if α < 3π/4
            a = @view st[2,1:3]
            b = @view st[3,1:3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
        else
            a = @view st[1:3,2]
            b = @view st[1:3,1]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
        end
    elseif -3π/8 < α < -π/8
        II_0 = δx⁺(δy⁻(II))
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        if α < -π/4
            a = @view st[2,1:3]
            b = @view st[1,1:3]
            T_1, T_2 = quadratic_interp(Ap, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
        else
            a = @view st[1:3,2]
            b = @view st[1:3,3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
        end
    elseif -7π/8 < α < -5π/8
        II_0 = δx⁻(δy⁻(II))
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        if α > -3π/4
            a = @view st[2,1:3]
            b = @view st[1,1:3]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
        else
            a = @view st[1:3,2]
            b = @view st[1:3,1]
            T_1, T_2 = quadratic_interp(Am, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
        end
    elseif -π/8 <= α <= π/8
        II_0 = δx⁺(II)
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        a = @view st[1:3,2]
        b = @view st[1:3,3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
    elseif 3π/8 <= α <= 5π/8
        II_0 = δy⁺(II)
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        a = @view st[2,1:3]
        b = @view st[3,1:3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
    elseif α >= 7π/8 || α <= -7π/8
        II_0 = δx⁻(II)
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        a = @view st[1:3,2]
        b = @view st[1:3,1]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.y/(2*dy[II]), P2.y/(2*dy[II]))
    elseif -5π/8 <= α <= -3π/8
        II_0 = δy⁻(II)
        st = static_stencil(temp, II_0, nx, ny, periodic_x, periodic_y)
        Ac, Ap, Am = Acpm(grid, II_0, II)
        a = @view st[2,1:3]
        b = @view st[1,1:3]
        T_1, T_2 = quadratic_interp(Ac, a, b, P1.x/(2*dx[II]), P2.x/(2*dx[II]))
    end
    return T_1, T_2
end

@inline function interpolated_temperature(grid, α, P1, temp, II, periodic_x, periodic_y)
    @unpack nx, ny, dx, dy = grid

    if II[1] == 1
        f = δy⁺
    elseif II[1] == grid.ny
        f = δy⁻
    else
        f = x->x
    end
    if II[2] == 1
        f = f ∘ δx⁺
    elseif II[2] == grid.nx
        f = f ∘ δx⁻
    else
        f = f ∘ (x->x)
    end

    A = B_BT(II, grid.x, grid.y, f)[1]
    T_1 = 0.
    st = static_stencil(temp, II, nx, ny, periodic_x, periodic_y)
    if π/4 <= α < 3π/4
        a = @view st[3,1:3]
        T_1 = quadratic_interp(A, a, P1.x/(2*dx[II]))
    elseif α >= 3π/4 || α < -3π/4
        a = @view st[1:3,1]
        T_1 = quadratic_interp(A, a, P1.y/(2*dy[II]))
    elseif -π/4 > α >= -3π/4
        a = @view st[1,1:3]
        T_1 = quadratic_interp(A, a, P1.x/(2*dx[II]))
    elseif -π/4 <= α < π/4
        a = @view st[1:3,3]
        T_1 = quadratic_interp(A, a, P1.y/(2*dy[II]))
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
