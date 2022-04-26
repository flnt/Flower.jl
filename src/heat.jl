function set_bc_bnds(D, H, BC)
    Dx = similar(D)
    Dy = similar(D)
    Dx .= D
    Dy .= D

    if BC.left.f == neumann
        Dx[:,1] .= H[:,1] .* BC.left.val
    else
        Dx[:,1] .= BC.left.val
    end
    if BC.bottom.f == neumann
        Dy[1,:] .= H[1,:] .* BC.bottom.val 
    else
        Dy[1,:] .= BC.bottom.val
    end
    if BC.right.f == neumann
        Dx[:,end] .= H[:,end] .* BC.right.val 
    else
        Dx[:,end] .= BC.right.val
    end
    if BC.top.f == neumann
        Dy[end,:] .= H[end,:] .* BC.top.val 
    else
        Dy[end,:] .= BC.top.val
    end

    return Dx, Dy
end

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

@inline function get_capacities(cap, II, Δ)
    ret = (cap[II,1]*Δ, cap[II,2]*Δ, cap[II,3]*Δ, cap[II,4]*Δ, cap[II,6]*Δ, cap[II,7]*Δ,
           cap[II,8]*Δ^2 + eps(0.01), cap[II,9]*Δ^2 + eps(0.01), cap[II,10]*Δ^2 + eps(0.01), cap[II,11]*Δ^2 + eps(0.01))
    return ret
end

@inline function set_lapl_bnd!(::Dirichlet, ::Dirichlet, L, B, W, C1, C2, n, Δ, b_indices, b2)
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Neumann, L, B, W, C1, C2, n, Δ, b_indices, b2)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] += B[II]*Δ * (C1[II]*Δ - C2[II]*Δ) / (W[II]*Δ^2+eps(0.01))
    end
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Periodic, L, B, W, C1, C2, n, Δ, b_indices, b_periodic)
    @inbounds @threads for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] += B[II]*Δ / (W[II]*Δ^2+eps(0.01)) * B[JJ]*Δ 
    end
    return nothing
end

function laplacian!(::Dirichlet, L, B, Dx, Dy, cap, n, Δ, BC, inside, empty, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]*Δ
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6]*Δ - A3) * Dx[δx⁺(II)]

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]*Δ
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]*Δ) * Dx[δx⁻(II)]
        
        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]*Δ
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7]*Δ - A4) * Dy[δy⁺(II)]
        
        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]*Δ
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]*Δ) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[1:end-1], b_top[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]*Δ
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6]*Δ - A3) * Dx[δx⁺(II)]
    end
    @inbounds @threads for II in vcat(b_bottom[2:end], b_right, b_top[2:end])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]*Δ
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]*Δ) * Dx[δx⁻(II)]
    end
    @inbounds @threads for II in vcat(b_left[1:end-1], b_bottom, b_right[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]*Δ
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7]*Δ - A4) * Dy[δy⁺(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end], b_right[2:end], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]*Δ
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]*Δ) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] = -4.0
    end

    set_lapl_bnd!(dir, BC.left.t, L, cap[:,:,6], cap[:,:,8], cap[:,:,6], cap[:,:,1], n, Δ, b_left, b_right)
    set_lapl_bnd!(dir, BC.bottom.t, L, cap[:,:,7], cap[:,:,9], cap[:,:,7], cap[:,:,2], n, Δ, b_bottom, b_top)
    set_lapl_bnd!(dir, BC.right.t, L, cap[:,:,6], cap[:,:,10], cap[:,:,6], cap[:,:,3], n, Δ, b_right, b_left)
    set_lapl_bnd!(dir, BC.top.t, L, cap[:,:,7], cap[:,:,11], cap[:,:,7], cap[:,:,4], n, Δ, b_top, b_bottom)
    
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Neumann, L, A, W, n, Δ, b_indices, b2)
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Dirichlet, L, A, W, n, Δ, b_indices, b2)
    @error ("Not implemented yet.\nTry Neumann or Periodic in the outer BCs.")
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Periodic, L, A, W, n, Δ, b_indices, b_periodic)
    @inbounds @threads for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] += (A[II]*Δ)^2 / (W[II]*Δ^2+eps(0.01))
    end
    return nothing
end

function laplacian!(::Neumann, L, B, Nx, Ny, cap, n, Δ, BC, inside, empty, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds L[pII,pII+n] = A3^2 / W3
        @inbounds L[pII,pII-n] = A1^2 / W1
        @inbounds L[pII,pII+1] = A4^2 / W4
        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[1:end-1], b_top[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds L[pII,pII+n] = A3^2 / W3
    end
    @inbounds @threads for II in vcat(b_bottom[2:end], b_right, b_top[2:end])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds L[pII,pII-n] = A1^2 / W1
    end
    @inbounds @threads for II in vcat(b_left[1:end-1], b_bottom, b_right[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds L[pII,pII+1] = A4^2 / W4
    end
    @inbounds @threads for II in vcat(b_left[2:end], b_right[2:end], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] = -4.0
    end

    set_lapl_bnd!(neu, BC.left.t, L, cap[:,:,1], cap[:,:,8], n, Δ, b_left, b_right)
    set_lapl_bnd!(neu, BC.bottom.t, L, cap[:,:,2], cap[:,:,9], n, Δ, b_bottom, b_top)
    set_lapl_bnd!(neu, BC.right.t, L, cap[:,:,3], cap[:,:,10], n, Δ, b_right, b_left)
    set_lapl_bnd!(neu, BC.top.t, L, cap[:,:,4], cap[:,:,11], n, Δ, b_top, b_bottom)

    return nothing
end

function crank_nicolson!(L, A, B, cap, τ, n, Δ, all_indices)
    @inbounds V = cap[:,:,5] .* Δ^2

    @inbounds A .= -L .* τ
    @inbounds B .= L .* τ
    @inbounds @threads for II in all_indices
        pII = lexicographic(II, n)
        @inbounds B[pII,pII] += V[II] * 2.
        @inbounds A[pII,pII] += V[II] * 2.
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
