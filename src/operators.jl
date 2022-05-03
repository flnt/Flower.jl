@inline function get_capacities(cap, II, Δ)
    @inbounds ret = (cap[II,1]*Δ, cap[II,2]*Δ, cap[II,3]*Δ, cap[II,4]*Δ, cap[II,6]*Δ, cap[II,7]*Δ,
           cap[II,8]*Δ^2 + eps(0.01), cap[II,9]*Δ^2 + eps(0.01), cap[II,10]*Δ^2 + eps(0.01), cap[II,11]*Δ^2 + eps(0.01))
    return ret
end

@inline function get_capacities_x(cap, II_1, II_2, Δ)
    @inbounds ret = cap[II_1,3] * Δ, cap[II_1,6] * Δ, cap[II_2,1] * Δ, cap[II_2,6] * Δ
    return ret
end

@inline function get_capacities_y(cap, II_1, II_2, Δ)
    @inbounds ret = cap[II_1,4] * Δ, cap[II_1,7] * Δ, cap[II_2,2] * Δ, cap[II_2,7] * Δ
    return ret
end

function set_bc_bnds(::Dirichlet, D, H, BC)
    Dx = copy(D)
    Dy = copy(D)

    if is_neumann(BC.left.t)
        @inbounds Dx[:,1] .= H[:,1] .* BC.left.val
    elseif is_dirichlet(BC.left.t)
        @inbounds Dx[:,1] .= BC.left.val
    end
    if is_neumann(BC.bottom.t)
        @inbounds Dy[1,:] .= H[1,:] .* BC.bottom.val 
    elseif is_dirichlet(BC.bottom.t)
        @inbounds Dy[1,:] .= BC.bottom.val
    end
    if is_neumann(BC.right.t)
        @inbounds Dx[:,end] .= H[:,end] .* BC.right.val 
    elseif is_dirichlet(BC.right.t)
        @inbounds Dx[:,end] .= BC.right.val
    end
    if is_neumann(BC.top.t)
        @inbounds Dy[end,:] .= H[end,:] .* BC.top.val 
    elseif is_dirichlet(BC.top.t)
        @inbounds Dy[end,:] .= BC.top.val
    end

    return Dx, Dy
end

@inline function set_lapl_bnd!(::Dirichlet, ::Dirichlet, L, B1, W, B2, n, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Neumann, L, B1, W, B2, n, Δ, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] += B1[II]*Δ * (B1[II]*Δ - B2[II]*Δ) / (W[II]*Δ^2+eps(0.01))
    end
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Periodic, L, B1, W, B2, n, Δ, b_indices, b_periodic)
    @inbounds @threads for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] += B1[II]*Δ / (W[II]*Δ^2+eps(0.01)) * B1[JJ]*Δ 
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

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]
    end
    @inbounds @threads for II in vcat(b_left, b_bottom[1:end-1], b_top[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]*Δ
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6]*Δ - A3) * Dx[δx⁺(II)]
    end
    @inbounds @threads for II in vcat(b_bottom[2:end], b_right, b_top[2:end])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]*Δ
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]*Δ) * Dx[δx⁻(II)]
    end
    @inbounds @threads for II in vcat(b_left[1:end-1], b_bottom, b_right[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]*Δ
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7]*Δ - A4) * Dy[δy⁺(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end], b_right[2:end], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]*Δ
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]*Δ) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] = -4.0
    end

    @inbounds _A1 = cap[:,:,1]
    @inbounds _A2 = cap[:,:,2]
    @inbounds _A3 = cap[:,:,3]
    @inbounds _A4 = cap[:,:,4]
    @inbounds _B1 = cap[:,:,6]
    @inbounds _B2 = cap[:,:,7]
    @inbounds _W1 = cap[:,:,8]
    @inbounds _W2 = cap[:,:,9]
    @inbounds _W3 = cap[:,:,10]
    @inbounds _W4 = cap[:,:,11]

    set_lapl_bnd!(dir, BC.left.t, L, _B1, _W1, _A1, n, Δ, b_left, b_right)
    set_lapl_bnd!(dir, BC.bottom.t, L, _B2, _W2, _A2, n, Δ, b_bottom, b_top)
    set_lapl_bnd!(dir, BC.right.t, L, _B1, _W3, _A3, n, Δ, b_right, b_left)
    set_lapl_bnd!(dir, BC.top.t, L, _B2, _W4, _A4, n, Δ, b_top, b_bottom)
    
    return nothing
end

function set_bc_bnds(::Neumann, Nx, Ny, BC)
    if is_neumann(BC.left.t)
        @inbounds Nx[:,1] .= BC.left.val
    end
    if is_neumann(BC.bottom.t)
        @inbounds Ny[1,:] .= BC.bottom.val 
    end
    if is_neumann(BC.right.t)
        @inbounds Nx[:,end] .= BC.right.val 
    end
    if is_neumann(BC.top.t)
        @inbounds Ny[end,:] .= BC.top.val 
    end

    return Nx, Ny
end

@inline function set_lapl_bnd!(::Neumann, ::Dirichlet, L, A, W, n, Δ, b_indices, b_periodic)
    @error ("Not implemented yet.\nTry Neumann or Periodic in the outer BCs.")
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Neumann, L, A, W, n, Δ, b_indices, b_periodic)
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

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]
    end
    @inbounds @threads for II in vcat(b_left, b_bottom[1:end-1], b_top[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+n] = A3^2 / W3
    end
    @inbounds @threads for II in vcat(b_bottom[2:end], b_right, b_top[2:end])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-n] = A1^2 / W1
    end
    @inbounds @threads for II in vcat(b_left[1:end-1], b_bottom, b_right[1:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+1] = A4^2 / W4
    end
    @inbounds @threads for II in vcat(b_left[2:end], b_right[2:end], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] = -4.0
    end

    @inbounds _A1 = cap[:,:,1]
    @inbounds _A2 = cap[:,:,2]
    @inbounds _A3 = cap[:,:,3]
    @inbounds _A4 = cap[:,:,4]
    @inbounds _W1 = cap[:,:,8]
    @inbounds _W2 = cap[:,:,9]
    @inbounds _W3 = cap[:,:,10]
    @inbounds _W4 = cap[:,:,11]

    set_lapl_bnd!(neu, BC.left.t, L, _A1, _W1, n, Δ, b_left, b_right)
    set_lapl_bnd!(neu, BC.bottom.t, L, _A2, _W2, n, Δ, b_bottom, b_top)
    set_lapl_bnd!(neu, BC.right.t, L, _A3, _W3, n, Δ, b_right, b_left)
    set_lapl_bnd!(neu, BC.top.t, L, _A4, _W4, n, Δ, b_top, b_bottom)

    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Dirichlet, fun, O, B1, B2, np, nuv, Δ, b_indices)
    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Neumann, fun, O, B1, B2, np, nuv, Δ, b_indices)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, np)
        pJJ = lexicographic(fun(II), nuv)
        @inbounds O[pII,pJJ] -= B1[fun(II)]*Δ - B2[fun(II)]*Δ
    end
    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Periodic, fun, O, B1, B2, np, nuv, Δ, b_indices)
    return nothing
end

@inline function divergence_boundaries(::Dirichlet, Ox, Oy, Dx, Dy, cap, capu, capv, n, Δ, BCu, BCv, b_left, b_bottom, b_right, b_top)
    Ax = hcat(cap[:,:,1], cap[:,end,3])
    Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

    mask = findall(iszero, Ax)
    capx = copy(capu)
    capx[mask,6] .= 0.
    mask = findall(iszero, Ay)
    capy = copy(capv)
    capy[mask,7] .= 0.

    set_div_bnd!(dir, BCu.left.t, x->x, Ox, capx[:,:,3], capx[:,:,6], n, n, Δ, b_left)
    set_div_bnd!(dir, BCv.bottom.t, x->x, Oy, capy[:,:,4], capy[:,:,7], n, n+1, Δ, b_bottom)
    set_div_bnd!(dir, BCu.right.t, δx⁺, Ox, capx[:,:,6], capx[:,:,1], n, n, Δ, b_right)
    set_div_bnd!(dir, BCv.top.t, δy⁺, Oy, capy[:,:,7], capy[:,:,2], n, n+1, Δ, b_top)
end

function divergence!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, cap, capu, capv, n, Δ, all_indices)
    Ax = hcat(cap[:,:,1], cap[:,end,3])
    Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

    mask = findall(iszero, Ax)
    capx = copy(capu)
    capx[mask,6] .= 0.
    mask = findall(iszero, Ay)
    capy = copy(capv)
    capy[mask,7] .= 0.

    Bx .= 0.0
    By .= 0.0
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        A3_1, B1_1, A1_2, B1_2 = get_capacities_x(capx, II, δx⁺(II), Δ)

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        A4_1, B2_1, A2_2, B2_2 = get_capacities_y(capy, II, δy⁺(II), Δ)
        
        @inbounds Ox[pII_1,pII_1] = -B1_1
        @inbounds Ox[pII_1,pII_2] = B1_2

        @inbounds Oy[pII_1,pJJ_1] = -B2_1
        @inbounds Oy[pII_1,pJJ_2] = B2_2

        @inbounds Bx[pII_1] += -(B1_2 - A1_2) * Dx[δx⁺(II)]
        @inbounds Bx[pII_1] += -(A3_1 - B1_1) * Dx[II]
        @inbounds By[pII_1] += -(B2_2 - A2_2) * Dy[δy⁺(II)]
        @inbounds By[pII_1] += -(A4_1 - B2_1) * Dy[II]
    end

    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Dirichlet, fun, O, B, nuv, np, Δ, b_indices, b_periodic)
    @error ("Not implemented yet.\nTry Neumann or Periodic in the outer BCs.")
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Neumann, fun, O, B, nuv, np, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Periodic, fun, O, B, nuv, np, Δ, b_indices, b_periodic)
    @inbounds @threads for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(JJ, np)
        @inbounds O[pII,pJJ] += B[II]
    end
    return nothing
end

function gradient!(::Neumann, Ox, Oy, Divx, Divy, cap, capu, capv, n, Δ, BC, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
    Ax = hcat(cap[:,:,1], cap[:,end,3])
    Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

    mask = findall(iszero, Ax)
    capx = copy(capu)
    capx[mask,6] .= 0.
    mask = findall(iszero, Ay)
    capy = copy(capv)
    capy[mask,7] .= 0.

    Ox .= -Divx'
    Oy .= -Divy'

    set_grad_bnd!(neu, BC.left.t, x->x, Ox, -capx[:,:,6], n, n, Δ, b_left_u, b_right_p)
    set_grad_bnd!(neu, BC.bottom.t, x->x, Oy, -capy[:,:,7], n+1, n, Δ, b_bottom_v, b_top_p)
    set_grad_bnd!(neu, BC.right.t, δx⁺, Ox, capx[:,:,6], n, n, Δ, b_right_u, b_left_p)
    set_grad_bnd!(neu, BC.top.t, δy⁺, Oy, capy[:,:,7], n+1, n, Δ, b_top_v, b_bottom_p)

    return nothing
end