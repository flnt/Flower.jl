function empty_laplacian(grid, O, empty, MIXED)
    @inbounds @threads for II in vcat(empty, MIXED)
        pII = lexicographic(II, grid.ny)
        if (sum(abs.(O[:,pII]))-4.0) <= 1e-8
            O[pII,pII] = 0.0
        end
    end
end

@inline function get_capacities(cap, II)
    @inbounds ret = (cap[II,1], cap[II,2], cap[II,3], cap[II,4], cap[II,6], cap[II,7],
           cap[II,8] + eps(0.01), cap[II,9] + eps(0.01), cap[II,10] + eps(0.01), cap[II,11] + eps(0.01))
    return ret
end

@inline function get_capacities_x(cap, II)
    @inbounds ret = cap[II,1], cap[II,3], cap[II,6]
    return ret
end

@inline function get_capacities_y(cap, II)
    @inbounds ret = cap[II,2], cap[II,4], cap[II,7]
    return ret
end

@inline function get_capacities_convection(cap, II)
    @inbounds ret = cap[II,1], cap[II,2], cap[II,3], cap[II,4], cap[II,6], cap[II,7]
    return ret
end

function set_bc_bnds(::Dirichlet, D, H, BC)
    Dx = copy(D)
    Dy = copy(D)

    if is_neumann(BC.left)
        @inbounds Dx[:,1] .= H[:,1] .* BC.left.val
    elseif is_dirichlet(BC.left)
        @inbounds Dx[:,1] .= BC.left.val
    end
    if is_neumann(BC.bottom)
        @inbounds Dy[1,:] .= H[1,:] .* BC.bottom.val 
    elseif is_dirichlet(BC.bottom)
        @inbounds Dy[1,:] .= BC.bottom.val
    end
    if is_neumann(BC.right)
        @inbounds Dx[:,end] .= H[:,end] .* BC.right.val 
    elseif is_dirichlet(BC.right)
        @inbounds Dx[:,end] .= BC.right.val
    end
    if is_neumann(BC.top)
        @inbounds Dy[end,:] .= H[end,:] .* BC.top.val 
    elseif is_dirichlet(BC.top)
        @inbounds Dy[end,:] .= BC.top.val
    end

    return Dx, Dy
end

@inline function set_lapl_bnd!(::Dirichlet, ::Dirichlet, L, B1, W, B2, n, b_indices, b_periodic)
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Neumann, L, B1, W, B2, n, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] += B1[II] * (B1[II] - B2[II]) / (W[II]+eps(0.01))
    end
    return nothing
end

@inline function set_lapl_bnd!(::Dirichlet, ::Periodic, L, B1, W, B2, n, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] = B1[II] / (W[II]+eps(0.01)) * B1[JJ]
        if abs(L[pII,pJJ]) < 1e-8
            L[pII,pII] = -4.0
        end
    end
    return nothing
end

function laplacian!(::Dirichlet, L, B, Dx, Dy, cap, n, BC, inside, empty, MIXED, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)
        
        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6] - A3) * Dx[δx⁺(II)]

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]) * Dx[δx⁻(II)]
        
        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7] - A4) * Dy[δy⁺(II)]
        
        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII] = -B1 * (B1/W3 + B1/W1) - B2 * (B2/W4 + B2/W2)

        @inbounds B[pII] += -B1 / W3 * (A3 - B1) * Dx[II]
        @inbounds B[pII] += B1 / W1 * (B1 - A1) * Dx[II]
        @inbounds B[pII] += -B2 / W4 * (A4 - B2) * Dy[II]
        @inbounds B[pII] += B2 / W2 * (B2 - A2) * Dy[II]
    end
    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6] - A3) * Dx[δx⁺(II)]
    end
    @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]) * Dx[δx⁻(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7] - A4) * Dy[δy⁺(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        if sum(abs.(L[:,pII])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
        end
    end
    @inbounds @threads for II in MIXED
        pII = lexicographic(II, n)
        if sum(abs.(L[:,pII])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
        end
    end

    @inbounds _A1 = @view cap[:,:,1]
    @inbounds _A2 = @view cap[:,:,2]
    @inbounds _A3 = @view cap[:,:,3]
    @inbounds _A4 = @view cap[:,:,4]
    @inbounds _B1 = @view cap[:,:,6]
    @inbounds _B2 = @view cap[:,:,7]
    @inbounds _W1 = @view cap[:,:,8]
    @inbounds _W2 = @view cap[:,:,9]
    @inbounds _W3 = @view cap[:,:,10]
    @inbounds _W4 = @view cap[:,:,11]

    set_lapl_bnd!(dir, BC.left, L, _B1, _W1, _A1, n, b_left, b_right)
    set_lapl_bnd!(dir, BC.bottom, L, _B2, _W2, _A2, n, b_bottom, b_top)
    set_lapl_bnd!(dir, BC.right, L, _B1, _W3, _A3, n, b_right, b_left)
    set_lapl_bnd!(dir, BC.top, L, _B2, _W4, _A4, n, b_top, b_bottom)
    
    return nothing
end

function set_bc_bnds(::Neumann, Nx, Ny, BC, dx, dy)
    if is_neumann(BC.left)
        @inbounds Nx[:,1] .= BC.left.val
    elseif is_dirichlet(BC.left)
        @inbounds Nx[:,1] .= -BC.left.val ./ dx[:,1] .* 2.0
    end
    if is_neumann(BC.bottom)
        @inbounds Ny[1,:] .= BC.bottom.val
    elseif is_dirichlet(BC.bottom)
        @inbounds Ny[1,:] .= -BC.bottom.val ./ dy[1,:] .* 2.0
    end
    if is_neumann(BC.right)
        @inbounds Nx[:,end] .= BC.right.val
    elseif is_dirichlet(BC.right)
        @inbounds Nx[:,end] .= BC.right.val ./ dx[:,end] .* 2.0
    end
    if is_neumann(BC.top)
        @inbounds Ny[end,:] .= BC.top.val
    elseif is_dirichlet(BC.top)
        @inbounds Ny[end,:] .= BC.top.val ./ dy[end,:] .* 2.0
    end

    return Nx, Ny
end

@inline function set_lapl_bnd!(::Neumann, ::Dirichlet, L, A, A1, A3, B1, dx, W, n, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds L[pII,pII] += (A3[II] - A1[II]) / dx[II] * 2.0
    end
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Neumann, L, A, A1, A3, B1, dx, W, n, b_indices, b_periodic)
    return nothing
end

@inline function set_lapl_bnd!(::Neumann, ::Periodic, L, A, A1, A3, B1, dx, W, n, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] = A[II]^2 / (W[II]+eps(0.01))
        # if abs(L[pII,pJJ]) < 1e-8
        if sum(abs.(L[pII,:])) <= 1e-10
            L[pII,pII] = -4.0
        end
    end
    return nothing
end

function laplacian!(::Neumann, L, B, Nx, Ny, HNx, HNy, cap, dx, dy, n, BC, inside, empty, MIXED, ns_vec, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds B[pII] += -A3 / W3 * ((cap[δx⁺(II),6] - A3) * HNx[δx⁺(II)] + (A3 - B1) * HNx[II])
        @inbounds B[pII] += A1 / W1 * ((B1 - A1) * HNx[II] + (A1 - cap[δx⁻(II),6]) * HNx[δx⁻(II)])
        @inbounds B[pII] += -A4 / W4 * ((cap[δy⁺(II),7] - A4) * HNy[δy⁺(II)] + (A4 - B2) * HNy[II])
        @inbounds B[pII] += A2 / W2 * ((B2 - A2) * HNy[II] + (A2 - cap[δy⁻(II),7]) * HNy[δy⁻(II)])

        @inbounds L[pII,pII+n] = A3^2 / W3
        @inbounds L[pII,pII-n] = A1^2 / W1
        @inbounds L[pII,pII+1] = A4^2 / W4
        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)
        
        @inbounds L[pII,pII] = -(A1^2 / W1 + A2^2 / W2 + A3^2 / W3 + A4^2 / W4)

        @inbounds B[pII] += -(A3 - B1) * Nx[II] - (B1 - A1) * Nx[II]
        @inbounds B[pII] += -(A4 - B2) * Ny[II] - (B2 - A2) * Ny[II]

        @inbounds B[pII] += -A3 / W3 * (A3 - B1) * HNx[II]
        @inbounds B[pII] += A1 / W1 * (B1 - A1) * HNx[II]
        @inbounds B[pII] += -A4 / W4 * (A4 - B2) * HNy[II]
        @inbounds B[pII] += A2 / W2 * (B2 - A2) * HNy[II]
    end
    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds B[pII] += -A3 / W3 * (cap[δx⁺(II),6] - A3) * HNx[δx⁺(II)]

        @inbounds L[pII,pII+n] = A3^2 / W3
    end
    @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds B[pII] += A1 / W1 * (A1 - cap[δx⁻(II),6]) * HNx[δx⁻(II)]

        @inbounds L[pII,pII-n] = A1^2 / W1
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds B[pII] += -A4 / W4 * (cap[δy⁺(II),7] - A4) * HNy[δy⁺(II)]

        @inbounds L[pII,pII+1] = A4^2 / W4
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II)

        @inbounds B[pII] += A2 / W2 * (A2 - cap[δy⁻(II),7]) * HNy[δy⁻(II)]

        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    ns_vec .= 1.
    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        if sum(abs.(L[:,pII])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
            @inbounds ns_vec[pII] = 0.
        end
    end
    @inbounds @threads for II in MIXED
        pII = lexicographic(II, n)
        if sum(abs.(L[:,pII])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
            @inbounds ns_vec[pII] = 0.
        end
    end
    ns_vec ./= norm(ns_vec)

    @inbounds _A1 = @view cap[:,:,1]
    @inbounds _A2 = @view cap[:,:,2]
    @inbounds _A3 = @view cap[:,:,3]
    @inbounds _A4 = @view cap[:,:,4]
    @inbounds _B1 = @view cap[:,:,6]
    @inbounds _B2 = @view cap[:,:,7]
    @inbounds _W1 = @view cap[:,:,8]
    @inbounds _W2 = @view cap[:,:,9]
    @inbounds _W3 = @view cap[:,:,10]
    @inbounds _W4 = @view cap[:,:,11]

    set_lapl_bnd!(neu, BC.left, L, _A1, _A3, _A1, _B1, dx, _W1, n, b_left, b_right)
    set_lapl_bnd!(neu, BC.bottom, L, _A2, _A4, _A2, _B2, dy, _W2, n, b_bottom, b_top)
    set_lapl_bnd!(neu, BC.right, L, _A3, _A1, _A3, _B1, dx, _W3, n, b_right, b_left)
    set_lapl_bnd!(neu, BC.top, L, _A4, _A2, _A4, _B2, dy, _W4, n, b_top, b_bottom)

    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Dirichlet, fun, O, B1, B2, np, nuv, b_indices)
    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Neumann, fun, O, B1, B2, np, nuv, b_indices)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, np)
        pJJ = lexicographic(fun(II), nuv)
        @inbounds O[pII,pJJ] += -(B1[II] - B2[II])
    end
    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Periodic, fun, O, B1, B2, np, nuv, b_indices)
    return nothing
end

@inline function divergence_boundaries(::Dirichlet, Ox, Oy, Dx, Dy, cap, n, BCu, BCv, b_left, b_bottom, b_right, b_top)
    set_div_bnd!(dir, BCu.left, x->x, Ox, view(cap,:,:,6), view(cap,:,:,1), n, n, b_left)
    set_div_bnd!(dir, BCv.bottom, x->x, Oy, view(cap,:,:,7), view(cap,:,:,2), n, n+1, b_bottom)
    set_div_bnd!(dir, BCu.right, δx⁺, Ox, view(cap,:,:,3), view(cap,:,:,6), n, n, b_right)
    set_div_bnd!(dir, BCv.top, δy⁺, Oy, view(cap,:,:,4), view(cap,:,:,7), n, n+1, b_top)
end

function divergence!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, cap, n, all_indices)
    Bx .= 0.0
    By .= 0.0
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        A1, A3, B1 = get_capacities_x(cap, II)

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        A2, A4, B2 = get_capacities_y(cap, II)
        
        @inbounds Ox[pII_1,pII_1] = -A1
        @inbounds Ox[pII_1,pII_2] = A3

        @inbounds Oy[pII_1,pJJ_1] = -A2
        @inbounds Oy[pII_1,pJJ_2] = A4

        @inbounds Bx[pII_1] += -(A3 - B1) * Dx[δx⁺(II)]
        @inbounds Bx[pII_1] += -(B1 - A1) * Dx[II]
        @inbounds By[pII_1] += -(A4 - B2) * Dy[δy⁺(II)]
        @inbounds By[pII_1] += -(B2 - A2) * Dy[II]
    end

    return nothing
end

function set_bc_bnds(::Neumann, HNx, HNy, BC, dx, dy, H)
    if is_neumann(BC.left)
        @inbounds HNx[:,1] .= H[:,1] .* BC.left.val
    elseif is_dirichlet(BC.left)
        @inbounds HNx[:,1] .= BC.left.val
    end
    if is_neumann(BC.bottom)
        @inbounds HNy[1,:] .= H[1,:] .* BC.bottom.val
    elseif is_dirichlet(BC.bottom)
        @inbounds HNy[1,:] .= BC.bottom.val
    end
    if is_neumann(BC.right)
        @inbounds HNx[:,end] .= H[:,end] .* BC.right.val
    elseif is_dirichlet(BC.right)
        @inbounds HNx[:,end] .= BC.right.val
    end
    if is_neumann(BC.top)
        @inbounds HNy[end,:] .= H[end,:] .* BC.top.val
    elseif is_dirichlet(BC.top)
        @inbounds HNy[end,:] .= BC.top.val
    end

    return HNx, HNy
end

@inline function set_grad_bnd!(::Neumann, ::Dirichlet, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(fun(II), np)
        @inbounds O[pII,pJJ] += B1[fun(II)] - A1[fun(II)]
    end
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Neumann, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Periodic, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(JJ, np)
        @inbounds O[pII,pJJ] = B[fun(II)]
    end
    return nothing
end

function gradient!(::Neumann, Ox, Oy, Bx, By, HNx, HNy, Divx, Divy, dcap, n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
    mat_T_op!(Ox, Divx, x->-x)
    mat_T_op!(Oy, Divy, x->-x)

    Bx .= 0.
    By .= 0.
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * HNx[II] + (A1 - dcap[δx⁻(II),6]) * HNx[δx⁻(II)])
        @inbounds By[pJJ] += -((B2 - A2) * HNy[II] + (A2 - dcap[δy⁻(II),7]) * HNy[δy⁻(II)])
    end

    @inbounds @threads for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * HNx[II])
    end
    @inbounds @threads for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((A3 - B1) * HNx[II])
    end
    @inbounds @threads for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * HNx[II] + (A1 - dcap[δx⁻(II),6]) * HNx[δx⁻(II)])
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((B2 - A2) * HNy[II])
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pJJ = lexicographic(δy⁺(II), n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((A4 - B2) * HNy[II])
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((B2 - A2) * HNy[II] + (A2 - dcap[δy⁻(II),7]) * HNy[δy⁻(II)])
    end

    @inbounds _A1 = @view dcap[:,:,1]
    @inbounds _A2 = @view dcap[:,:,2]
    @inbounds _A3 = @view dcap[:,:,3]
    @inbounds _A4 = @view dcap[:,:,4]
    @inbounds _B1 = @view dcap[:,:,6]
    @inbounds _B2 = @view dcap[:,:,7]

    set_grad_bnd!(neu, BC.left, x->x, Ox, -_A1, _B1, _A1, n, n, b_left_u, b_right_p)
    set_grad_bnd!(neu, BC.bottom, x->x, Oy, -_A2, _B2, _A2, n+1, n, b_bottom_v, b_top_p)
    set_grad_bnd!(neu, BC.right, δx⁻, Ox, _A3, _A3, _B1, n, n, b_right_u, b_left_p)
    set_grad_bnd!(neu, BC.top, δy⁻, Oy, _A4, _A4, _B2, n+1, n, b_top_v, b_bottom_p)

    return nothing
end

@inline function set_div_bnd!(::Neumann, ::Dirichlet, fun, O, B1, B2, np, nuv, b_indices)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, np)
        pJJ = lexicographic(fun(II), nuv)
        @inbounds O[pII,pJJ] += (B1[II] - B2[II])
    end

    return nothing
end

@inline function set_div_bnd!(::Neumann, ::Neumann, fun, O, B1, B2, np, nuv, b_indices)
    return nothing
end

@inline function set_div_bnd!(::Neumann, ::Periodic, fun, O, B1, B2, np, nuv, b_indices)
    return nothing
end

@inline function divergence_boundaries(::Neumann, Ox, Oy, Dx, Dy, cap, n, BCu, BCv, b_left, b_bottom, b_right, b_top)
    set_div_bnd!(neu, BCu.left, x->x, Ox, view(cap,:,:,6), view(cap,:,:,1), n, n, b_left)
    set_div_bnd!(neu, BCv.bottom, x->x, Oy, view(cap,:,:,7), view(cap,:,:,2), n, n+1, b_bottom)
    set_div_bnd!(neu, BCu.right, δx⁺, Ox, view(cap,:,:,3), view(cap,:,:,6), n, n, b_right)
    set_div_bnd!(neu, BCv.top, δy⁺, Oy, view(cap,:,:,4), view(cap,:,:,7), n, n+1, b_top)
end

function divergence!(::Neumann, Ox, Oy, Bx, By, NHx, NHy, cap, n, all_indices)
    Bx .= 0.0
    By .= 0.0
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        
        A1, A3, B1 = get_capacities_x(cap, II)
        A2, A4, B2 = get_capacities_y(cap, II)

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        
        @inbounds Ox[pII_1,pII_1] = -B1
        @inbounds Ox[pII_1,pII_2] = B1

        @inbounds Oy[pII_1,pJJ_1] = -B2
        @inbounds Oy[pII_1,pJJ_2] = B2

        @inbounds Bx[pII_1] += -(A3 - B1) * NHx[δx⁺(II)]
        @inbounds Bx[pII_1] += -(B1 - A1) * NHx[II]
        @inbounds By[pII_1] += -(A4 - B2) * NHy[δy⁺(II)]
        @inbounds By[pII_1] += -(B2 - A2) * NHy[II]
    end

    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Dirichlet, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Neumann, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(fun(II), np)
        @inbounds O[pII,pJJ] += -(B1[fun(II)] - A1[fun(II)])
    end

    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Periodic, fun, O, B, B1, A1, nuv, np, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(JJ, np)
        @inbounds O[pII,pJJ] = B[JJ]
    end
    return nothing
end

function gradient!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, Divx, Divy, dcap, n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
    mat_T_op!(Ox, Divx, x->-x)
    mat_T_op!(Oy, Divy, x->-x)

    Bx .= 0.
    By .= 0.
    @inbounds @threads for II in @view all_indices[2:end,2:end]
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * Dx[II] + (A1 - dcap[δx⁻(II),6]) * Dx[δx⁻(II)])
        @inbounds By[pJJ] += -((B2 - A2) * Dy[II] + (A2 - dcap[δy⁻(II),7]) * Dy[δy⁻(II)])
    end

    @inbounds @threads for II in @view all_indices[:,1]
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * Dx[II])
    end
    @inbounds @threads for II in @view all_indices[:,end]
        pII = lexicographic(δx⁺(II), n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((A3 - B1) * Dx[II])
    end
    @inbounds @threads for II in @views vcat(all_indices[1,2:end], all_indices[end,2:end])
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds Bx[pII] += -((B1 - A1) * Dx[II] + (A1 - dcap[δx⁻(II),6]) * Dx[δx⁻(II)])
    end

    @inbounds @threads for II in @view all_indices[1,:]
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((B2 - A2) * Dy[II])
    end
    @inbounds @threads for II in @view all_indices[end,:]
        pJJ = lexicographic(δy⁺(II), n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((A4 - B2) * Dy[II])
    end
    @inbounds @threads for II in @views vcat(all_indices[2:end,1], all_indices[2:end,end])
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(dcap, II)
        
        @inbounds By[pJJ] += -((B2 - A2) * Dy[II] + (A2 - dcap[δy⁻(II),7]) * Dy[δy⁻(II)])
    end

    @inbounds _A1 = @view dcap[:,:,1]
    @inbounds _A2 = @view dcap[:,:,2]
    @inbounds _A3 = @view dcap[:,:,3]
    @inbounds _A4 = @view dcap[:,:,4]
    @inbounds _B1 = @view dcap[:,:,6]
    @inbounds _B2 = @view dcap[:,:,7]

    set_grad_bnd!(dir, BC.left, x->x, Ox, -_B1, _B1, _A1, n, n, b_left_u, b_right_p)
    set_grad_bnd!(dir, BC.bottom, x->x, Oy, -_B2, _B2, _A2, n+1, n, b_bottom_v, b_top_p)
    set_grad_bnd!(dir, BC.right, δx⁻, Ox, _B1, _A3, _B1, n, n, b_right_u, b_left_p)
    set_grad_bnd!(dir, BC.top, δy⁻, Oy, _B2, _A4, _B2, n+1, n, b_top_v, b_bottom_p)

    return nothing
end

function uv_to_p!(Ox, Oy, cap, dx, dy, n, all_indices)
    @inbounds @threads for II in all_indices
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)

        @inbounds Ox[pII,pII] = 0.5
        @inbounds Ox[pII,pII+n] = 0.5

        @inbounds Oy[pII,pJJ] = 0.5
        @inbounds Oy[pII,pJJ+1] = 0.5
    end

    return nothing
end

function harmonic_average(W4, W3)
    # Harmonic average of volume capacities
    if W4 < 1e-8 || W3 < 1e-8
        Ŵ =  W4 + W3
    else
        Ŵ = 2 * W4 * W3 / (W4 + W3)
    end

    return Ŵ
end

# Didn't implement the inhomogeneous BC term for the moment
# since we only need this operator to compute forces inside
# the domain and we have the no-slip condition at the wall
function strain_rate!(::Dirichlet, O11, O12_x, O12_y, O22, cap_x, cap_y, n, all_indices, inside)
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)

        JJ = δx⁺(II)
        pJJ = lexicographic(JJ, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_x, II)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_x, JJ)

        @inbounds O11[pII,pII] = -B1_1 / W3_1
        @inbounds O11[pII,pJJ] = B1_2 / W3_1

        JJ = δy⁺(II)
        _pII = lexicographic(II, n+1)
        pJJ = lexicographic(JJ, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_y, II)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_y, JJ)

        @inbounds O22[pII,_pII] = -B2_1 / W4_1
        @inbounds O22[pII,pJJ] = B2_2 / W4_1
    end
    @inbounds @threads for II in all_indices[1:end-1,2:end]
        pII = lexicographic(II, n)

        JJ_1 = II
        JJ_2 = δy⁺(II)
        pJJ_1 = lexicographic(JJ_1, n)
        pJJ_2 = lexicographic(JJ_2, n)
        A1_1_x, A2_1_x, A3_1_x, A4_1_x, B1_1_x, B2_1_x, W1_1_x, W2_1_x, W3_1_x, W4_1_x = get_capacities(cap_x, JJ_1)
        A1_2_x, A2_2_x, A3_2_x, A4_2_x, B1_2_x, B2_2_x, W1_2_x, W2_2_x, W3_2_x, W4_2_x = get_capacities(cap_x, JJ_2)

        @inbounds O12_x[pII,pJJ_1] = -B2_1_x / 2W4_1_x
        @inbounds O12_x[pII,pJJ_2] = B2_2_x / 2W4_1_x

        KK_1 = δy⁺(δx⁻(II))
        KK_2 = δy⁺(II)
        pKK_1 = lexicographic(KK_1, n+1)
        pKK_2 = lexicographic(KK_2, n+1)
        A1_1_y, A2_1_y, A3_1_y, A4_1_y, B1_1_y, B2_1_y, W1_1_y, W2_1_y, W3_1_y, W4_1_y = get_capacities(cap_y, KK_1)
        A1_2_y, A2_2_y, A3_2_y, A4_2_y, B1_2_y, B2_2_y, W1_2_y, W2_2_y, W3_2_y, W4_2_y = get_capacities(cap_y, KK_2)

        Ŵ = harmonic_average(W4_1_x, W3_1_y)

        @inbounds O12_y[pII,pKK_1] = -B1_1_y / 2Ŵ
        @inbounds O12_y[pII,pKK_2] = B1_2_y / 2Ŵ
    end

    return nothing
end

function set_bc_bnds(::Dirichlet, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
    Dx = copy(Du)
    Dy = copy(Dv)

    if is_neumann(BC_u.left)
        @inbounds Dx[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
    elseif is_dirichlet(BC_u.left)
        @inbounds Dx[:,1] .= BC_u.left.val
    elseif is_periodic(BC_u.left)
        @inbounds Dx[:,1] .= u[:,end]
    end
    if is_neumann(BC_v.bottom)
        @inbounds Dy[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
    elseif is_dirichlet(BC_v.bottom)
        @inbounds Dy[1,:] .= BC_v.bottom.val
    elseif is_periodic(BC_v.bottom)
        @inbounds Dy[1,:] .= v[end,:]
    end
    if is_neumann(BC_u.right)
        @inbounds Dx[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
    elseif is_dirichlet(BC_u.right)
        @inbounds Dx[:,end] .= BC_u.right.val
    elseif is_periodic(BC_u.right)
        @inbounds Dx[:,end] .= u[:,1]
    end
    if is_neumann(BC_v.top)
        @inbounds Dy[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
    elseif is_dirichlet(BC_v.top)
        @inbounds Dy[end,:] .= BC_v.top.val
    elseif is_periodic(BC_v.top)
        @inbounds Dy[end,:] .= v[1,:]
    end

    return Dx, Dy
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
    return nothing
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Neumann, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
    @inbounds for II in b_indices
        pII = lexicographic(II, n)
        @inbounds O[pII,pII] += -0.5 * ((A2[II] - B1[II]) * D[fun(II)] + (B1[II] - A1[II]) * D[II])
    end
    return nothing
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Periodic, O, fun, A1, A2, B1, D, n, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pII,pJJ] += -0.5 * ((A2[II] - B1[II]) * D[fun(II)] + (B1[II] - A1[II]) * D[II])
    end
    return nothing
end

function scalar_convection!(::Dirichlet, O, B, u, v, Dx, Dy, Du, Dv, cap, n, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    O .= 0.0
    @inbounds for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
        u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

        @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)
        @inbounds O[pII,pII+n] = 0.5 * A3 * u3
        @inbounds O[pII,pII-n] = -0.5 * A1 * u1
        @inbounds O[pII,pII+1] = 0.5 * A4 * v4
        @inbounds O[pII,pII-1] = -0.5 * A2 * v2

        @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

        @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
    end

    @inbounds for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
        u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

        @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)

        @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

        @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
    end
    @inbounds for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

        @inbounds O[pII,pII+n] = 0.5 * A3 * u[δx⁺(II)]
    end
    @inbounds for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

        @inbounds O[pII,pII-n] = -0.5 * A1 * u[II]
    end
    @inbounds for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

        @inbounds O[pII,pII+1] = 0.5 * A4 * v[δy⁺(II)]
    end
    @inbounds for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)

        @inbounds O[pII,pII-1] = -0.5 * A2 * v[II]
    end

    if is_periodic(BC.left) && is_periodic(BC.right)
        @inbounds for (II,JJ) in zip(b_right, b_left)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
            @inbounds O[pII,pJJ] = 0.5 * A3 * u[δx⁺(II)]
        end
        @inbounds for (II,JJ) in zip(b_left, b_right)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
            @inbounds O[pII,pJJ] = -0.5 * A1 * u[II]
        end
    end
    if is_periodic(BC.bottom) && is_periodic(BC.top)
        @inbounds for (II,JJ) in zip(b_top, b_bottom)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
            @inbounds O[pII,pJJ] = 0.5 * A4 * v[δy⁺(II)]
        end
        @inbounds for (II,JJ) in zip(b_bottom, b_top)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II)
    
            @inbounds O[pII,pJJ] = -0.5 * A2 * v[II]
        end
    end

    @inbounds _A1 = @view cap[:,:,1]
    @inbounds _A2 = @view cap[:,:,2]
    @inbounds _A3 = @view cap[:,:,3]
    @inbounds _A4 = @view cap[:,:,4]
    @inbounds _B1 = @view cap[:,:,6]
    @inbounds _B2 = @view cap[:,:,7]

    set_sca_conv_bnd!(dir, BC.left, O, δx⁺, _A1, _A3, _B1, Du, n, b_left, b_right)
    set_sca_conv_bnd!(dir, BC.bottom, O, δy⁺, _A2, _A4, _B2, Dv, n, b_bottom, b_top)
    set_sca_conv_bnd!(dir, BC.right, O, δx⁺, _A1, _A3, _B1, Du, n, b_right, b_left)
    set_sca_conv_bnd!(dir, BC.top, O, δy⁺, _A2, _A4, _B2, Dv, n, b_top, b_bottom)

    return nothing
end

function set_bc_bnds(::Dirichlet, ::Union{Type{GridFCx},Type{GridFCy}}, Du_x, Du_y, Dv_x, Dv_y, Hu, Hv, u, v, BC_u, BC_v)

    if is_neumann(BC_u.left)
        @inbounds Du_x[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
        # @inbounds Du_x[:,2] .= u[:,2]
    elseif is_dirichlet(BC_u.left)
        # @inbounds Du_x[:,1] .= BC_u.left.val
        # @inbounds Du_x[:,2] .= BC_u.left.val
    elseif is_periodic(BC_u.left)
        @inbounds Du_x[:,1] .= u[:,end]
        @inbounds Du_x[:,2] .= u[:,2]
    end
    if is_neumann(BC_u.bottom)
        @inbounds Du_y[1,:] .= u[1,:] .+ Hu[1,:] .* BC_u.bottom.val
        # @inbounds Du_y[2,:] .= u[2,:]
    elseif is_dirichlet(BC_u.bottom)
        # @inbounds Du_y[1,:] .= BC_u.bottom.val
        # @inbounds Du_y[2,:] .= BC_u.bottom.val
    elseif is_periodic(BC_u.bottom)
        @inbounds Du_y[1,:] .= u[end,:]
        @inbounds Du_y[2,:] .= u[2,:]
    end
    if is_neumann(BC_u.right)
        @inbounds Du_x[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
        # @inbounds Du_x[:,end-1] .= u[:,end-1]
    elseif is_dirichlet(BC_u.right)
        # @inbounds Du_x[:,end] .= BC_u.right.val
        # @inbounds Du_x[:,end-1] .= BC_u.right.val
    elseif is_periodic(BC_u.right)
        @inbounds Du_x[:,end] .= u[:,1]
        @inbounds Du_x[:,end-1] .= u[:,end-1]
    end
    if is_neumann(BC_u.top)
        @inbounds Du_y[end,:] .= u[end,:] .+ Hu[end,:] .* BC_u.top.val
        # @inbounds Du_y[end-1,:] .= u[end-1,:]
    elseif is_dirichlet(BC_u.top)
        # @inbounds Du_y[end,:] .= BC_u.top.val
        # @inbounds Du_y[end-1,:] .= BC_u.top.val
    elseif is_periodic(BC_u.top)
        @inbounds Du_y[end,:] .= u[1,:]
        @inbounds Du_y[end-1,:] .= u[end-1,:]
    end

    if is_neumann(BC_v.left)
        @inbounds Dv_x[:,1] .= v[:,1] .+ Hv[:,1] .* BC_v.left.val 
        # @inbounds Dv_x[:,2] .= v[:,2]
    elseif is_dirichlet(BC_v.left)
        # @inbounds Dv_x[:,1] .= BC_v.left.val
        # @inbounds Dv_x[:,2] .= BC_v.left.val
    elseif is_periodic(BC_v.left)
        @inbounds Dv_x[:,1] .= v[:,end]
        @inbounds Dv_x[:,2] .= v[:,2]
    end
    if is_neumann(BC_v.bottom)
        @inbounds Dv_y[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
        # @inbounds Dv_y[2,:] .= v[2,:]
    elseif is_dirichlet(BC_v.bottom)
        # @inbounds Dv_y[1,:] .= BC_v.bottom.val
        # @inbounds Dv_y[2,:] .= BC_v.bottom.val
    elseif is_periodic(BC_v.bottom)
        @inbounds Dv_y[1,:] .= v[end,:]
        @inbounds Dv_y[2,:] .= v[2,:]
    end
    if is_neumann(BC_v.right)
        @inbounds Dv_x[:,end] .= v[:,end] .+ Hv[:,end] .* BC_v.right.val 
        # @inbounds Dv_x[:,end-1] .= v[:,end-1]
    elseif is_dirichlet(BC_v.right)
        # @inbounds Dv_x[:,end] .= BC_v.right.val
        # @inbounds Dv_x[:,end-1] .= BC_v.right.val
    elseif is_periodic(BC_v.right)
        @inbounds Dv_x[:,end] .= v[:,1]
        @inbounds Dv_x[:,end-1] .= v[:,end-1]
    end
    if is_neumann(BC_v.top)
        @inbounds Dv_y[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
        # @inbounds Dv_y[end-1,:] .= v[end-1,:]
    elseif is_dirichlet(BC_v.top)
        # @inbounds Dv_y[end,:] .= BC_v.top.val
        # @inbounds Dv_y[end-1,:] .= BC_v.top.val
    elseif is_periodic(BC_v.top)
        @inbounds Dv_y[end,:] .= v[1,:]
        @inbounds Dv_y[end-1,:] .= v[end-1,:]
    end

    return nothing
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, b_indices, b_periodic)
    return nothing
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Neumann, O, fun_cap, fun1, fun2, fun3, fun4, A3, B1, A1, A4, B2, A2, Du, Dv, n, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds O[pII,pII] += -0.25 * (A3[fun_cap(II)] - B1[fun_cap(II)]) * Du[fun1(II)]
        @inbounds O[pII,pII] += -0.25 * (B1[fun_cap(II)] - A1[fun_cap(II)]) * Du[fun2(II)]
        @inbounds O[pII,pII] += -0.25 * (A4[fun_cap(II)] - B2[fun_cap(II)]) * Dv[fun3(II)]
        @inbounds O[pII,pII] += -0.25 * (B2[fun_cap(II)] - A2[fun_cap(II)]) * Dv[fun4(II)]
    end
    return nothing
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Periodic, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, b_indices, b_periodic)
    @inbounds for (II,JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pII,pJJ] = -0.25 * (A1[fun_cap(II)] - B1[fun_cap(II)]) * D[fun1(II)]
        @inbounds O[pII,pJJ] = -0.25 * (B1[fun_cap(II)] - A2[fun_cap(II)]) * D[fun2(II)]
    end
    return nothing
end

# function fill_inside_conv!(::Type{GridFCx}, O, B, u, v, Du_x, Du_y, Dv_y, cap, ny, II)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)

#     Auim1, Aui, Auip1 = A1_1 * u[δx⁻(II)], A1_2 * u[II], A3_2 * u[δx⁺(II)]
#     Avim1jm1, Avim1jp1 = A2_1 * v[δx⁻(II)], A4_1 * v[δy⁺(δx⁻(II))]
#     Avip1jm1, Avip1jp1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

#     Au1 = 0.5 * (Auim1 + Aui)
#     Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
#     Au3 = 0.5 * (Aui + Auip1)
#     Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

#     @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
#     @inbounds O[pII,pII+ny] = 0.5 * Au3
#     @inbounds O[pII,pII-ny] = -0.5 * Au1
#     @inbounds O[pII,pII+1] = 0.5 * Au4
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - B1_1) * Du_x[II]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[δx⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

#     @inbounds B[pII] += -0.25 * Du_x[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Du_x[II] * (B1_2 - B1_1) * Du_x[II]
#     @inbounds B[pII] += -0.25 * Du_x[II] * (B1_1 - A1_1) * Du_x[δx⁻(II)]

#     @inbounds B[pII] += -0.25 * Du_y[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
#     @inbounds B[pII] += -0.25 * Du_y[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
#     @inbounds B[pII] += -0.25 * Du_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Du_y[II] * (B2_2 - A2_2) * Dv_y[II]

#     return nothing
# end

# function vec_convx_1!(II, O, B, u, Du, Dv, cap, ny)
#     pII = lexicographic(II, ny)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Aui, Auip1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

#     Au3 = 0.5 * (Aui + Auip1)

#     @inbounds O[pII,pII] += 0.5 * Au3
#     @inbounds O[pII,pII+ny] = 0.5 * Au3

#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

#     @inbounds B[pII] += -0.25 * Du[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Du[II] * (B1_2 - A1_2) * Du[II]

#     @inbounds B[pII] += -0.25 * Du[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Du[II] * (B2_2 - A2_2) * Dv[II]

#     return nothing
# end

# function vec_convx_2!(II, O, B, u, Du, Dv, cap, ny)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
#     Auim1, Aui = A1_1 * u[δx⁻(II)], A3_1 * u[II]

#     Au1 = 0.5 * (Auim1 + Aui)

#     @inbounds O[pII,pII] += -0.5 * Au1
#     @inbounds O[pII,pII-ny] = -0.5 * Au1

#     @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[II]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δx⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δx⁻(II)]

#     @inbounds B[pII] += -0.25 * Du[II] * (A3_1 - B1_1) * Du[II]
#     @inbounds B[pII] += -0.25 * Du[II] * (B1_1 - A1_1) * Du[δx⁻(II)]

#     @inbounds B[pII] += -0.25 * Du[II] * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
#     @inbounds B[pII] += -0.25 * Du[II] * (B2_1 - A2_1) * Dv[δx⁻(II)]

#     return nothing
# end

# function vec_convx_3!(II, O, v, cap, ny)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]
#     Avip1jp1 = A4_2 * v[δy⁺(II)]

#     Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

#     @inbounds O[pII,pII] += 0.5 * Au4
#     @inbounds O[pII,pII+1] = 0.5 * Au4

#     return nothing
# end

# function vec_convx_4!(II, O, v, cap, ny)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avim1jm1 = A2_1 * v[δx⁻(II)]
#     Avip1jm1 = A2_2 * v[II]

#     Au2 = 0.5 * (Avim1jm1 + Avip1jm1)

#     @inbounds O[pII,pII] += -0.5 * Au2
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     return nothing
# end

# function vec_convx_5!(II, O, v, cap, n, ny, BC)
#     pII = lexicographic(II, ny)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avip1jm1 = A2_2 * v[II]

#     Au2 = 0.5 * Avip1jm1

#     if is_periodic(BC.left)
#         JJ = II + CartesianIndex(0, n-1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#         Avim1jm1 = A2_1 * v[JJ]
#         Au2 += 0.5 * Avim1jm1
#     end

#     @inbounds O[pII,pII] += -0.5 * Au2
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     return nothing
# end

# function vec_convx_6!(II, O, v, cap, n, ny, BC)
#     pII = lexicographic(II, ny)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avip1jp1 = A4_2 * v[δy⁺(II)]

#     Au4 = 0.5 * Avip1jp1

#     if is_periodic(BC.left)
#         JJ = II + CartesianIndex(0, n-1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#         Avim1jp1 = A4_1 * v[δy⁺(JJ)]
#         Au4 += 0.5 * Avim1jp1
#     end

#     @inbounds O[pII,pII] += 0.5 * Au4
#     @inbounds O[pII,pII+1] = 0.5 * Au4

#     return nothing
# end

# function vec_convx_7!(II, O, v, cap, n, ny, BC)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
#     Avim1jm1 = A2_1 * v[δx⁻(II)]

#     Au2 = 0.5 * Avim1jm1

#     if is_periodic(BC.right)
#         JJ = II + CartesianIndex(0, -n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#         Avip1jm1 = A2_2 * v[JJ]
#         Au2 += 0.5 * Avip1jm1
#     end

#     @inbounds O[pII,pII] += -0.5 * Au2
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     return nothing
# end

# function vec_convx_8!(II, O, v, cap, n, ny, BC)
#     pII = lexicographic(II, ny)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
#     Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]

#     Au4 = 0.5 * Avim1jp1

#     if is_periodic(BC.right)
#         JJ = II + CartesianIndex(0, -n)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#         Avip1jp1 = A4_2 * v[δy⁺(JJ)]
#         Au4 += 0.5 * Avip1jp1
#     end

#     @inbounds O[pII,pII] += 0.5 * Au4
#     @inbounds O[pII,pII+1] = 0.5 * Au4

#     return nothing
# end

# function vector_convection!(::Dirichlet, ::Type{GridFCx}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
#     B .= 0.0
#     @inbounds @threads for II in inside
#         fill_inside_conv!(GridFCx, O, B, u, v, Du_x, Du_y, Dv_y, cap, ny, II)
#     end

#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, ny)
#         @inbounds O[pII,pII] = 0.0
#     end
#     bnds = (b_left, b_bottom[2:end-1], b_top[2:end-1])
#     bc = ((Du_x, Dv_x), (Du_y, Dv_y), (Du_y, Dv_y))
#     for (bnd, (Du, Dv)) in zip(bnds, bc)
#         @inbounds @threads for II in bnd
#             vec_convx_1!(II, O, B, u, Du, Dv, cap, ny)
#         end
#     end
#     bnds = (b_bottom[2:end-1], b_right, b_top[2:end-1])
#     bc = ((Du_y, Dv_y), (Du_x, Dv_x), (Du_y, Dv_y))
#     for (bnd, (Du, Dv)) in zip(bnds, bc)
#         @inbounds @threads for II in bnd
#             vec_convx_2!(II, O, B, u, Du, Dv, cap, ny)
#         end
#     end
#     @inbounds @threads for II in b_bottom[2:end-1]
#         vec_convx_3!(II, O, v, cap, ny)
#     end
#     @inbounds @threads for II in b_top[2:end-1]
#         vec_convx_4!(II, O, v, cap, ny)
#     end
#     @inbounds @threads for II in b_left[2:end]
#         vec_convx_5!(II, O, v, cap, n, ny, BC)
#     end
#     @inbounds @threads for II in b_left[1:end-1]
#         vec_convx_6!(II, O, v, cap, n, ny, BC)
#     end
#     @inbounds @threads for II in b_right[2:end]
#         vec_convx_7!(II, O, v, cap, n, ny, BC)
#     end
#     @inbounds @threads for II in b_right[1:end-1]
#         vec_convx_8!(II, O, v, cap, n, ny, BC)
#     end

#     if is_periodic(BC.left) && is_periodic(BC.right)
#         @inbounds for (II, JJ) in zip(b_left, b_right)
#             pII = lexicographic(II, ny)
#             pJJ = lexicographic(JJ, ny)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(JJ))
            
#             Auim1, Aui = A1_1 * u[JJ], A3_1 * u[II]

#             Au1 = 0.5 * (Auim1 + Aui)

#             @inbounds O[pII,pII] += -0.5 * Au1
#             @inbounds O[pII,pJJ] = -0.5 * Au1

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[II]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[JJ]

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_x[δy⁺(δx⁻(JJ))]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_x[δx⁻(JJ)]

#             @inbounds B[pII] += -0.25 * Du_x[II] * (A3_1 - B1_1) * Du_x[II]
#             @inbounds B[pII] += -0.25 * Du_x[II] * (B1_1 - A1_1) * Du_x[JJ]

#             @inbounds B[pII] += -0.25 * Du_x[II] * (A4_1 - B2_1) * Dv_x[δy⁺(δx⁻(JJ))]
#             @inbounds B[pII] += -0.25 * Du_x[II] * (B2_1 - A2_1) * Dv_x[δx⁻(JJ)]
#         end
#         @inbounds for (II, JJ) in zip(b_right, b_left)
#             pII = lexicographic(II, ny)
#             pJJ = lexicographic(JJ, ny)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
#             Aui, Auip1 = A1_2 * u[II], A3_2 * u[JJ]

#             Au3 = 0.5 * (Aui + Auip1)

#             @inbounds O[pII,pII] += 0.5 * Au3
#             @inbounds O[pII,pJJ] = 0.5 * Au3

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[JJ]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[II]

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_x[δy⁺(δx⁻(II))]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_x[δx⁻(II)]

#             @inbounds B[pII] += -0.25 * Du_x[II] * (A3_2 - B1_2) * Du_x[JJ]
#             @inbounds B[pII] += -0.25 * Du_x[II] * (B1_2 - A1_2) * Du_x[II]

#             @inbounds B[pII] += -0.25 * Du_x[II] * (A4_2 - B2_2) * Dv_x[δy⁺(δx⁻(II))]
#             @inbounds B[pII] += -0.25 * Du_x[II] * (B2_2 - A2_2) * Dv_x[δx⁻(II)]
#         end
#     end
#     if is_periodic(BC.bottom) && is_periodic(BC.top)
#         @inbounds for (II,JJ) in zip(b_bottom[2:end-1], b_top[2:end-1])
#             pII = lexicographic(II, ny)
#             pJJ = lexicographic(JJ, ny)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Avim1jm1 = A2_1 * v[δx⁻(II)]
#             Avip1jm1 = A2_2 * v[II]
    
#             Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    
#             @inbounds O[pII,pII] += -0.5 * Au2
#             @inbounds O[pII,pJJ] = -0.5 * Au2
#         end
#         @inbounds for (II,JJ) in zip(b_top[2:end-1], b_bottom[2:end-1])
#             pII = lexicographic(II, ny)
#             pJJ = lexicographic(JJ, ny)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Avim1jp1 = A4_1 * v[δx⁻(δy⁺(JJ))]
#             Avip1jp1 = A4_2 * v[δy⁺(JJ)]

#             Au4 = 0.5 * (Avim1jp1 + Avip1jp1)
    
#             @inbounds O[pII,pII] += 0.5 * Au4
#             @inbounds O[pII,pJJ] = 0.5 * Au4
#         end

#         ii = b_left[1]
#         pii = lexicographic(ii, ny)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Avip1jm1 = A2_2 * v[ii]

#         Au2 = 0.5 * Avip1jm1

#         if is_periodic(BC.left)
#             JJ = ii + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jm1 = A2_1 * v[δx⁻(JJ)]
#             Au2 += 0.5 * Avim1jm1
#         end

#         JJ = ii + CartesianIndex(ny-1, 0)
#         pJJ = lexicographic(JJ, ny)
#         @inbounds O[pii,pii] += -0.5 * Au2
#         @inbounds O[pii,pJJ] = -0.5 * Au2
        
#         ii = b_left[end]
#         pii = lexicographic(ii, ny)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Avip1jp1 = A4_2 * v[δy⁺(ii)]

#         Au4 = 0.5 * Avip1jp1

#         if is_periodic(BC.left)
#             JJ = ii + CartesianIndex(0, n)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
#             Au4 += 0.5 * Avim1jp1
#         end

#         JJ = ii + CartesianIndex(-ny+1, 0)
#         pJJ = lexicographic(JJ, ny)
#         @inbounds O[pii,pii] += 0.5 * Au4
#         @inbounds O[pii,pJJ] = 0.5 * Au4

#         ii = b_right[1]
#         pii = lexicographic(ii, ny)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
#         Avim1jm1 = A2_1 * v[δx⁻(ii)]

#         Au2 = 0.5 * Avim1jm1

#         if is_periodic(BC.right)
#             JJ = ii + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jm1 = A2_2 * v[JJ]
#             Au2 += 0.5 * Avip1jm1
#         end

#         JJ = ii + CartesianIndex(ny-1, 0)
#         pJJ = lexicographic(JJ, ny)
#         @inbounds O[pii,pii] += -0.5 * Au2
#         @inbounds O[pii,pJJ] = -0.5 * Au2
        
#         ii = b_right[end]
#         pii = lexicographic(ii, ny)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
#         Avim1jp1 = A4_1 * v[δy⁺(δx⁻(ii))]

#         Au4 = 0.5 * Avim1jp1

#         if is_periodic(BC.right)
#             JJ = ii + CartesianIndex(0, -n)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Avip1jp1 = A4_2 * v[δy⁺(JJ)]
#             Au4 += 0.5 * Avip1jp1
#         end
        
#         JJ = ii + CartesianIndex(-ny+1, 0)
#         pJJ = lexicographic(JJ, ny)
#         @inbounds O[pii,pii] += 0.5 * Au4
#         @inbounds O[pii,pJJ] = 0.5 * Au4
#     end

#     return nothing
# end

# function fill_inside_conv!(::Type{GridFCy}, O, B, u, v, Du, Dv_x, Dv_y, cap, ny, II)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avim1, Avi, Avip1 = A2_1 * v[δy⁻(II)], A2_2 * v[II], A4_2 * v[δy⁺(II)]
#     Auim1jm1, Auip1jm1 = A1_1 * u[δy⁻(II)], A3_1 * u[δx⁺(δy⁻(II))]
#     Auim1jp1, Auip1jp1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

#     Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
#     Au2 = 0.5 * (Avim1 + Avi)
#     Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
#     Au4 = 0.5 * (Avi + Avip1)

#     @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
#     @inbounds O[pII,pII+ny+1] = 0.5 * Au3
#     @inbounds O[pII,pII-ny-1] = -0.5 * Au1
#     @inbounds O[pII,pII+1] = 0.5 * Au4
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - B2_1) * Dv_y[II]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δy⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]
#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

#     @inbounds B[pII] += -0.25 * Dv_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv_y[II] * (B2_2 - B2_1) * Dv_y[II]
#     @inbounds B[pII] += -0.25 * Dv_y[II] * (B2_1 - A2_1) * Dv_y[δy⁻(II)]

#     @inbounds B[pII] += -0.25 * Dv_x[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#     @inbounds B[pII] += -0.25 * Dv_x[II] * (B1_1 - A1_1) * Du[δy⁻(II)]
#     @inbounds B[pII] += -0.25 * Dv_x[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv_x[II] * (B1_2 - A1_2) * Du[II]
# end

# function vec_convy_1!(II, O, B, v, Du, Dv, cap, ny)
#     pII = lexicographic(II, ny+1)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Avi, Avip1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

#     Au4 = 0.5 * (Avi + Avip1)

#     @inbounds O[pII,pII] += 0.5 * Au4
#     @inbounds O[pII,pII+1] = 0.5 * Au4

#     @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

#     @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

#     @inbounds B[pII] += -0.25 * Dv[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv[II] * (B2_2 - A2_2) * Dv[II]

#     @inbounds B[pII] += -0.25 * Dv[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
#     @inbounds B[pII] += -0.25 * Dv[II] * (B1_2 - A1_2) * Du[II]

#     return nothing
# end

# function vec_convy_2!(II, O, B, v, Du, Dv, cap, ny)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
#     Avim1, Avi = A2_1 * v[δy⁻(II)], A4_1 * v[II]

#     Au2 = 0.5 * (Avim1 + Avi)

#     @inbounds O[pII,pII] += -0.5 * Au2
#     @inbounds O[pII,pII-1] = -0.5 * Au2

#     @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[II]
#     @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δy⁻(II)]

#     @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#     @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]

#     @inbounds B[pII] += -0.25 * Dv[II] * (A4_1 - B2_1) * Dv[II]
#     @inbounds B[pII] += -0.25 * Dv[II] * (B2_1 - A2_1) * Dv[δy⁻(II)]

#     @inbounds B[pII] += -0.25 * Dv[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
#     @inbounds B[pII] += -0.25 * Dv[II] * (B1_1 - A1_1) * Du[δy⁻(II)]

#     return nothing
# end

# function vec_convy_3!(II, O, u, cap, ny)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
#     Auip1jp1 = A3_2 * u[δx⁺(II)]

#     Au3 = 0.5 * (Auip1jm1 + Auip1jp1)

#     @inbounds O[pII,pII] += 0.5 * Au3
#     @inbounds O[pII,pII+ny+1] = 0.5 * Au3

#     return nothing
# end

# function vec_convy_4!(II, O, u, cap, ny)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Auim1jm1 = A1_1 * u[δy⁻(II)]
#     Auim1jp1 = A1_2 * u[II]

#     Au1 = 0.5 * (Auim1jm1 + Auim1jp1)

#     @inbounds O[pII,pII] += -0.5 * Au1
#     @inbounds O[pII,pII-ny-1] = -0.5 * Au1

#     return nothing
# end

# function vec_convy_5!(II, O, u, cap, ny, BC)
#     pII = lexicographic(II, ny+1)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Auim1jp1 = A1_2 * u[II]

#     Au1 = 0.5 * Auim1jp1

#     if is_periodic(BC.bottom)
#         JJ = II + CartesianIndex(ny-1, 0)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#         Auim1jm1 = A1_1 * u[JJ]
#         Au1 += 0.5 * Auim1jm1
#     end

#     @inbounds O[pII,pII] += -0.5 * Au1
#     @inbounds O[pII,pII-ny-1] = -0.5 * Au1

#     return nothing
# end

# function vec_convy_6!(II, O, u, cap, ny, BC)
#     pII = lexicographic(II, ny+1)
#     A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
#     Auip1jp1 = A3_2 * u[δx⁺(II)]

#     Au3 = 0.5 * Auip1jp1

#     if is_periodic(BC.bottom)
#         JJ = II + CartesianIndex(ny-1, 0)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#         Auip1jm1 = A3_1 * u[δx⁺(JJ)]
#         Au3 += 0.5 * Auip1jm1
#     end

#     @inbounds O[pII,pII] += 0.5 * Au3
#     @inbounds O[pII,pII+ny+1] = 0.5 * Au3

#     return nothing
# end

# function vec_convy_7!(II, O, u, cap, ny, BC)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
#     Auim1jm1 = A1_1 * u[δy⁻(II)]

#     Au1 = 0.5 * Auim1jm1

#     if is_periodic(BC.top)
#         JJ = II + CartesianIndex(-ny, 0)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#         Auim1jp1 = A1_2 * u[JJ]
#         Au1 += 0.5 * Auim1jp1
#     end

#     @inbounds O[pII,pII] += -0.5 * Au1
#     @inbounds O[pII,pII-ny-1] = -0.5 * Au1

#     return nothing
# end

# function vec_convy_8!(II, O, u, cap, ny, BC)
#     pII = lexicographic(II, ny+1)
#     A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
#     Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]

#     Au3 = 0.5 * Auip1jm1

#     if is_periodic(BC.top)
#         JJ = II + CartesianIndex(-ny, 0)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#         Auip1jp1 = A3_2 * u[δx⁺(JJ)]
#         Au3 += 0.5 * Auip1jp1
#     end

#     @inbounds O[pII,pII] += 0.5 * Au3
#     @inbounds O[pII,pII+ny+1] = 0.5 * Au3

#     return nothing
# end

# function vector_convection!(::Dirichlet, ::Type{GridFCy}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
#     B .= 0.0
#     @inbounds @threads for II in inside
#         fill_inside_conv!(GridFCy, O, B, u, v, Du_x, Dv_x, Dv_y, cap, ny, II)
#     end

#     @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
#         pII = lexicographic(II, ny+1)
#         @inbounds O[pII,pII] = 0.0
#     end
#     bnds = (b_left[2:end-1], b_bottom, b_right[2:end-1])
#     bc = ((Du_x, Dv_x), (Du_y, Dv_y), (Du_x, Dv_x))
#     for (bnd, (Du, Dv)) in zip(bnds, bc)
#         @inbounds @threads for II in bnd
#             vec_convy_1!(II, O, B, v, Du, Dv, cap, ny)
#         end
#     end
#     bnds = (b_left[2:end-1], b_right[2:end-1], b_top)
#     bc = ((Du_x, Dv_x), (Du_x, Dv_x), (Du_y, Dv_y))
#     for (bnd, (Du, Dv)) in zip(bnds, bc)
#         @inbounds @threads for II in bnd
#             vec_convy_2!(II, O, B, v, Du, Dv, cap, ny)
#         end
#     end
#     @inbounds @threads for II in b_left[2:end-1]
#         vec_convy_3!(II, O, u, cap, ny)
#     end
#     @inbounds @threads for II in b_right[2:end-1]
#         vec_convy_4!(II, O, u, cap, ny)
#     end
#     @inbounds @threads for II in b_bottom[2:end]
#         vec_convy_5!(II, O, u, cap, ny, BC)
#     end
#     @inbounds @threads for II in b_bottom[1:end-1]
#         vec_convy_6!(II, O, u, cap, ny, BC)
#     end
#     @inbounds @threads for II in b_top[2:end]
#         vec_convy_7!(II, O, u, cap, ny, BC)
#     end
#     @inbounds @threads for II in b_top[1:end-1]
#         vec_convy_8!(II, O, u, cap, ny, BC)
#     end

#     if is_periodic(BC.bottom) && is_periodic(BC.top)
#         @inbounds for (II, JJ) in zip(b_bottom, b_top)
#             pII = lexicographic(II, ny+1)
#             pJJ = lexicographic(JJ, ny+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(JJ))
            
#             Avim1, Avi = A2_1 * v[JJ], A4_1 * v[II]

#             Au2 = 0.5 * (Avim1 + Avi)

#             @inbounds O[pII,pII] += -0.5 * Au2
#             @inbounds O[pII,pJJ] = -0.5 * Au2

#             @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[II]
#             @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[JJ]

#             @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_y[δx⁺(δy⁻(JJ))]
#             @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_y[δy⁻(JJ)]

#             @inbounds B[pII] += -0.25 * Dv_y[II] * (A4_1 - B2_1) * Dv_y[II]
#             @inbounds B[pII] += -0.25 * Dv_y[II] * (B2_1 - A2_1) * Dv_y[JJ]

#             @inbounds B[pII] += -0.25 * Dv_y[II] * (A3_1 - B1_1) * Du_y[δx⁺(δy⁻(JJ))]
#             @inbounds B[pII] += -0.25 * Dv_y[II] * (B1_1 - A1_1) * Du_y[δy⁻(JJ)]
#         end
#         @inbounds for (II, JJ) in zip(b_top, b_bottom)
#             pII = lexicographic(II, ny+1)
#             pJJ = lexicographic(JJ, ny+1)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
#             Avi, Avip1 = A2_2 * v[II], A4_2 * v[JJ]

#             Au4 = 0.5 * (Avi + Avip1)

#             @inbounds O[pII,pII] += 0.5 * Au4
#             @inbounds O[pII,pJJ] = 0.5 * Au4

#             @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[JJ]
#             @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

#             @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_y[δx⁺(δy⁻(II))]
#             @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_y[δy⁻(II)]

#             @inbounds B[pII] += -0.25 * Dv_y[II] * (A4_2 - B2_2) * Dv_y[JJ]
#             @inbounds B[pII] += -0.25 * Dv_y[II] * (B2_2 - A2_2) * Dv_y[II]

#             @inbounds B[pII] += -0.25 * Dv_y[II] * (A3_2 - B1_2) * Du_y[δx⁺(δy⁻(II))]
#             @inbounds B[pII] += -0.25 * Dv_y[II] * (B1_2 - A1_2) * Du_y[δy⁻(II)]
#         end
#     end
#     if is_periodic(BC.left) && is_periodic(BC.right)
#         @inbounds for (II,JJ) in zip(b_left[2:end-1], b_right[2:end-1])
#             pII = lexicographic(II, ny+1)
#             pJJ = lexicographic(JJ, ny+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Auim1jm1 = A1_1 * u[δy⁻(II)]
#             Auim1jp1 = A1_2 * u[II]
    
#             Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    
#             @inbounds O[pII,pII] += -0.5 * Au1
#             @inbounds O[pII,pJJ] = -0.5 * Au1
#         end
#         @inbounds for (II,JJ) in zip(b_right[2:end-1], b_left[2:end-1])
#             pII = lexicographic(II, ny+1)
#             pJJ = lexicographic(JJ, ny+1)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
#             Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
#             Auip1jp1 = A3_2 * u[δx⁺(II)]
    
#             Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    
#             @inbounds O[pII,pII] += 0.5 * Au3
#             @inbounds O[pII,pJJ] = 0.5 * Au3
#         end

#         ii = b_bottom[1]
#         pii = lexicographic(ii, ny+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Auim1jp1 = A1_2 * u[ii]

#         Au1 = 0.5 * Auim1jp1

#         if is_periodic(BC.bottom)
#             JJ = ii + CartesianIndex(ny-1, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auim1jm1 = A1_1 * u[JJ]
#             Au1 += 0.5 * Auim1jm1
#         end

#         JJ = ii + CartesianIndex(0, n-1)
#         pJJ = lexicographic(JJ, ny+1)
#         @inbounds O[pii,pii] += -0.5 * Au1
#         @inbounds O[pii,pJJ] = -0.5 * Au1
        
#         ii = b_bottom[end]
#         pii = lexicographic(ii, ny+1)
#         A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
#         Auip1jp1 = A3_2 * u[δx⁺(ii)]

#         Au3 = 0.5 * Auip1jp1

#         if is_periodic(BC.bottom)
#             JJ = ii + CartesianIndex(ny-1, 0)
#             A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
#             Auip1jm1 = A3_1 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jm1
#         end

#         JJ = ii + CartesianIndex(0, -n+1)
#         pJJ = lexicographic(JJ, ny+1)
#         @inbounds O[pii,pii] += 0.5 * Au3
#         @inbounds O[pii,pJJ] = 0.5 * Au3
        
#         ii = b_top[1]
#         pii = lexicographic(ii, ny+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
#         Auim1jm1 = A1_1 * u[δy⁻(ii)]

#         Au1 = 0.5 * Auim1jm1

#         if is_periodic(BC.top)
#             JJ = ii + CartesianIndex(-ny, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auim1jp1 = A1_2 * u[JJ]
#             Au1 += 0.5 * Auim1jp1
#         end

#         JJ = ii + CartesianIndex(0, n-1)
#         pJJ = lexicographic(JJ, ny+1)
#         @inbounds O[pii,pii] += -0.5 * Au1
#         @inbounds O[pii,pJJ] = -0.5 * Au1
        
#         ii = b_top[end]
#         pii = lexicographic(ii, ny+1)
#         A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
#         Auip1jm1 = A3_1 * u[δx⁺(δy⁻(ii))]

#         Au3 = 0.5 * Auip1jm1

#         if is_periodic(BC.top)
#             JJ = ii + CartesianIndex(-ny, 0)
#             A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
#             Auip1jp1 = A3_2 * u[δx⁺(JJ)]
#             Au3 += 0.5 * Auip1jp1
#         end

#         JJ = ii + CartesianIndex(0, -n+1)
#         pJJ = lexicographic(JJ, ny+1)
#         @inbounds O[pii,pii] += 0.5 * Au3
#         @inbounds O[pii,pJJ] = 0.5 * Au3
#     end

#     return nothing
# end

function fill_inside_conv!(::Type{GridFCx}, O, B, u, v, Du_x, Du_y, Dv_y, cap, ny, II)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)

    Auim1, Aui, Auip1 = A1_1 * u[δx⁻(II)], A1_2 * u[II], A3_2 * u[δx⁺(II)]
    Avim1jm1, Avim1jp1 = A2_1 * v[δx⁻(II)], A4_1 * v[δy⁺(δx⁻(II))]
    Avip1jm1, Avip1jp1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

    Au1 = 0.5 * (Auim1 + Aui)
    Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    Au3 = 0.5 * (Aui + Auip1)
    Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

    @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
    @inbounds O[pII,pII+ny] = 0.5 * Au3
    @inbounds O[pII,pII-ny] = -0.5 * Au1
    @inbounds O[pII,pII+1] = 0.5 * Au4
    @inbounds O[pII,pII-1] = -0.5 * Au2

    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - B1_1) * Du_x[II]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[δx⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

    @inbounds B[pII] += -0.25 * u[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
    @inbounds B[pII] += -0.25 * u[II] * (B1_2 - B1_1) * Du_x[II]
    @inbounds B[pII] += -0.25 * u[II] * (B1_1 - A1_1) * Du_x[δx⁻(II)]

    @inbounds B[pII] += -0.25 * u[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
    @inbounds B[pII] += -0.25 * u[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
    @inbounds B[pII] += -0.25 * u[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds B[pII] += -0.25 * u[II] * (B2_2 - A2_2) * Dv_y[II]

    return nothing
end

function vec_convx_1!(II, O, B, u, Du, Dv, cap, ny)
    pII = lexicographic(II, ny)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Aui, Auip1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

    Au3 = 0.5 * (Aui + Auip1)

    @inbounds O[pII,pII] += 0.5 * Au3
    @inbounds O[pII,pII+ny] = 0.5 * Au3

    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

    @inbounds B[pII] += -0.25 * u[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds B[pII] += -0.25 * u[II] * (B1_2 - A1_2) * Du[II]

    @inbounds B[pII] += -0.25 * u[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
    @inbounds B[pII] += -0.25 * u[II] * (B2_2 - A2_2) * Dv[II]

    return nothing
end

function vec_convx_2!(II, O, B, u, Du, Dv, cap, ny)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
    Auim1, Aui = A1_1 * u[δx⁻(II)], A3_1 * u[II]

    Au1 = 0.5 * (Auim1 + Aui)

    @inbounds O[pII,pII] += -0.5 * Au1
    @inbounds O[pII,pII-ny] = -0.5 * Au1

    @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[II]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δx⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δx⁻(II)]

    @inbounds B[pII] += -0.25 * u[II] * (A3_1 - B1_1) * Du[II]
    @inbounds B[pII] += -0.25 * u[II] * (B1_1 - A1_1) * Du[δx⁻(II)]

    @inbounds B[pII] += -0.25 * u[II] * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
    @inbounds B[pII] += -0.25 * u[II] * (B2_1 - A2_1) * Dv[δx⁻(II)]

    return nothing
end

function vec_convx_3!(II, O, v, cap, ny)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]
    Avip1jp1 = A4_2 * v[δy⁺(II)]

    Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

    @inbounds O[pII,pII] += 0.5 * Au4
    @inbounds O[pII,pII+1] = 0.5 * Au4

    return nothing
end

function vec_convx_4!(II, O, v, cap, ny)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avim1jm1 = A2_1 * v[δx⁻(II)]
    Avip1jm1 = A2_2 * v[II]

    Au2 = 0.5 * (Avim1jm1 + Avip1jm1)

    @inbounds O[pII,pII] += -0.5 * Au2
    @inbounds O[pII,pII-1] = -0.5 * Au2

    return nothing
end

function vec_convx_5!(II, O, v, cap, n, ny, BC)
    pII = lexicographic(II, ny)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avip1jm1 = A2_2 * v[II]

    Au2 = 0.5 * Avip1jm1

    if is_periodic(BC.left)
        JJ = II + CartesianIndex(0, n-1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
        Avim1jm1 = A2_1 * v[JJ]
        Au2 += 0.5 * Avim1jm1
    end

    @inbounds O[pII,pII] += -0.5 * Au2
    @inbounds O[pII,pII-1] = -0.5 * Au2

    return nothing
end

function vec_convx_6!(II, O, v, cap, n, ny, BC)
    pII = lexicographic(II, ny)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avip1jp1 = A4_2 * v[δy⁺(II)]

    Au4 = 0.5 * Avip1jp1

    if is_periodic(BC.left)
        JJ = II + CartesianIndex(0, n-1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
        Avim1jp1 = A4_1 * v[δy⁺(JJ)]
        Au4 += 0.5 * Avim1jp1
    end

    @inbounds O[pII,pII] += 0.5 * Au4
    @inbounds O[pII,pII+1] = 0.5 * Au4

    return nothing
end

function vec_convx_7!(II, O, v, cap, n, ny, BC)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
    Avim1jm1 = A2_1 * v[δx⁻(II)]

    Au2 = 0.5 * Avim1jm1

    if is_periodic(BC.right)
        JJ = II + CartesianIndex(0, -n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
        Avip1jm1 = A2_2 * v[JJ]
        Au2 += 0.5 * Avip1jm1
    end

    @inbounds O[pII,pII] += -0.5 * Au2
    @inbounds O[pII,pII-1] = -0.5 * Au2

    return nothing
end

function vec_convx_8!(II, O, v, cap, n, ny, BC)
    pII = lexicographic(II, ny)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
    
    Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]

    Au4 = 0.5 * Avim1jp1

    if is_periodic(BC.right)
        JJ = II + CartesianIndex(0, -n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
        Avip1jp1 = A4_2 * v[δy⁺(JJ)]
        Au4 += 0.5 * Avip1jp1
    end

    @inbounds O[pII,pII] += 0.5 * Au4
    @inbounds O[pII,pII+1] = 0.5 * Au4

    return nothing
end

function vector_convection!(::Dirichlet, ::Type{GridFCx}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(GridFCx, O, B, u, v, Du_x, Du_y, Dv_y, cap, ny, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, ny)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left, b_bottom[2:end-1], b_top[2:end-1])
    bc = ((Du_x, Dv_x), (Du_y, Dv_y), (Du_y, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convx_1!(II, O, B, u, Du, Dv, cap, ny)
        end
    end
    bnds = (b_bottom[2:end-1], b_right, b_top[2:end-1])
    bc = ((Du_y, Dv_y), (Du_x, Dv_x), (Du_y, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convx_2!(II, O, B, u, Du, Dv, cap, ny)
        end
    end
    @inbounds @threads for II in b_bottom[2:end-1]
        vec_convx_3!(II, O, v, cap, ny)
    end
    @inbounds @threads for II in b_top[2:end-1]
        vec_convx_4!(II, O, v, cap, ny)
    end
    @inbounds @threads for II in b_left[2:end]
        vec_convx_5!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_left[1:end-1]
        vec_convx_6!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_right[2:end]
        vec_convx_7!(II, O, v, cap, n, ny, BC)
    end
    @inbounds @threads for II in b_right[1:end-1]
        vec_convx_8!(II, O, v, cap, n, ny, BC)
    end

    if is_periodic(BC.left) && is_periodic(BC.right)
        @inbounds for (II, JJ) in zip(b_left, b_right)
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(JJ))
            
            Auim1, Aui = A1_1 * u[JJ], A3_1 * u[II]

            Au1 = 0.5 * (Auim1 + Aui)

            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[II]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[JJ]

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_x[δy⁺(δx⁻(JJ))]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_x[δx⁻(JJ)]

            @inbounds B[pII] += -0.25 * u[II] * (A3_1 - B1_1) * Du_x[II]
            @inbounds B[pII] += -0.25 * u[II] * (B1_1 - A1_1) * Du_x[JJ]

            @inbounds B[pII] += -0.25 * u[II] * (A4_1 - B2_1) * Dv_x[δy⁺(δx⁻(JJ))]
            @inbounds B[pII] += -0.25 * u[II] * (B2_1 - A2_1) * Dv_x[δx⁻(JJ)]
        end
        @inbounds for (II, JJ) in zip(b_right, b_left)
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
            Aui, Auip1 = A1_2 * u[II], A3_2 * u[JJ]

            Au3 = 0.5 * (Aui + Auip1)

            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[JJ]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[II]

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_x[δy⁺(δx⁻(II))]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_x[δx⁻(II)]

            @inbounds B[pII] += -0.25 * u[II] * (A3_2 - B1_2) * Du_x[JJ]
            @inbounds B[pII] += -0.25 * u[II] * (B1_2 - A1_2) * Du_x[II]

            @inbounds B[pII] += -0.25 * u[II] * (A4_2 - B2_2) * Dv_x[δy⁺(δx⁻(II))]
            @inbounds B[pII] += -0.25 * u[II] * (B2_2 - A2_2) * Dv_x[δx⁻(II)]
        end
    end
    if is_periodic(BC.bottom) && is_periodic(BC.top)
        @inbounds for (II,JJ) in zip(b_bottom[2:end-1], b_top[2:end-1])
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Avim1jm1 = A2_1 * v[δx⁻(II)]
            Avip1jm1 = A2_2 * v[II]
    
            Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    
            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2
        end
        @inbounds for (II,JJ) in zip(b_top[2:end-1], b_bottom[2:end-1])
            pII = lexicographic(II, ny)
            pJJ = lexicographic(JJ, ny)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Avim1jp1 = A4_1 * v[δx⁻(δy⁺(JJ))]
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]

            Au4 = 0.5 * (Avim1jp1 + Avip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4
        end

        ii = b_left[1]
        pii = lexicographic(ii, ny)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Avip1jm1 = A2_2 * v[ii]

        Au2 = 0.5 * Avip1jm1

        if is_periodic(BC.left)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Avim1jm1 = A2_1 * v[δx⁻(JJ)]
            Au2 += 0.5 * Avim1jm1
        end

        JJ = ii + CartesianIndex(ny-1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_left[end]
        pii = lexicographic(ii, ny)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Avip1jp1 = A4_2 * v[δy⁺(ii)]

        Au4 = 0.5 * Avip1jp1

        if is_periodic(BC.left)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
            Au4 += 0.5 * Avim1jp1
        end

        JJ = ii + CartesianIndex(-ny+1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4

        ii = b_right[1]
        pii = lexicographic(ii, ny)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
        Avim1jm1 = A2_1 * v[δx⁻(ii)]

        Au2 = 0.5 * Avim1jm1

        if is_periodic(BC.right)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Avip1jm1 = A2_2 * v[JJ]
            Au2 += 0.5 * Avip1jm1
        end

        JJ = ii + CartesianIndex(ny-1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_right[end]
        pii = lexicographic(ii, ny)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii))
        
        Avim1jp1 = A4_1 * v[δy⁺(δx⁻(ii))]

        Au4 = 0.5 * Avim1jp1

        if is_periodic(BC.right)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]
            Au4 += 0.5 * Avip1jp1
        end
        
        JJ = ii + CartesianIndex(-ny+1, 0)
        pJJ = lexicographic(JJ, ny)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4
    end

    return nothing
end

function fill_inside_conv!(::Type{GridFCy}, O, B, u, v, Du, Dv_x, Dv_y, cap, ny, II)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avim1, Avi, Avip1 = A2_1 * v[δy⁻(II)], A2_2 * v[II], A4_2 * v[δy⁺(II)]
    Auim1jm1, Auip1jm1 = A1_1 * u[δy⁻(II)], A3_1 * u[δx⁺(δy⁻(II))]
    Auim1jp1, Auip1jp1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

    Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    Au2 = 0.5 * (Avim1 + Avi)
    Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    Au4 = 0.5 * (Avi + Avip1)

    @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
    @inbounds O[pII,pII+ny+1] = 0.5 * Au3
    @inbounds O[pII,pII-ny-1] = -0.5 * Au1
    @inbounds O[pII,pII+1] = 0.5 * Au4
    @inbounds O[pII,pII-1] = -0.5 * Au2

    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - B2_1) * Dv_y[II]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δy⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]
    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

    @inbounds B[pII] += -0.25 * v[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds B[pII] += -0.25 * v[II] * (B2_2 - B2_1) * Dv_y[II]
    @inbounds B[pII] += -0.25 * v[II] * (B2_1 - A2_1) * Dv_y[δy⁻(II)]

    @inbounds B[pII] += -0.25 * v[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
    @inbounds B[pII] += -0.25 * v[II] * (B1_1 - A1_1) * Du[δy⁻(II)]
    @inbounds B[pII] += -0.25 * v[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds B[pII] += -0.25 * v[II] * (B1_2 - A1_2) * Du[II]
end

function vec_convy_1!(II, O, B, v, Du, Dv, cap, ny)
    pII = lexicographic(II, ny+1)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Avi, Avip1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

    Au4 = 0.5 * (Avi + Avip1)

    @inbounds O[pII,pII] += 0.5 * Au4
    @inbounds O[pII,pII+1] = 0.5 * Au4

    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

    @inbounds B[pII] += -0.25 * v[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
    @inbounds B[pII] += -0.25 * v[II] * (B2_2 - A2_2) * Dv[II]

    @inbounds B[pII] += -0.25 * v[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
    @inbounds B[pII] += -0.25 * v[II] * (B1_2 - A1_2) * Du[II]

    return nothing
end

function vec_convy_2!(II, O, B, v, Du, Dv, cap, ny)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
    Avim1, Avi = A2_1 * v[δy⁻(II)], A4_1 * v[II]

    Au2 = 0.5 * (Avim1 + Avi)

    @inbounds O[pII,pII] += -0.5 * Au2
    @inbounds O[pII,pII-1] = -0.5 * Au2

    @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[II]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δy⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]

    @inbounds B[pII] += -0.25 * v[II] * (A4_1 - B2_1) * Dv[II]
    @inbounds B[pII] += -0.25 * v[II] * (B2_1 - A2_1) * Dv[δy⁻(II)]

    @inbounds B[pII] += -0.25 * v[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
    @inbounds B[pII] += -0.25 * v[II] * (B1_1 - A1_1) * Du[δy⁻(II)]

    return nothing
end

function vec_convy_3!(II, O, u, cap, ny)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
    Auip1jp1 = A3_2 * u[δx⁺(II)]

    Au3 = 0.5 * (Auip1jm1 + Auip1jp1)

    @inbounds O[pII,pII] += 0.5 * Au3
    @inbounds O[pII,pII+ny+1] = 0.5 * Au3

    return nothing
end

function vec_convy_4!(II, O, u, cap, ny)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Auim1jm1 = A1_1 * u[δy⁻(II)]
    Auim1jp1 = A1_2 * u[II]

    Au1 = 0.5 * (Auim1jm1 + Auim1jp1)

    @inbounds O[pII,pII] += -0.5 * Au1
    @inbounds O[pII,pII-ny-1] = -0.5 * Au1

    return nothing
end

function vec_convy_5!(II, O, u, cap, ny, BC)
    pII = lexicographic(II, ny+1)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Auim1jp1 = A1_2 * u[II]

    Au1 = 0.5 * Auim1jp1

    if is_periodic(BC.bottom)
        JJ = II + CartesianIndex(ny-1, 0)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
        Auim1jm1 = A1_1 * u[JJ]
        Au1 += 0.5 * Auim1jm1
    end

    @inbounds O[pII,pII] += -0.5 * Au1
    @inbounds O[pII,pII-ny-1] = -0.5 * Au1

    return nothing
end

function vec_convy_6!(II, O, u, cap, ny, BC)
    pII = lexicographic(II, ny+1)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
    
    Auip1jp1 = A3_2 * u[δx⁺(II)]

    Au3 = 0.5 * Auip1jp1

    if is_periodic(BC.bottom)
        JJ = II + CartesianIndex(ny-1, 0)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
        Auip1jm1 = A3_1 * u[δx⁺(JJ)]
        Au3 += 0.5 * Auip1jm1
    end

    @inbounds O[pII,pII] += 0.5 * Au3
    @inbounds O[pII,pII+ny+1] = 0.5 * Au3

    return nothing
end

function vec_convy_7!(II, O, u, cap, ny, BC)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
    Auim1jm1 = A1_1 * u[δy⁻(II)]

    Au1 = 0.5 * Auim1jm1

    if is_periodic(BC.top)
        JJ = II + CartesianIndex(-ny, 0)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
        Auim1jp1 = A1_2 * u[JJ]
        Au1 += 0.5 * Auim1jp1
    end

    @inbounds O[pII,pII] += -0.5 * Au1
    @inbounds O[pII,pII-ny-1] = -0.5 * Au1

    return nothing
end

function vec_convy_8!(II, O, u, cap, ny, BC)
    pII = lexicographic(II, ny+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
    
    Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]

    Au3 = 0.5 * Auip1jm1

    if is_periodic(BC.top)
        JJ = II + CartesianIndex(-ny, 0)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
        Auip1jp1 = A3_2 * u[δx⁺(JJ)]
        Au3 += 0.5 * Auip1jp1
    end

    @inbounds O[pII,pII] += 0.5 * Au3
    @inbounds O[pII,pII+ny+1] = 0.5 * Au3

    return nothing
end

function vector_convection!(::Dirichlet, ::Type{GridFCy}, O, B, u, v, Du_x, Du_y, Dv_x, Dv_y, cap, n, ny, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(GridFCy, O, B, u, v, Du_x, Dv_x, Dv_y, cap, ny, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, ny+1)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left[2:end-1], b_bottom, b_right[2:end-1])
    bc = ((Du_x, Dv_x), (Du_y, Dv_y), (Du_x, Dv_x))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convy_1!(II, O, B, v, Du, Dv, cap, ny)
        end
    end
    bnds = (b_left[2:end-1], b_right[2:end-1], b_top)
    bc = ((Du_x, Dv_x), (Du_x, Dv_x), (Du_y, Dv_y))
    for (bnd, (Du, Dv)) in zip(bnds, bc)
        @inbounds @threads for II in bnd
            vec_convy_2!(II, O, B, v, Du, Dv, cap, ny)
        end
    end
    @inbounds @threads for II in b_left[2:end-1]
        vec_convy_3!(II, O, u, cap, ny)
    end
    @inbounds @threads for II in b_right[2:end-1]
        vec_convy_4!(II, O, u, cap, ny)
    end
    @inbounds @threads for II in b_bottom[2:end]
        vec_convy_5!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_bottom[1:end-1]
        vec_convy_6!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_top[2:end]
        vec_convy_7!(II, O, u, cap, ny, BC)
    end
    @inbounds @threads for II in b_top[1:end-1]
        vec_convy_8!(II, O, u, cap, ny, BC)
    end

    if is_periodic(BC.bottom) && is_periodic(BC.top)
        @inbounds for (II, JJ) in zip(b_bottom, b_top)
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(JJ))
            
            Avim1, Avi = A2_1 * v[JJ], A4_1 * v[II]

            Au2 = 0.5 * (Avim1 + Avi)

            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[II]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[JJ]

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_y[δx⁺(δy⁻(JJ))]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_y[δy⁻(JJ)]

            @inbounds B[pII] += -0.25 * v[II] * (A4_1 - B2_1) * Dv_y[II]
            @inbounds B[pII] += -0.25 * v[II] * (B2_1 - A2_1) * Dv_y[JJ]

            @inbounds B[pII] += -0.25 * v[II] * (A3_1 - B1_1) * Du_y[δx⁺(δy⁻(JJ))]
            @inbounds B[pII] += -0.25 * v[II] * (B1_1 - A1_1) * Du_y[δy⁻(JJ)]
        end
        @inbounds for (II, JJ) in zip(b_top, b_bottom)
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            
            Avi, Avip1 = A2_2 * v[II], A4_2 * v[JJ]

            Au4 = 0.5 * (Avi + Avip1)

            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[JJ]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_y[δx⁺(δy⁻(II))]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_y[δy⁻(II)]

            @inbounds B[pII] += -0.25 * v[II] * (A4_2 - B2_2) * Dv_y[JJ]
            @inbounds B[pII] += -0.25 * v[II] * (B2_2 - A2_2) * Dv_y[II]

            @inbounds B[pII] += -0.25 * v[II] * (A3_2 - B1_2) * Du_y[δx⁺(δy⁻(II))]
            @inbounds B[pII] += -0.25 * v[II] * (B1_2 - A1_2) * Du_y[δy⁻(II)]
        end
    end
    if is_periodic(BC.left) && is_periodic(BC.right)
        @inbounds for (II,JJ) in zip(b_left[2:end-1], b_right[2:end-1])
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Auim1jm1 = A1_1 * u[δy⁻(II)]
            Auim1jp1 = A1_2 * u[II]
    
            Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    
            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1
        end
        @inbounds for (II,JJ) in zip(b_right[2:end-1], b_left[2:end-1])
            pII = lexicographic(II, ny+1)
            pJJ = lexicographic(JJ, ny+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II))
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II)
            
            Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
            Auip1jp1 = A3_2 * u[δx⁺(II)]
    
            Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3
        end

        ii = b_bottom[1]
        pii = lexicographic(ii, ny+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Auim1jp1 = A1_2 * u[ii]

        Au1 = 0.5 * Auim1jp1

        if is_periodic(BC.bottom)
            JJ = ii + CartesianIndex(ny-1, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Auim1jm1 = A1_1 * u[JJ]
            Au1 += 0.5 * Auim1jm1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_bottom[end]
        pii = lexicographic(ii, ny+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii)
        
        Auip1jp1 = A3_2 * u[δx⁺(ii)]

        Au3 = 0.5 * Auip1jp1

        if is_periodic(BC.bottom)
            JJ = ii + CartesianIndex(ny-1, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ)
            Auip1jm1 = A3_1 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jm1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
        
        ii = b_top[1]
        pii = lexicographic(ii, ny+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
        Auim1jm1 = A1_1 * u[δy⁻(ii)]

        Au1 = 0.5 * Auim1jm1

        if is_periodic(BC.top)
            JJ = ii + CartesianIndex(-ny, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Auim1jp1 = A1_2 * u[JJ]
            Au1 += 0.5 * Auim1jp1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_top[end]
        pii = lexicographic(ii, ny+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii))
        
        Auip1jm1 = A3_1 * u[δx⁺(δy⁻(ii))]

        Au3 = 0.5 * Auip1jm1

        if is_periodic(BC.top)
            JJ = ii + CartesianIndex(-ny, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ)
            Auip1jp1 = A3_2 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jp1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, ny+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
    end

    return nothing
end
