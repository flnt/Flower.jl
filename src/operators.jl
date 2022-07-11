@inline function get_capacities(cap, II, Δ)
    @inbounds ret = (cap[II,1]*Δ, cap[II,2]*Δ, cap[II,3]*Δ, cap[II,4]*Δ, cap[II,6]*Δ, cap[II,7]*Δ,
           cap[II,8]*Δ^2 + eps(0.01), cap[II,9]*Δ^2 + eps(0.01), cap[II,10]*Δ^2 + eps(0.01), cap[II,11]*Δ^2 + eps(0.01))
    return ret
end

# @inline function get_capacities_x(cap, II_1, II_2, Δ)
#     @inbounds ret = cap[II_1,3] * Δ, cap[II_1,6] * Δ, cap[II_2,1] * Δ, cap[II_2,6] * Δ
#     return ret
# end

# @inline function get_capacities_y(cap, II_1, II_2, Δ)
#     @inbounds ret = cap[II_1,4] * Δ, cap[II_1,7] * Δ, cap[II_2,2] * Δ, cap[II_2,7] * Δ
#     return ret
# end

@inline function get_capacities_x(cap, II, Δ)
    @inbounds ret = cap[II,1] * Δ, cap[II,3] * Δ, cap[II,6] * Δ
    return ret
end

@inline function get_capacities_y(cap, II, Δ)
    @inbounds ret = cap[II,2] * Δ, cap[II,4] * Δ, cap[II,7] * Δ
    return ret
end

@inline function get_capacities_convection(cap, II, Δ)
    @inbounds ret = cap[II,1]*Δ, cap[II,2]*Δ, cap[II,3]*Δ, cap[II,4]*Δ, cap[II,6]*Δ, cap[II,7]*Δ
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
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] = B1[II]*Δ / (W[II]*Δ^2+eps(0.01)) * B1[JJ]*Δ
    end
    return nothing
end

function laplacian!(::Dirichlet, L, B, Dx, Dy, cap, n, Δ, BC, inside, empty, MIXED, b_left, b_bottom, b_right, b_top)
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
    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+n] = B1 / W3 * cap[δx⁺(II),6]*Δ
        @inbounds B[pII] += -B1 / W3 * (cap[δx⁺(II),6]*Δ - A3) * Dx[δx⁺(II)]
    end
    @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-n] = B1 / W1 * cap[δx⁻(II),6]*Δ
        @inbounds B[pII] += B1 / W1 * (A1 - cap[δx⁻(II),6]*Δ) * Dx[δx⁻(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+1] = B2 / W4 * cap[δy⁺(II),7]*Δ
        @inbounds B[pII] += -B2 / W4 * (cap[δy⁺(II),7]*Δ - A4) * Dy[δy⁺(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-1] = B2 / W2 * cap[δy⁻(II),7]*Δ
        @inbounds B[pII] += B2 / W2 * (A2 - cap[δy⁻(II),7]*Δ) * Dy[δy⁻(II)]
    end

    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        if sum(abs.(L[pII,:])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
        end
    end
    @inbounds @threads for II in MIXED
        pII = lexicographic(II, n)
        if sum(abs.(L[pII,:])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
        end
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
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds L[pII,pJJ] = (A[II]*Δ)^2 / (W[II]*Δ^2+eps(0.01))
    end
    return nothing
end

function laplacian!(::Neumann, L, B, Nx, Ny, cap, n, Δ, BC, inside, empty, ns_vec, b_left, b_bottom, b_right, b_top)
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
    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+n] = A3^2 / W3
    end
    @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-n] = A1^2 / W1
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII+1] = A4^2 / W4
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2, W1, W2, W3, W4 = get_capacities(cap, II, Δ)

        @inbounds L[pII,pII-1] = A2^2 / W2
    end

    ns_vec .= 1.
    @inbounds @threads for II in empty
        pII = lexicographic(II, n)
        if sum(abs.(L[pII,:])) <= 1e-10
            @inbounds L[pII,pII] = -4.0
            @inbounds ns_vec[pII] = 0.
        end
    end
    ns_vec ./= norm(ns_vec)

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
        @inbounds O[pII,pJJ] += -(B1[II]*Δ - B2[II]*Δ)
    end
    return nothing
end

@inline function set_div_bnd!(::Dirichlet, ::Periodic, fun, O, B1, B2, np, nuv, Δ, b_indices)
    return nothing
end

# @inline function divergence_boundaries(::Dirichlet, Ox, Oy, Dx, Dy, cap, capu, capv, n, Δ, BCu, BCv, b_left, b_bottom, b_right, b_top)
#     Ax = hcat(cap[:,:,1], cap[:,end,3])
#     Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

#     mask = findall(iszero, Ax)
#     capx = copy(capu)
#     capx[mask,6] .= 0.
#     mask = findall(iszero, Ay)
#     capy = copy(capv)
#     capy[mask,7] .= 0.

#     set_div_bnd!(dir, BCu.left.t, x->x, Ox, capx[:,:,3], capx[:,:,6], n, n, Δ, b_left)
#     set_div_bnd!(dir, BCv.bottom.t, x->x, Oy, capy[:,:,4], capy[:,:,7], n, n+1, Δ, b_bottom)
#     set_div_bnd!(dir, BCu.right.t, δx⁺, Ox, capx[:,:,6], capx[:,:,1], n, n, Δ, b_right)
#     set_div_bnd!(dir, BCv.top.t, δy⁺, Oy, capy[:,:,7], capy[:,:,2], n, n+1, Δ, b_top)
# end

@inline function divergence_boundaries(::Dirichlet, Ox, Oy, Dx, Dy, cap, capu, capv, n, Δ, BCu, BCv, b_left, b_bottom, b_right, b_top)
    set_div_bnd!(dir, BCu.left.t, x->x, Ox, cap[:,:,6], cap[:,:,1], n, n, Δ, b_left)
    set_div_bnd!(dir, BCv.bottom.t, x->x, Oy, cap[:,:,7], cap[:,:,2], n, n+1, Δ, b_bottom)
    set_div_bnd!(dir, BCu.right.t, δx⁺, Ox, cap[:,:,3], cap[:,:,6], n, n, Δ, b_right)
    set_div_bnd!(dir, BCv.top.t, δy⁺, Oy, cap[:,:,4], cap[:,:,7], n, n+1, Δ, b_top)
end

# function divergence!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, cap, capu, capv, n, Δ, all_indices)
#     Ax = hcat(cap[:,:,1], cap[:,end,3])
#     Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

#     mask = findall(iszero, Ax)
#     capx = copy(capu)
#     capx[mask,6] .= 0.
#     mask = findall(iszero, Ay)
#     capy = copy(capv)
#     capy[mask,7] .= 0.

#     Bx .= 0.0
#     By .= 0.0
#     @inbounds @threads for II in all_indices
#         pII_1 = lexicographic(II, n)
#         pII_2 = lexicographic(δx⁺(II), n)
#         A3_1, B1_1, A1_2, B1_2 = get_capacities_x(capx, II, δx⁺(II), Δ)

#         pJJ_1 = lexicographic(II, n+1)
#         pJJ_2 = lexicographic(δy⁺(II), n+1)
#         A4_1, B2_1, A2_2, B2_2 = get_capacities_y(capy, II, δy⁺(II), Δ)
        
#         @inbounds Ox[pII_1,pII_1] = -B1_1
#         @inbounds Ox[pII_1,pII_2] = B1_2

#         @inbounds Oy[pII_1,pJJ_1] = -B2_1
#         @inbounds Oy[pII_1,pJJ_2] = B2_2

#         @inbounds Bx[pII_1] += -(B1_2 - A1_2) * Dx[δx⁺(II)]
#         @inbounds Bx[pII_1] += -(A3_1 - B1_1) * Dx[II]
#         @inbounds By[pII_1] += -(B2_2 - A2_2) * Dy[δy⁺(II)]
#         @inbounds By[pII_1] += -(A4_1 - B2_1) * Dy[II]
#     end

#     return nothing
# end

function divergence!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, cap, capu, capv, n, Δ, all_indices)
    Bx .= 0.0
    By .= 0.0
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        A1, A3, B1 = get_capacities_x(cap, II, Δ)

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        A2, A4, B2 = get_capacities_y(cap, II, Δ)
        
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

@inline function set_grad_bnd!(::Neumann, ::Dirichlet, O, B, nuv, np, Δ, b_indices, b_periodic)
    @error ("Not implemented yet.\nTry Neumann or Periodic in the outer BCs.")
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Neumann, O, B, nuv, np, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_grad_bnd!(::Neumann, ::Periodic, O, B, nuv, np, Δ, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, nuv)
        pJJ = lexicographic(JJ, np)
        @inbounds O[pII,pJJ] = B[II]*Δ
    end
    return nothing
end

# function gradient!(::Neumann, Ox, Oy, Divx, Divy, cap, capu, capv, n, Δ, BC, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
#     Ax = hcat(cap[:,:,1], cap[:,end,3])
#     Ay = vcat(cap[:,:,2], cap[newaxis,end,:,4])

#     mask = findall(iszero, Ax)
#     capx = copy(capu)
#     capx[mask,6] .= 0.
#     mask = findall(iszero, Ay)
#     capy = copy(capv)
#     capy[mask,7] .= 0.

#     Ox .= -Divx'
#     Oy .= -Divy'

#     set_grad_bnd!(neu, BC.left.t, x->x, Ox, -capx[:,:,6], n, n, Δ, b_left_u, b_right_p)
#     set_grad_bnd!(neu, BC.bottom.t, x->x, Oy, -capy[:,:,7], n+1, n, Δ, b_bottom_v, b_top_p)
#     set_grad_bnd!(neu, BC.right.t, δx⁺, Ox, capx[:,:,6], n, n, Δ, b_right_u, b_left_p)
#     set_grad_bnd!(neu, BC.top.t, δy⁺, Oy, capy[:,:,7], n+1, n, Δ, b_top_v, b_bottom_p)

#     return nothing
# end

function gradient!(::Neumann, Ox, Oy, Divx, Divy, cap, capu, capv, n, Δ, BC, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
    mat_T_op!(Ox, Divx, x->-x)
    mat_T_op!(Oy, Divy, x->-x)

    set_grad_bnd!(neu, BC.left.t, Ox, -cap[:,:,1], n, n, Δ, b_left_u, b_right_p)
    set_grad_bnd!(neu, BC.bottom.t, Oy, -cap[:,:,2], n+1, n, Δ, b_bottom_v, b_top_p)
    set_grad_bnd!(neu, BC.right.t, Ox, cap[:,:,3], n, n, Δ, b_right_u, b_left_p)
    set_grad_bnd!(neu, BC.top.t, Oy, cap[:,:,4], n+1, n, Δ, b_top_v, b_bottom_p)

    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Dirichlet, fun, O, B1, B2, nuv, np, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Neumann, fun, O, B1, B2, nuv, np, Δ, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, np)
        pJJ = lexicographic(fun(II), nuv)
        @inbounds O[pJJ,pII] += -(B1[II] - B2[II])*Δ
    end
    return nothing
end

@inline function set_grad_bnd!(::Dirichlet, ::Periodic, fun, O, B1, B2, nuv, np, Δ, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, np)
        pJJ = lexicographic(JJ, nuv)
        @inbounds O[pJJ,pII] = -(B1[II] - B2[II])*Δ
    end
    return nothing
end

function gradient!(::Dirichlet, Ox, Oy, Bx, By, Dx, Dy, cap, n, Δ, BC, all_indices, b_left, b_bottom, b_right, b_top)
    Bx .= 0.0
    By .= 0.0
    @inbounds @threads for II in all_indices
        pII_1 = lexicographic(II, n)
        pII_2 = lexicographic(δx⁺(II), n)
        A1, A3, B1 = get_capacities_x(cap, II, Δ)

        pJJ_1 = lexicographic(II, n+1)
        pJJ_2 = lexicographic(δy⁺(II), n+1)
        A2, A4, B2 = get_capacities_y(cap, II, Δ)
        
        @inbounds Ox[pII_1,pII_1] = B1
        @inbounds Ox[pII_2,pII_1] = -B1

        @inbounds Oy[pJJ_1,pII_1] = B2
        @inbounds Oy[pJJ_2,pII_1] = -B2

        @inbounds Bx[pII_1] += -(B1 - A1) * Dx[II]
        @inbounds Bx[pII_2] += -(A3 - B1) * Dx[II]
        @inbounds By[pJJ_1] += -(B2 - A2) * Dy[II]
        @inbounds By[pJJ_2] += -(A4 - B2) * Dy[II]
    end

    set_grad_bnd!(dir, BC.left.t, x->x, Ox, cap[:,:,6], cap[:,:,1], n, n, Δ, b_left, b_right)
    set_grad_bnd!(dir, BC.bottom.t, x->x, Oy, cap[:,:,7], cap[:,:,2], n+1, n, Δ, b_bottom, b_top)
    set_grad_bnd!(dir, BC.right.t, δx⁺, Ox, cap[:,:,3], cap[:,:,6], n, n, Δ, b_right, b_left)
    set_grad_bnd!(dir, BC.top.t, δy⁺, Oy, cap[:,:,4], cap[:,:,7], n+1, n, Δ, b_top, b_bottom)

    return nothing
end

function face_to_cell_gradient!(::Dirichlet, Ox, Oy, Gx, Gy, cap, n, all_indices)
    @inbounds @threads for II in all_indices
        pII = lexicographic(II, n)
        pJJ = lexicographic(II, n+1)
        A1, A2, A3, A4, V, B1, B2 = cap[II,1:7]

        @inbounds Ox[pII,pII] = B1 - A1
        @inbounds Ox[pII,pII+n] = A3 - B1

        @inbounds Oy[pII,pJJ] = B2 - A2
        @inbounds Oy[pII,pJJ+1] = A4 - B2
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
function strain_rate!(::Dirichlet, O11, O12_x, O12_y, O22, cap_x, cap_y, n, Δ, all_indices, inside)
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)

        JJ = δx⁺(II)
        pJJ = lexicographic(JJ, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_x, II, Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_x, JJ, Δ)

        @inbounds O11[pII,pII] = -B1_1 / W3_1
        @inbounds O11[pII,pJJ] = B1_2 / W3_1

        JJ = δy⁺(II)
        _pII = lexicographic(II, n+1)
        pJJ = lexicographic(JJ, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1, W1_1, W2_1, W3_1, W4_1 = get_capacities(cap_y, II, Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2, W1_2, W2_2, W3_2, W4_2 = get_capacities(cap_y, JJ, Δ)

        @inbounds O22[pII,_pII] = -B2_1 / W4_1
        @inbounds O22[pII,pJJ] = B2_2 / W4_1
    end
    @inbounds @threads for II in all_indices[1:end-1,2:end]
        pII = lexicographic(II, n)

        JJ_1 = II
        JJ_2 = δy⁺(II)
        pJJ_1 = lexicographic(JJ_1, n)
        pJJ_2 = lexicographic(JJ_2, n)
        A1_1_x, A2_1_x, A3_1_x, A4_1_x, B1_1_x, B2_1_x, W1_1_x, W2_1_x, W3_1_x, W4_1_x = get_capacities(cap_x, JJ_1, Δ)
        A1_2_x, A2_2_x, A3_2_x, A4_2_x, B1_2_x, B2_2_x, W1_2_x, W2_2_x, W3_2_x, W4_2_x = get_capacities(cap_x, JJ_2, Δ)

        @inbounds O12_x[pII,pJJ_1] = -B2_1_x / 2W4_1_x
        @inbounds O12_x[pII,pJJ_2] = B2_2_x / 2W4_1_x

        KK_1 = δy⁺(δx⁻(II))
        KK_2 = δy⁺(II)
        pKK_1 = lexicographic(KK_1, n+1)
        pKK_2 = lexicographic(KK_2, n+1)
        A1_1_y, A2_1_y, A3_1_y, A4_1_y, B1_1_y, B2_1_y, W1_1_y, W2_1_y, W3_1_y, W4_1_y = get_capacities(cap_y, KK_1, Δ)
        A1_2_y, A2_2_y, A3_2_y, A4_2_y, B1_2_y, B2_2_y, W1_2_y, W2_2_y, W3_2_y, W4_2_y = get_capacities(cap_y, KK_2, Δ)

        Ŵ = harmonic_average(W4_1_x, W3_1_y)

        @inbounds O12_y[pII,pKK_1] = -B1_1_y / 2Ŵ
        @inbounds O12_y[pII,pKK_2] = B1_2_y / 2Ŵ
    end

    return nothing
end

function set_bc_bnds(::Dirichlet, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
    Dx = copy(Du)
    Dy = copy(Dv)

    if is_neumann(BC_u.left.t)
        @inbounds Dx[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
    elseif is_dirichlet(BC_u.left.t)
        @inbounds Dx[:,1] .= BC_u.left.val
    elseif is_periodic(BC_u.left.t)
        @inbounds Dx[:,1] .= u[:,end]
    end
    if is_neumann(BC_v.bottom.t)
        @inbounds Dy[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
    elseif is_dirichlet(BC_v.bottom.t)
        @inbounds Dy[1,:] .= BC_v.bottom.val
    elseif is_periodic(BC_v.bottom.t)
        @inbounds Dy[1,:] .= v[end,:]
    end
    if is_neumann(BC_u.right.t)
        @inbounds Dx[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
    elseif is_dirichlet(BC_u.right.t)
        @inbounds Dx[:,end] .= BC_u.right.val
    elseif is_periodic(BC_u.right.t)
        @inbounds Dx[:,end] .= u[:,1]
    end
    if is_neumann(BC_v.top.t)
        @inbounds Dy[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
    elseif is_dirichlet(BC_v.top.t)
        @inbounds Dy[end,:] .= BC_v.top.val
    elseif is_periodic(BC_v.top.t)
        @inbounds Dy[end,:] .= v[1,:]
    end

    return Dx, Dy
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun, A1, A2, B1, D, n, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Neumann, O, fun, A1, A2, B1, D, n, Δ, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds O[pII,pII] += -0.5 * ((A2[II] - B1[II])*Δ * D[fun(II)] + (B1[II] - A1[II])*Δ * D[II])
    end
    return nothing
end

@inline function set_sca_conv_bnd!(::Dirichlet, ::Periodic, O, fun, A1, A2, B1, D, n, Δ, b_indices, b_periodic)
    @inbounds for (II, JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pII,pJJ] += -0.5 * ((A2[II] - B1[II])*Δ * D[fun(II)] + (B1[II] - A1[II])*Δ * D[II])
    end
    return nothing
end

function scalar_convection!(::Dirichlet, O, B, u, v, Dx, Dy, Du, Dv, cap, n, Δ, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
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

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
        u1, v2, u3, v4 = u[II], v[II], u[δx⁺(II)], v[δy⁺(II)]

        @inbounds O[pII,pII] = 0.5 * (A3 * u3 - A1 * u1 + A4 * v4 - A2 * v2)

        @inbounds O[pII,pII] += -0.5 * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds O[pII,pII] += -0.5 * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])

        @inbounds B[pII] += -0.5 * Dx[II] * ((A3 - B1) * Du[δx⁺(II)] + (B1 - A1) * Du[II])
        @inbounds B[pII] += -0.5 * Dy[II] * ((A4 - B2) * Dv[δy⁺(II)] + (B2 - A2) * Dv[II])
    end
    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)

        @inbounds O[pII,pII+n] = 0.5 * A3 * u[δx⁺(II)]
    end
    @inbounds @threads for II in vcat(b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)

        @inbounds O[pII,pII-n] = -0.5 * A1 * u[II]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_bottom, b_right[2:end-1])
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)

        @inbounds O[pII,pII+1] = 0.5 * A4 * v[δy⁺(II)]
    end
    @inbounds @threads for II in vcat(b_left[2:end-1], b_right[2:end-1], b_top)
        pII = lexicographic(II, n)
        A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)

        @inbounds O[pII,pII-1] = -0.5 * A2 * v[II]
    end

    if is_periodic(BC.left.t) && is_periodic(BC.right.t)
        @inbounds for (II,JJ) in zip(b_right, b_left)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
    
            @inbounds O[pII,pJJ] = 0.5 * A3 * u[δx⁺(II)]
        end
        @inbounds for (II,JJ) in zip(b_left, b_right)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
    
            @inbounds O[pII,pJJ] = -0.5 * A1 * u[II]
        end
    end
    if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
        @inbounds for (II,JJ) in zip(b_top, b_bottom)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
    
            @inbounds O[pII,pJJ] = 0.5 * A4 * v[δy⁺(II)]
        end
        @inbounds for (II,JJ) in zip(b_bottom, b_top)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1, A2, A3, A4, B1, B2 = get_capacities_convection(cap, II, Δ)
    
            @inbounds O[pII,pJJ] = -0.5 * A2 * v[II]
        end
    end

    @inbounds _A1 = cap[:,:,1]
    @inbounds _A2 = cap[:,:,2]
    @inbounds _A3 = cap[:,:,3]
    @inbounds _A4 = cap[:,:,4]
    @inbounds _B1 = cap[:,:,6]
    @inbounds _B2 = cap[:,:,7]

    set_sca_conv_bnd!(dir, BC.left.t, O, δx⁺, _A1, _A3, _B1, Du, n, Δ, b_left, b_right)
    set_sca_conv_bnd!(dir, BC.bottom.t, O, δy⁺, _A2, _A4, _B2, Dv, n, Δ, b_bottom, b_top)
    set_sca_conv_bnd!(dir, BC.right.t, O, δx⁺, _A1, _A3, _B1, Du, n, Δ, b_right, b_left)
    set_sca_conv_bnd!(dir, BC.top.t, O, δy⁺, _A2, _A4, _B2, Dv, n, Δ, b_top, b_bottom)

    return nothing
end

function set_bc_bnds(::Dirichlet, ::Union{GridFCx,GridFCy}, Du, Dv, Hu, Hv, u, v, BC_u, BC_v)
    Du1_x = copy(Du)
    Du1_y = copy(Du)
    Du2_x = copy(Du)
    Du2_y = copy(Du)
    Dv_x = copy(Dv)
    Dv_y = copy(Dv)

    if is_neumann(BC_u.left.t)
        @inbounds Du1_x[:,1] .= u[:,1] .+ Hu[:,1] .* BC_u.left.val
        @inbounds Du1_x[:,2] .= u[:,2]
        @inbounds Du2_x[:,1] .= Hu[:,1] .* BC_u.left.val
        @inbounds Du2_x[:,2] .= u[:,2]
    elseif is_dirichlet(BC_u.left.t)
        @inbounds Du1_x[:,1] .= BC_u.left.val
        @inbounds Du1_x[:,2] .= BC_u.left.val
        @inbounds Du2_x[:,1] .= BC_u.left.val
        @inbounds Du2_x[:,2] .= BC_u.left.val
    elseif is_periodic(BC_u.left.t)
        @inbounds Du1_x[:,1] .= u[:,end]
        @inbounds Du1_x[:,2] .= u[:,2]
        @inbounds Du2_x[:,1] .= 0.0
        @inbounds Du2_x[:,2] .= u[:,2]
    end
    if is_neumann(BC_u.bottom.t)
        @inbounds Du1_y[1,:] .= u[1,:] .+ Hu[1,:] .* BC_u.bottom.val
        @inbounds Du1_y[2,:] .= u[2,:]
        @inbounds Du2_y[1,:] .= Hu[1,:] .* BC_u.bottom.val
        @inbounds Du2_y[2,:] .= u[2,:]
    elseif is_dirichlet(BC_u.bottom.t)
        @inbounds Du1_y[1,:] .= BC_u.bottom.val
        @inbounds Du1_y[2,:] .= BC_u.bottom.val
        @inbounds Du2_y[1,:] .= BC_u.bottom.val
        @inbounds Du2_y[2,:] .= BC_u.bottom.val
    elseif is_periodic(BC_u.bottom.t)
        @inbounds Du1_y[1,:] .= u[end,:]
        @inbounds Du1_y[2,:] .= u[2,:]
        @inbounds Du2_y[1,:] .= 0.0
        @inbounds Du2_y[2,:] .= u[2,:]
    end
    if is_neumann(BC_u.right.t)
        @inbounds Du1_x[:,end] .= u[:,end] .+ Hu[:,end] .* BC_u.right.val 
        @inbounds Du1_x[:,end-1] .= u[:,end-1]
        @inbounds Du2_x[:,end] .= Hu[:,end] .* BC_u.right.val
        @inbounds Du2_x[:,end-1] .= u[:,end-1]
    elseif is_dirichlet(BC_u.right.t)
        @inbounds Du1_x[:,end] .= BC_u.right.val
        @inbounds Du1_x[:,end-1] .= BC_u.right.val
        @inbounds Du2_x[:,end] .= BC_u.right.val
        @inbounds Du2_x[:,end-1] .= BC_u.right.val
    elseif is_periodic(BC_u.right.t)
        @inbounds Du1_x[:,end] .= u[:,1]
        @inbounds Du1_x[:,end-1] .= u[:,end-1]
        @inbounds Du2_x[:,end] .= 0.0
        @inbounds Du2_x[:,end-1] .= u[:,end-1]
    end
    if is_neumann(BC_u.top.t)
        @inbounds Du1_y[end,:] .= u[end,:] .+ Hu[end,:] .* BC_u.top.val
        @inbounds Du1_y[end-1,:] .= u[end-1,:]
        @inbounds Du2_y[end,:] .= Hu[end,:] .* BC_u.top.val
        @inbounds Du2_y[end-1,:] .= u[end-1,:]
    elseif is_dirichlet(BC_u.top.t)
        @inbounds Du1_y[end,:] .= BC_u.top.val
        @inbounds Du1_y[end-1,:] .= BC_u.top.val
        @inbounds Du2_y[end,:] .= BC_u.top.val
        @inbounds Du2_y[end-1,:] .= BC_u.top.val
    elseif is_periodic(BC_u.top.t)
        @inbounds Du1_y[end,:] .= u[1,:]
        @inbounds Du1_y[end-1,:] .= u[end-1,:]
        @inbounds Du2_y[end,:] .= 0.0
        @inbounds Du2_y[end-1,:] .= u[end-1,:]
    end

    if is_neumann(BC_v.left.t)
        @inbounds Dv_x[:,1] .= v[:,1] .+ Hv[:,1] .* BC_v.left.val 
        @inbounds Dv_x[:,2] .= v[:,2]
    elseif is_dirichlet(BC_v.left.t)
        @inbounds Dv_x[:,1] .= BC_v.left.val
        @inbounds Dv_x[:,2] .= BC_v.left.val
    elseif is_periodic(BC_v.left.t)
        @inbounds Dv_x[:,1] .= v[:,end]
        @inbounds Dv_x[:,2] .= v[:,2]
    end
    if is_neumann(BC_v.bottom.t)
        @inbounds Dv_y[1,:] .= v[1,:] .+ Hv[1,:] .* BC_v.bottom.val 
        @inbounds Dv_y[2,:] .= v[2,:]
    elseif is_dirichlet(BC_v.bottom.t)
        @inbounds Dv_y[1,:] .= BC_v.bottom.val
        @inbounds Dv_y[2,:] .= BC_v.bottom.val
    elseif is_periodic(BC_v.bottom.t)
        @inbounds Dv_y[1,:] .= v[end,:]
        @inbounds Dv_y[2,:] .= v[2,:]
    end
    if is_neumann(BC_v.right.t)
        @inbounds Dv_x[:,end] .= v[:,end] .+ Hv[:,end] .* BC_v.right.val 
        @inbounds Dv_x[:,end-1] .= v[:,end-1]
    elseif is_dirichlet(BC_v.right.t)
        @inbounds Dv_x[:,end] .= BC_v.right.val
        @inbounds Dv_x[:,end-1] .= BC_v.right.val
    elseif is_periodic(BC_v.right.t)
        @inbounds Dv_x[:,end] .= v[:,1]
        @inbounds Dv_x[:,end-1] .= v[:,end-1]
    end
    if is_neumann(BC_v.top.t)
        @inbounds Dv_y[end,:] .= v[end,:] .+ Hv[end,:] .* BC_v.top.val 
        @inbounds Dv_y[end-1,:] .= v[end-1,:]
    elseif is_dirichlet(BC_v.top.t)
        @inbounds Dv_y[end,:] .= BC_v.top.val
        @inbounds Dv_y[end-1,:] .= BC_v.top.val
    elseif is_periodic(BC_v.top.t)
        @inbounds Dv_y[end-1,:] .= v[end-1,:]
    end

    return Du1_x, Du1_y, Du2_x, Du2_y, Dv_x, Dv_y
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Dirichlet, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, Δ, b_indices, b_periodic)
    return nothing
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Neumann, O, fun_cap, fun1, fun2, fun3, fun4, A3, B1, A1, A4, B2, A2, Du, Dv, n, Δ, b_indices, b_periodic)
    @inbounds @threads for II in b_indices
        pII = lexicographic(II, n)
        @inbounds O[pII,pII] += -0.25 * (A3[fun_cap(II)] - B1[fun_cap(II)])*Δ * Du[fun1(II)]
        @inbounds O[pII,pII] += -0.25 * (B1[fun_cap(II)] - A1[fun_cap(II)])*Δ * Du[fun2(II)]
        @inbounds O[pII,pII] += -0.25 * (A4[fun_cap(II)] - B2[fun_cap(II)])*Δ * Dv[fun3(II)]
        @inbounds O[pII,pII] += -0.25 * (B2[fun_cap(II)] - A2[fun_cap(II)])*Δ * Dv[fun4(II)]
    end
    return nothing
end

@inline function set_vec_conv_bnd!(::Dirichlet, ::Periodic, O, fun_cap, fun1, fun2, fun3, fun4, A1, A2, A3, A4, B1, B2, Du, Dv, n, Δ, b_indices, b_periodic)
    @inbounds for (II,JJ) in zip(b_indices, b_periodic)
        pII = lexicographic(II, n)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pII,pJJ] = -0.25 * (A1[fun_cap(II)] - B1[fun_cap(II)])*Δ * D[fun1(II)]
        @inbounds O[pII,pJJ] = -0.25 * (B1[fun_cap(II)] - A2[fun_cap(II)])*Δ * D[fun2(II)]
    end
    return nothing
end

function fill_inside_conv!(::GridFCx, O, B, u, v, Du1_x, Du1_y, Dv_y, cap, n, Δ, II)
    pII = lexicographic(II, n)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)

    Auim1, Aui, Auip1 = A1_1 * u[δx⁻(II)], A1_2 * u[II], A3_2 * u[δx⁺(II)]
    Avim1jm1, Avim1jp1 = A2_1 * v[δx⁻(II)], A4_1 * v[δy⁺(δx⁻(II))]
    Avip1jm1, Avip1jp1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

    Au1 = 0.5 * (Auim1 + Aui)
    Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    Au3 = 0.5 * (Aui + Auip1)
    Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

    @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
    @inbounds O[pII,pII+n] = 0.5 * Au3
    @inbounds O[pII,pII-n] = -0.5 * Au1
    @inbounds O[pII,pII+1] = 0.5 * Au4
    @inbounds O[pII,pII-1] = -0.5 * Au2

    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - B1_1) * Du1_x[II]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv_y[II]

    # @inbounds B[pII] += -0.25 * Du2_x[II] * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
    # @inbounds B[pII] += -0.25 * Du2_x[II] * (B1_2 - B1_1) * Du1_x[II]
    # @inbounds B[pII] += -0.25 * Du2_x[II] * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

    # @inbounds B[pII] += -0.25 * Du2_y[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
    # @inbounds B[pII] += -0.25 * Du2_y[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
    # @inbounds B[pII] += -0.25 * Du2_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    # @inbounds B[pII] += -0.25 * Du2_y[II] * (B2_2 - A2_2) * Dv_y[II]

    @inbounds B[pII] += -0.25 * Du1_x[II] * (A3_2 - B1_2) * Du1_x[δx⁺(II)]
    @inbounds B[pII] += -0.25 * Du1_x[II] * (B1_2 - B1_1) * Du1_x[II]
    @inbounds B[pII] += -0.25 * Du1_x[II] * (B1_1 - A1_1) * Du1_x[δx⁻(II)]

    @inbounds B[pII] += -0.25 * Du1_y[II] * (A4_1 - B2_1) * Dv_y[δx⁻(δy⁺(II))]
    @inbounds B[pII] += -0.25 * Du1_y[II] * (B2_1 - A2_1) * Dv_y[δx⁻(II)]
    @inbounds B[pII] += -0.25 * Du1_y[II] * (A4_2 - B2_2) * Dv_y[δy⁺(II)]
    @inbounds B[pII] += -0.25 * Du1_y[II] * (B2_2 - A2_2) * Dv_y[II]

    return nothing
end

function vector_convection!(::Dirichlet, ::GridFCx, O, B, u, v, Du1_x, Du1_y, Du2_x, Du2_y, Dv_x, Dv_y, cap, n, Δ, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(gfcx, O, B, u, v, Du1_x, Du1_y, Dv_y, cap, n, Δ, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left, b_bottom[2:end-1], b_top[2:end-1])
    bc = ((Du1_x, Du1_x, Dv_x), (Du1_y, Du1_y, Dv_y), (Du1_y, Du1_y, Dv_y))
    for (bnd, (Du1, Du2, Dv)) in zip(bnds, bc)
        @inbounds for II in bnd
            pII = lexicographic(II, n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
            
            Aui, Auip1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

            Au3 = 0.5 * (Aui + Auip1)

            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pII+n] = 0.5 * Au3

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1[δx⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du1[II]

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

            @inbounds B[pII] += -0.25 * Du2[II] * (A3_2 - B1_2) * Du1[δx⁺(II)]
            @inbounds B[pII] += -0.25 * Du2[II] * (B1_2 - A1_2) * Du1[II]

            @inbounds B[pII] += -0.25 * Du2[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
            @inbounds B[pII] += -0.25 * Du2[II] * (B2_2 - A2_2) * Dv[II]
        end
    end
    bnds = (b_bottom[2:end-1], b_right, b_top[2:end-1])
    bc = ((Du1_y, Du1_y, Dv_y), (Du1_x, Du1_x, Dv_x), (Du1_y, Du1_y, Dv_y))
    for (bnd, (Du1, Du2, Dv)) in zip(bnds, bc)
        @inbounds for II in bnd
            pII = lexicographic(II, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
            
            Auim1, Aui = A1_1 * u[δx⁻(II)], A3_1 * u[II]

            Au1 = 0.5 * (Auim1 + Aui)

            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pII-n] = -0.5 * Au1

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du1[II]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1[δx⁻(II)]

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[δx⁻(II)]

            @inbounds B[pII] += -0.25 * Du2[II] * (A3_1 - B1_1) * Du1[II]
            @inbounds B[pII] += -0.25 * Du2[II] * (B1_1 - A1_1) * Du1[δx⁻(II)]

            @inbounds B[pII] += -0.25 * Du2[II] * (A4_1 - B2_1) * Dv[δx⁻(δy⁺(II))]
            @inbounds B[pII] += -0.25 * Du2[II] * (B2_1 - A2_1) * Dv[δx⁻(II)]
        end
    end
    @inbounds for II in b_bottom[2:end-1]
        pII = lexicographic(II, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]
        Avip1jp1 = A4_2 * v[δy⁺(II)]

        Au4 = 0.5 * (Avim1jp1 + Avip1jp1)

        @inbounds O[pII,pII] += 0.5 * Au4
        @inbounds O[pII,pII+1] = 0.5 * Au4
    end
    @inbounds for II in b_top[2:end-1]
        pII = lexicographic(II, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Avim1jm1 = A2_1 * v[δx⁻(II)]
        Avip1jm1 = A2_2 * v[II]

        Au2 = 0.5 * (Avim1jm1 + Avip1jm1)

        @inbounds O[pII,pII] += -0.5 * Au2
        @inbounds O[pII,pII-1] = -0.5 * Au2
    end
    @inbounds for II in b_left[2:end]
        pII = lexicographic(II, n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Avip1jm1 = A2_2 * v[II]

        Au2 = 0.5 * Avip1jm1

        if is_periodic(BC.left.t)
            JJ = II + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Avim1jm1 = A2_1 * v[δx⁻(JJ)]
            Au2 += 0.5 * Avim1jm1
        end

        @inbounds O[pII,pII] += -0.5 * Au2
        @inbounds O[pII,pII-1] = -0.5 * Au2
    end
    @inbounds for II in b_left[1:end-1]
        pII = lexicographic(II, n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Avip1jp1 = A4_2 * v[δy⁺(II)]

        Au4 = 0.5 * Avip1jp1

        if is_periodic(BC.left.t)
            JJ = II + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
            Au4 += 0.5 * Avim1jp1
        end

        @inbounds O[pII,pII] += 0.5 * Au4
        @inbounds O[pII,pII+1] = 0.5 * Au4
    end
    @inbounds for II in b_right[2:end]
        pII = lexicographic(II, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
        
        Avim1jm1 = A2_1 * v[δx⁻(II)]

        Au2 = 0.5 * Avim1jm1

        if is_periodic(BC.right.t)
            JJ = II + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Avip1jm1 = A2_2 * v[JJ]
            Au2 += 0.5 * Avip1jm1
        end

        @inbounds O[pII,pII] += -0.5 * Au2
        @inbounds O[pII,pII-1] = -0.5 * Au2
    end
    @inbounds for II in b_right[1:end-1]
        pII = lexicographic(II, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
        
        Avim1jp1 = A4_1 * v[δy⁺(δx⁻(II))]

        Au4 = 0.5 * Avim1jp1

        if is_periodic(BC.right.t)
            JJ = II + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]
            Au4 += 0.5 * Avip1jp1
        end

        @inbounds O[pII,pII] += 0.5 * Au4
        @inbounds O[pII,pII+1] = 0.5 * Au4
    end

    if is_periodic(BC.left.t) && is_periodic(BC.right.t)
        (Du1, Du2, Dv) = (Du1_x, Du1_x, Dv_x)
        @inbounds for (II,JJ) in zip(b_left, b_right)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(JJ), Δ)
            
            Auim1, Aui = A1_1 * u[JJ], A3_1 * u[II]

            Au1 = 0.5 * (Auim1 + Aui)

            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du1[II]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du1[JJ]

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv[δy⁺(JJ)]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv[JJ]

            @inbounds B[pII] += -0.25 * Du2[II] * (A3_1 - B1_1) * Du1[II]
            @inbounds B[pII] += -0.25 * Du2[II] * (B1_1 - A1_1) * Du1[JJ]

            @inbounds B[pII] += -0.25 * Du2[II] * (A4_1 - B2_1) * Dv[δy⁺(JJ)]
            @inbounds B[pII] += -0.25 * Du2[II] * (B2_1 - A2_1) * Dv[JJ]
        end
        (Du1, Du2, Dv) = (Du1_x, Du1_x, Dv_x)
        @inbounds for (II, JJ) in zip(b_right, b_left)
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            
            Aui, Auip1 = A1_2 * u[II], A3_2 * u[JJ]

            Au3 = 0.5 * (Aui + Auip1)

            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du1[JJ]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du1[II]

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv[δy⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv[II]

            @inbounds B[pII] += -0.25 * Du2[II] * (A3_2 - B1_2) * Du1[JJ]
            @inbounds B[pII] += -0.25 * Du2[II] * (B1_2 - A1_2) * Du1[II]

            @inbounds B[pII] += -0.25 * Du2[II] * (A4_2 - B2_2) * Dv[δy⁺(II)]
            @inbounds B[pII] += -0.25 * Du2[II] * (B2_2 - A2_2) * Dv[II]
        end
    end
    if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
        @inbounds for (II,JJ) in zip(b_bottom[2:end-1], b_top[2:end-1])
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(II), Δ)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
            
            Avim1jm1 = A2_1 * v[δx⁻(II)]
            Avip1jm1 = A2_2 * v[II]
    
            Au2 = 0.5 * (Avim1jm1 + Avip1jm1)
    
            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2
        end
        @inbounds for (II,JJ) in zip(b_top[2:end-1], b_bottom[2:end-1])
            pII = lexicographic(II, n)
            pJJ = lexicographic(JJ, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(δx⁻(II)), Δ)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, δy⁻(II), Δ)
            
            Avim1jp1 = A4_1 * v[δx⁻(JJ)]
            Avip1jp1 = A4_2 * v[JJ]

            Au4 = 0.5 * (Avim1jp1 + Avip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4
        end

        ii = b_left[1]
        pii = lexicographic(ii, n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii, Δ)
        
        Avip1jm1 = A2_2 * v[ii]

        Au2 = 0.5 * Avip1jm1

        if is_periodic(BC.left.t)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Avim1jm1 = A2_1 * v[δx⁻(JJ)]
            Au2 += 0.5 * Avim1jm1
        end

        JJ = ii + CartesianIndex(n-1, 0)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_left[end]
        pii = lexicographic(ii, n)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii, Δ)
        
        Avip1jp1 = A4_2 * v[δy⁺(ii)]

        Au4 = 0.5 * Avip1jp1

        if is_periodic(BC.left.t)
            JJ = ii + CartesianIndex(0, n)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Avim1jp1 = A4_1 * v[δy⁺(δx⁻(JJ))]
            Au4 += 0.5 * Avim1jp1
        end

        JJ = ii + CartesianIndex(-n+1, 0)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4

        ii = b_right[1]
        pii = lexicographic(ii, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii), Δ)
        
        Avim1jm1 = A2_1 * v[δx⁻(ii)]

        Au2 = 0.5 * Avim1jm1

        if is_periodic(BC.right.t)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Avip1jm1 = A2_2 * v[JJ]
            Au2 += 0.5 * Avip1jm1
        end

        JJ = ii + CartesianIndex(n-1, 0)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pii,pii] += -0.5 * Au2
        @inbounds O[pii,pJJ] = -0.5 * Au2
        
        ii = b_right[end]
        pii = lexicographic(ii, n)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δx⁻(ii), Δ)
        
        Avim1jp1 = A4_1 * v[δy⁺(δx⁻(ii))]

        Au4 = 0.5 * Avim1jp1

        if is_periodic(BC.right.t)
            JJ = ii + CartesianIndex(0, -n)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Avip1jp1 = A4_2 * v[δy⁺(JJ)]
            Au4 += 0.5 * Avip1jp1
        end
        
        JJ = ii + CartesianIndex(-n+1, 0)
        pJJ = lexicographic(JJ, n)
        @inbounds O[pii,pii] += 0.5 * Au4
        @inbounds O[pii,pJJ] = 0.5 * Au4
    end

    # @inbounds _A1 = cap[:,:,1]
    # @inbounds _A2 = cap[:,:,2]
    # @inbounds _A3 = cap[:,:,3]
    # @inbounds _A4 = cap[:,:,4]
    # @inbounds _B1 = cap[:,:,6]
    # @inbounds _B2 = cap[:,:,7]

    # set_vec_conv_bnd!(dir, BC.left.t, O, x->x, δx⁺, x->x,  δy⁺, x->x, _A3, _B1, _A1, _A4, _B2, _A2, Du1_x, Dv_x, n, Δ, b_left, b_right)
    # set_vec_conv_bnd!(dir, BC.bottom.t, O, x->x,  δx⁺, x->x, δy⁺, x->x, _A3, _B1, _A1, _A4, _B2, _A2, Du1_y, Dv_y, n, Δ, b_bottom[1:end-1], b_top[1:end-1])
    # set_vec_conv_bnd!(dir, BC.bottom.t, O, δx⁻, x->x, δx⁻, (δx⁻ ∘ δy⁺), δx⁻, _A3, _B1, _A1, _A4, _B2, _A2, Du1_y, Dv_y, n, Δ, b_bottom[2:end], b_top[2:end])
    # set_vec_conv_bnd!(dir, BC.right.t, O, δx⁻, x->x, δx⁻, (δx⁻ ∘ δy⁺), δx⁻, _A3, _B1, _A1, _A4, _B2, _A2, Du1_x, Dv_x, n, Δ, b_right, b_left)
    # set_vec_conv_bnd!(dir, BC.top.t, O, x->x, δx⁺, x->x, δy⁺, x->x, _A3, _B1, _A1, _A4, _B2, _A2, Du1_y, Dv_y, n, Δ, b_top[1:end-1], b_bottom[1:end-1])
    # set_vec_conv_bnd!(dir, BC.top.t, O, δx⁻, x->x, δx⁻, (δx⁻ ∘ δy⁺), δx⁻, _A3, _B1, _A1, _A4, _B2, _A2, Du1_y, Dv_y, n, Δ, b_top[2:end], b_bottom[2:end])

    return nothing
end

function fill_inside_conv!(::GridFCy, O, B, u, v, Du_x, Dv1_x, Dv1_y, cap, n, Δ, II)
    pII = lexicographic(II, n+1)
    A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
    A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
    
    Avim1, Avi, Avip1 = A2_1 * v[δy⁻(II)], A2_2 * v[II], A4_2 * v[δy⁺(II)]
    Auim1jm1, Auip1jm1 = A1_1 * u[δy⁻(II)], A3_1 * u[δx⁺(δy⁻(II))]
    Auim1jp1, Auip1jp1 = A1_2 * u[II], A3_2 * u[δx⁺(II)]

    Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    Au2 = 0.5 * (Avim1 + Avi)
    Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    Au4 = 0.5 * (Avi + Avip1)

    @inbounds O[pII,pII] = 0.5 * (Au3 - Au1 + Au4 - Au2)
    @inbounds O[pII,pII+n+1] = 0.5 * Au3
    @inbounds O[pII,pII-n-1] = -0.5 * Au1
    @inbounds O[pII,pII+1] = 0.5 * Au4
    @inbounds O[pII,pII-1] = -0.5 * Au2

    @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B2_2 - B2_1) * Dv1_y[II]
    @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

    @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
    @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du_x[δy⁻(II)]
    @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du_x[δx⁺(II)]
    @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du_x[II]

    # @inbounds B[pII] += -0.25 * Dv2_y[II] * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
    # @inbounds B[pII] += -0.25 * Dv2_y[II] * (B2_2 - B2_1) * Dv1_y[II]
    # @inbounds B[pII] += -0.25 * Dv2_y[II] * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

    # @inbounds B[pII] += -0.25 * Dv2_x[II] * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
    # @inbounds B[pII] += -0.25 * Dv2_x[II] * (B1_1 - A1_1) * Du_x[δy⁻(II)]
    # @inbounds B[pII] += -0.25 * Dv2_x[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
    # @inbounds B[pII] += -0.25 * Dv2_x[II] * (B1_2 - A1_2) * Du_x[II]

    @inbounds B[pII] += -0.25 * Dv1_y[II] * (A4_2 - B2_2) * Dv1_y[δy⁺(II)]
    @inbounds B[pII] += -0.25 * Dv1_y[II] * (B2_2 - B2_1) * Dv1_y[II]
    @inbounds B[pII] += -0.25 * Dv1_y[II] * (B2_1 - A2_1) * Dv1_y[δy⁻(II)]

    @inbounds B[pII] += -0.25 * Dv1_x[II] * (A3_1 - B1_1) * Du_x[δy⁻(δx⁺(II))]
    @inbounds B[pII] += -0.25 * Dv1_x[II] * (B1_1 - A1_1) * Du_x[δy⁻(II)]
    @inbounds B[pII] += -0.25 * Dv1_x[II] * (A3_2 - B1_2) * Du_x[δx⁺(II)]
    @inbounds B[pII] += -0.25 * Dv1_x[II] * (B1_2 - A1_2) * Du_x[II]
end

function vector_convection!(::Dirichlet, ::GridFCy, O, B, u, v, Du_x, Du_y, Dv1_x, Dv1_y, Dv2_x, Dv2_y, cap, n, Δ, BC, inside, b_left, b_bottom, b_right, b_top)
    B .= 0.0
    @inbounds @threads for II in inside
        fill_inside_conv!(gfcy, O, B, u, v, Du_x, Dv1_x, Dv1_y, cap, n, Δ, II)
    end

    @inbounds @threads for II in vcat(b_left, b_bottom[2:end-1], b_right, b_top[2:end-1])
        pII = lexicographic(II, n+1)
        @inbounds O[pII,pII] = 0.0
    end
    bnds = (b_left[2:end-1], b_bottom, b_right[2:end-1])
    bc = ((Du_x, Dv1_x, Dv1_x), (Du_y, Dv1_y, Dv1_y), (Du_x, Dv1_x, Dv1_x))
    for (bnd, (Du, Dv1, Dv2)) in zip(bnds, bc)
        @inbounds for II in bnd
            pII = lexicographic(II, n+1)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
            
            Avi, Avip1 = A2_2 * v[II], A4_2 * v[δy⁺(II)]

            Au4 = 0.5 * (Avi + Avip1)

            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pII+1] = 0.5 * Au4

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1[δy⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv1[II]

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A4_2 - B2_2) * Dv1[δy⁺(II)]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B2_2 - A2_2) * Dv1[II]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B1_2 - A1_2) * Du[II]
        end
    end
    bnds = (b_left[2:end-1], b_right[2:end-1], b_top)
    bc = ((Du_x, Dv1_x, Dv1_x), (Du_x, Dv1_x, Dv1_x), (Du_y, Dv1_y, Dv1_y))
    for (bnd, (Du, Dv1, Dv2)) in zip(bnds, bc)
        @inbounds for II in bnd
            pII = lexicographic(II, n+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
            
            Avim1, Avi = A2_1 * v[δy⁻(II)], A4_1 * v[II]

            Au2 = 0.5 * (Avim1 + Avi)

            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pII-1] = -0.5 * Au2

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv1[II]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1[δy⁻(II)]

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[δy⁻(II)]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A4_1 - B2_1) * Dv1[II]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B2_1 - A2_1) * Dv1[δy⁻(II)]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A3_1 - B1_1) * Du[δy⁻(δx⁺(II))]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B1_1 - A1_1) * Du[δy⁻(II)]
        end
    end
    @inbounds for II in b_left[2:end-1]
        pII = lexicographic(II, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
        Auip1jp1 = A3_2 * u[δx⁺(II)]

        Au3 = 0.5 * (Auip1jm1 + Auip1jp1)

        @inbounds O[pII,pII] += 0.5 * Au3
        @inbounds O[pII,pII+n+1] = 0.5 * Au3
    end
    @inbounds for II in b_right[2:end-1]
        pII = lexicographic(II, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Auim1jm1 = A1_1 * u[δy⁻(II)]
        Auim1jp1 = A1_2 * u[II]

        Au1 = 0.5 * (Auim1jm1 + Auim1jp1)

        @inbounds O[pII,pII] += -0.5 * Au1
        @inbounds O[pII,pII-n-1] = -0.5 * Au1
    end
    @inbounds for II in b_bottom[2:end]
        pII = lexicographic(II, n+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Auim1jp1 = A1_2 * u[II]

        Au1 = 0.5 * Auim1jp1

        if is_periodic(BC.bottom.t)
            JJ = II + CartesianIndex(n, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Auim1jm1 = A1_1 * u[JJ]
            Au1 += 0.5 * Auim1jm1
        end

        @inbounds O[pII,pII] += -0.5 * Au1
        @inbounds O[pII,pII-n-1] = -0.5 * Au1
    end
    @inbounds for II in b_bottom[1:end-1]
        pII = lexicographic(II, n+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
        
        Auip1jp1 = A3_2 * u[δx⁺(II)]

        Au3 = 0.5 * Auip1jp1

        if is_periodic(BC.bottom.t)
            JJ = II + CartesianIndex(n, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Auip1jm1 = A3_1 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jm1
        end

        @inbounds O[pII,pII] += 0.5 * Au3
        @inbounds O[pII,pII+n+1] = 0.5 * Au3
    end
    @inbounds for II in b_top[2:end]
        pII = lexicographic(II, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
        
        Auim1jm1 = A1_1 * u[δy⁻(II)]

        Au1 = 0.5 * Auim1jm1

        if is_periodic(BC.top.t)
            JJ = II + CartesianIndex(-n, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Auim1jp1 = A1_2 * u[JJ]
            Au1 += 0.5 * Auim1jp1
        end

        @inbounds O[pII,pII] += -0.5 * Au1
        @inbounds O[pII,pII-n-1] = -0.5 * Au1
    end
    @inbounds for II in b_top[1:end-1]
        pII = lexicographic(II, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
        
        Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]

        Au3 = 0.5 * Auip1jm1

        if is_periodic(BC.top.t)
            JJ = II + CartesianIndex(n, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Auip1jp1 = A3_2 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jp1
        end

        @inbounds O[pII,pII] += 0.5 * Au3
        @inbounds O[pII,pII+n+1] = 0.5 * Au3
    end

    if is_periodic(BC.bottom.t) && is_periodic(BC.top.t)
        (Du, Dv1, Dv2) = (Du_y, Dv1_y, Dv1_y)
        @inbounds for (II, JJ) in zip(b_bottom, b_top)
            pII = lexicographic(II, n+1)
            pJJ = lexicographic(JJ, n+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(JJ), Δ)
            
            Avim1, Avi = A2_1 * v[JJ], A4_1 * v[II]

            Au2 = 0.5 * (Avim1 + Avi)

            @inbounds O[pII,pII] += -0.5 * Au2
            @inbounds O[pII,pJJ] = -0.5 * Au2

            @inbounds O[pII,pII] += -0.25 * (A4_1 - B2_1) * Dv1[II]
            @inbounds O[pII,pII] += -0.25 * (B2_1 - A2_1) * Dv1[JJ]

            @inbounds O[pII,pII] += -0.25 * (A3_1 - B1_1) * Du[δx⁺(JJ)]
            @inbounds O[pII,pII] += -0.25 * (B1_1 - A1_1) * Du[JJ]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A4_1 - B2_1) * Dv1[II]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B2_1 - A2_1) * Dv1[JJ]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A3_1 - B1_1) * Du[δx⁺(JJ)]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B1_1 - A1_1) * Du[JJ]
        end
        (Du, Dv1, Dv2) = (Du_y, Dv1_y, Dv1_y)
        @inbounds for (II, JJ) in zip(b_top, b_bottom)
            pII = lexicographic(II, n+1)
            pJJ = lexicographic(JJ, n+1)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            
            Avi, Avip1 = A2_2 * v[II], A4_2 * v[JJ]

            Au4 = 0.5 * (Avi + Avip1)

            @inbounds O[pII,pII] += 0.5 * Au4
            @inbounds O[pII,pJJ] = 0.5 * Au4

            @inbounds O[pII,pII] += -0.25 * (A4_2 - B2_2) * Dv1[JJ]
            @inbounds O[pII,pII] += -0.25 * (B2_2 - A2_2) * Dv1[II]

            @inbounds O[pII,pII] += -0.25 * (A3_2 - B1_2) * Du[δx⁺(II)]
            @inbounds O[pII,pII] += -0.25 * (B1_2 - A1_2) * Du[II]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A4_2 - B2_2) * Dv1[JJ]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B2_2 - A2_2) * Dv1[II]

            @inbounds B[pII] += -0.25 * Dv2[II] * (A3_2 - B1_2) * Du[δx⁺(II)]
            @inbounds B[pII] += -0.25 * Dv2[II] * (B1_2 - A1_2) * Du[II]
        end
    end
    if is_periodic(BC.left.t) && is_periodic(BC.right.t)
        @inbounds for (II,JJ) in zip(b_left[2:end-1], b_right[2:end-1])
            pII = lexicographic(II, n+1)
            pJJ = lexicographic(JJ, n+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
            
            Auim1jm1 = A1_1 * u[δy⁻(II)]
            Auim1jp1 = A1_2 * u[II]
    
            Au1 = 0.5 * (Auim1jm1 + Auim1jp1)
    
            @inbounds O[pII,pII] += -0.5 * Au1
            @inbounds O[pII,pJJ] = -0.5 * Au1
        end
        @inbounds for (II,JJ) in zip(b_right[2:end-1], b_left[2:end-1])
            pII = lexicographic(II, n+1)
            pJJ = lexicographic(JJ, n+1)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(II), Δ)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, II, Δ)
            
            Auip1jm1 = A3_1 * u[δx⁺(δy⁻(II))]
            Auip1jp1 = A3_2 * u[δx⁺(II)]
    
            Au3 = 0.5 * (Auip1jm1 + Auip1jp1)
    
            @inbounds O[pII,pII] += 0.5 * Au3
            @inbounds O[pII,pJJ] = 0.5 * Au3
        end

        ii = b_bottom[1]
        pii = lexicographic(ii, n+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii, Δ)
        
        Auim1jp1 = A1_2 * u[ii]

        Au1 = 0.5 * Auim1jp1

        if is_periodic(BC.bottom.t)
            JJ = ii + CartesianIndex(n, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Auim1jm1 = A1_1 * u[JJ]
            Au1 += 0.5 * Auim1jm1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, n+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_bottom[end]
        pii = lexicographic(ii, n+1)
        A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, ii, Δ)
        
        Auip1jp1 = A3_2 * u[δx⁺(ii)]

        Au3 = 0.5 * Auip1jp1

        if is_periodic(BC.bottom.t)
            JJ = II + CartesianIndex(n, 0)
            A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, JJ, Δ)
            Auip1jm1 = A3_1 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jm1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, n+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
        
        ii = b_top[1]
        pii = lexicographic(ii, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii), Δ)
        
        Auim1jm1 = A1_1 * u[δy⁻(ii)]

        Au1 = 0.5 * Auim1jm1

        if is_periodic(BC.top.t)
            JJ = II + CartesianIndex(-n, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Auim1jp1 = A1_2 * u[JJ]
            Au1 += 0.5 * Auim1jp1
        end

        JJ = ii + CartesianIndex(0, n-1)
        pJJ = lexicographic(JJ, n+1)
        @inbounds O[pii,pii] += -0.5 * Au1
        @inbounds O[pii,pJJ] = -0.5 * Au1
        
        ii = b_top[end]
        pii = lexicographic(ii, n+1)
        A1_1, A2_1, A3_1, A4_1, B1_1, B2_1 = get_capacities_convection(cap, δy⁻(ii), Δ)
        
        Auip1jm1 = A3_1 * u[δx⁺(δy⁻(ii))]

        Au3 = 0.5 * Auip1jm1

        if is_periodic(BC.top.t)
            JJ = II + CartesianIndex(n, 0)
            A1_2, A2_2, A3_2, A4_2, B1_2, B2_2 = get_capacities_convection(cap, JJ, Δ)
            Auip1jp1 = A3_2 * u[δx⁺(JJ)]
            Au3 += 0.5 * Auip1jp1
        end

        JJ = ii + CartesianIndex(0, -n+1)
        pJJ = lexicographic(JJ, n+1)
        @inbounds O[pii,pii] += 0.5 * Au3
        @inbounds O[pii,pJJ] = 0.5 * Au3
    end

    # @inbounds _A1 = cap[:,:,1]
    # @inbounds _A2 = cap[:,:,2]
    # @inbounds _A3 = cap[:,:,3]
    # @inbounds _A4 = cap[:,:,4]
    # @inbounds _B1 = cap[:,:,6]
    # @inbounds _B2 = cap[:,:,7]

    # set_vec_conv_bnd!(dir, BC.left.t, O, x->x, δy⁺, x->x, δx⁺, x->x, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_x, Du_x, n+1, Δ, b_left[1:end-1], b_right[1:end-1])
    # set_vec_conv_bnd!(dir, BC.left.t, O, δy⁻, x->x, δy⁻, (δy⁻ ∘ δx⁺), δy⁻, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_x, Du_x, n+1, Δ, b_left[2:end], b_right[2:end])
    # set_vec_conv_bnd!(dir, BC.bottom.t, O, x->x, δy⁺, x->x, δx⁺, x->x, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_y, Du_y, n+1, Δ, b_bottom, b_top)
    # set_vec_conv_bnd!(dir, BC.right.t, O, x->x, δy⁺, x->x, δx⁺, x->x, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_x, Du_x, n+1, Δ, b_right[1:end-1], b_left[1:end-1])
    # set_vec_conv_bnd!(dir, BC.right.t, O, δy⁻, x->x, δy⁻, (δy⁻ ∘ δx⁺), δy⁻, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_x, Du_x, n+1, Δ, b_right[2:end], b_left[2:end])
    # set_vec_conv_bnd!(dir, BC.top.t, O, δy⁻, x->x, δy⁻, (δy⁻ ∘ δx⁺), δy⁻, _A4, _B2, _A2, _A3, _B1, _A1, Dv1_y, Du_y, n+1, Δ, b_top, b_bottom)

    return nothing
end
