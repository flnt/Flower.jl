"""
    BC_LS!(num, grid, A, B, rhs, BC)

Update levelset matrices to apply inhomogeneous Neumann boundary conditions in presence of
contact lines.
"""
function BC_LS!(grid, A, B, rhs, BC, n_ext)
    @unpack x, y, nx, ny, dx, dy, ind, u = grid
    @unpack b_left, b_bottom, b_right, b_top, cl = ind
    @unpack left, bottom, right, top = BC

    boundaries_idx = [b_left, b_bottom, b_right, b_top]
    boundaries_per = [b_right, b_top, b_left, b_bottom]
    boundaries_t = [left, bottom, right, top]

    for (i, (idx, per)) in enumerate(zip(boundaries_idx, boundaries_per))
        if is_neumann(boundaries_t[i])
            for (II, JJ) in zip(idx[1], idx[2])
                pII = lexicographic(II, grid.ny)
                pJJ = lexicographic(JJ, grid.ny)

                dist = sqrt((x[JJ] - x[II])^2 + (y[JJ] - y[II])^2)
                ϵb = n_ext * dist

                A[pII,:] .= 0.0
                A[pII,pII] = 1.0
                A[pII,pJJ] = -1.0
                B[pII,pII] = 0.0
                if II in cl && u[II] < ϵb
                    # Gradually change the contact angle
                    Δθe = 5 * π / 180
                    old = u[II] - u[JJ]
                    if abs(old) > dist
                        old = sign(old) * dist
                    end
                    θe_old = acos(old / dist)
                    new = boundaries_t[i].val * bell_function2(u[II], ϵb)
                    if abs(new) > dist
                        new = sign(new) * dist
                    end
                    θe_new = acos(new / dist)
                    
                    if abs(θe_new - θe_old) > Δθe
                        new = dist * cos(θe_old + sign(old - new) * Δθe)
                    end

                    rhs[pII] = new
                end
            end
        end
    end

    return nothing
end


"""
    locate_contact_line!(grid)

Intersection between mixed cells and borders of the domain.

At some point it should include the intersection of two levelsets, for contact lines with
non-planar shapes.
"""
function locate_contact_line!(grid)
    @unpack ind = grid

    boundaries_idx = [ind.b_left, ind.b_bottom, ind.b_right, ind.b_top]
    ind.cl = []
    for idx in boundaries_idx
        append!(ind.cl, intersect(ind.MIXED, idx[1]))
    end

    return nothing
end

"""
    locate_contact_line!(num, grid, grid_u, grid_v, u, per_x, per_y)

Intersection between mixed cells and borders of the domain.

At some point it should include the intersection of two levelsets, for contact lines with
non-planar shapes.
"""
function locate_contact_line!(num, grid, grid_u, grid_v, u, per_x, per_y)
    @unpack ind = grid

    update_ls_data(num, grid, grid_u, grid_v, u, grid.κ, per_x, per_y)

    boundaries_idx = [ind.b_left, ind.b_bottom, ind.b_right, ind.b_top]
    ind.cl = []
    for idx in boundaries_idx
        append!(ind.cl, intersect(ind.MIXED, idx[1]))
    end

    return nothing
end

"""
    extend_contact_line(grid, n)

Extend every contact line point by `n` points parallel to the boundary so that the
inhomogeneous Neumann boundary conditions are also applied on them.
"""
function extend_contact_line!(grid, n_ext)
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top, cl = ind

    boundaries = [b_left, b_bottom, b_right, b_top]
    a = [1, 2, 1, 2]
    b = [ny, nx, ny, nx]
    c = [δy⁻, δx⁻, δy⁻, δx⁻]
    d = [δy⁺, δx⁺, δy⁺, δx⁺]

    _cl = []
    for (boundary, ai, bi, ci, di) in zip(boundaries, a, b, c, d)
        for II in cl
            if II in boundary[1]
                JJm = II
                JJp = II
                for j in 1:n_ext
                    if II[ai] >= 1 + j
                        JJm = ci(JJm)
                        push!(_cl, JJm)
                    end
                    if II[ai] <= bi - j
                        JJp = di(JJp)
                        push!(_cl, JJp)
                    end
                end
            end
        end
    end
    append!(cl, _cl)

    return nothing
end


"""
    bell_function(grid)

Define a bell function to be applied at contact lines when using the generalized Navier BC.
"""
bell_function(grid, ϵ) = (1.0 .- tanh(grid.u ./ ϵ).^2) ./ ϵ

"""
    bell_function2(grid)

Define a bell function to be applied at contact lines with the inhomogeneous Neumann BC.
"""
bell_function2(u, ϵ) = ((1.0 + cos(π * u / ϵ)) / 2.0)^2

"""
    dynamic_contact_angle(grid)

Compute dynamic contact angle to be applied when usung the generalized Navier BC.
"""
function dynamic_contact_angle(grid)
    @unpack ind, κ = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    θapp = zeros(grid)

    boundaries = [b_left, b_bottom, b_right, b_top]

    θd = θapp .+ 1.5 .* κ .+ sqrt(1 + hy) ./ sin(θapp)
   
    return θd
end