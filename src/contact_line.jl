"""
    BC_LS!(num, grid, A, B, rhs, BC)

Update levelset matrices to apply inhomogeneous Neumann boundary conditions in presence of
contact lines. 

Outside, the contact angle asimptotically converges to an angle of 90°. Inside, the contact
angle converges to an angle of 0° if the imposed contact angle at the contact line is
smaller than 90° and to an angle of 180° if the imposed contact angle is bigger than 90°.
"""
function BC_LS!(grid, A, B, rhs, BC, n_ext)
    @unpack x, y, nx, ny, dx, dy, ind, u = grid
    @unpack all_indices, b_left, b_bottom, b_right, b_top, cl = ind
    @unpack left, bottom, right, top = BC

    π2 = π / 2.0

    boundaries_idx = [b_left[1], b_bottom[1], b_right[1], b_top[1]]

    left2 = vcat(all_indices[2,2], all_indices[2:end-1,2], all_indices[end-1,2])
    bottom2 = vcat(all_indices[2,2], all_indices[2,2:end-1], all_indices[2,end-1])
    right2 = vcat(all_indices[2,end-1], all_indices[2:end-1,end-1], all_indices[end-1,end-1])
    top2 = vcat(all_indices[end-1,2], all_indices[end-1,2:end-1], all_indices[end-1,end-1])
    boundaries2 = [left2, bottom2, right2, top2]

    left3 = vcat(all_indices[3,3], all_indices[2:end-1,3], all_indices[end-2,3])
    bottom3 = vcat(all_indices[3,3], all_indices[3,2:end-1], all_indices[3,end-2])
    right3 = vcat(all_indices[3,end-2], all_indices[2:end-1,end-2], all_indices[end-2,end-2])
    top3 = vcat(all_indices[end-2,3], all_indices[end-2,2:end-1], all_indices[end-2,end-2])
    boundaries3 = [left3, bottom3, right3, top3]

    boundaries_t = [left, bottom, right, top]

    direction = [y, x, y, x]

    for (i, (idx, idx2, idx3, xy)) in enumerate(zip(boundaries_idx, boundaries2, boundaries3, direction))
        pks, _ = findminima(abs.(u[idx]))
        if is_neumann(boundaries_t[i])
            for (II, JJ) in zip(idx, idx2)
                pII = lexicographic(II, grid.ny)
                pJJ = lexicographic(JJ, grid.ny)

                dist = sqrt((x[JJ] - x[II])^2 + (y[JJ] - y[II])^2)
                ϵb = n_ext * dist

                A[pII,:] .= 0.0
                A[pII,pII] = 1.0
                A[pII,pJJ] = -1.0
                B[pII,pII] = 0.0
            end
        elseif is_neumann_cl(boundaries_t[i]) && maximum(u[idx]) > 0.0 && minimum(u[idx]) < 0.0 && length(pks) >= 2
            pks1 = idx[pks[1]]
            pkse = idx[pks[end]]

            # Gradually update the contact angle
            Δθe = 90.0 * π / 180

            # Find current contact angle
            dist = sqrt((x[idx2[pks[1]]] - x[pks1])^2 + (y[idx2[pks[1]]] - y[pks1])^2)
            old = u[pks1] - u[idx2[pks[1]]]
            # Levelset difference between two consecutive points might be bigger
            # than the distance between them if it's not reinitialized often enough
            if abs(old) > dist
                old = sign(old) * dist
            end
            θe_old = acos(old / dist)

            # Compute new contact angle
            if abs(boundaries_t[i].θe - θe_old) > Δθe
                θe = θe_old + sign(boundaries_t[i].θe - θe_old) * Δθe
            else
                θe = boundaries_t[i].θe
            end

            # distance between the center of the drop and the contact line
            d = abs(xy[pks1] + u[pks1] - (xy[pkse] + u[pkse])) / 2.0

            for (II, JJ) in zip(idx[2:end-1], idx2[2:end-1])
                pII = lexicographic(II, grid.ny)
                pJJ = lexicographic(JJ, grid.ny)

                A[pII,:] .= 0.0
                A[pII,pII] = 1.0
                A[pII,pJJ] = -1.0
                B[pII,pII] = 0.0

                # Compute levelset angle at a distance u[II] from the contact line
                if θe < π2
                    newθ = atan(tan(θe) * (1.0 - u[II] / d))
                else
                    newθ = π - atan(tan(π - θe) * (1.0 - u[II] / d))
                end

                rhs[pII] = dist * cos(newθ)
            end
        elseif is_neumann_cl(boundaries_t[i]) || is_neumann_inh(boundaries_t[i])
            for (II, JJ, KK) in zip(idx, idx2, idx3)
                pII = lexicographic(II, grid.ny)
                pJJ = lexicographic(JJ, grid.ny)

                A[pII,:] .= 0.0
                A[pII,pII] = 1.0
                A[pII,pJJ] = -1.0
                B[pII,pII] = 0.0

                rhs[pII] = u[JJ] - u[KK]
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