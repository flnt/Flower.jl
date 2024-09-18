"""
    BC_LS!(num, cl, grid, A, B, rhs, BC)

Update levelset matrices to apply inhomogeneous Neumann boundary conditions in presence of
contact lines. 

Outside, the contact angle asymptotically converges to an angle of 90°. Inside, the contact
angle converges to an angle of 0° if the imposed contact angle at the contact line is
smaller than 90° and to an angle of 180° if the imposed contact angle is bigger than 90°.
"""
function BC_LS!(grid, u, A, B, rhs, BC)
    @unpack x, y, nx, ny, dx, dy, ind = grid
    @unpack all_indices, b_left, b_bottom, b_right, b_top = ind
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
    npts = [ny, nx, ny, nx]

    for (i, (idx, idx2, idx3, xy, nxy)) in enumerate(zip(boundaries_idx, boundaries2, boundaries3, direction, npts))
        pks, _ = findminima(abs.(u[idx]))
        if is_neumann(boundaries_t[i])
            for (II, JJ) in zip(idx, idx2)
                pII = lexicographic(II, grid.ny)
                pJJ = lexicographic(JJ, grid.ny)

                A[pII,:] .= 0.0
                A[pII,pII] = 1.0
                A[pII,pJJ] = -1.0
                B[pII,:] .= 0.0
            end
        elseif is_neumann_cl(boundaries_t[i]) && maximum(u[idx]) > 0.0 && minimum(u[idx]) < 0.0 && length(pks) >= 2
            pks1 = idx[pks[1]]
            pkse = idx[pks[end]]

            # Gradually update the contact angle
            Δθe = 2.0 * π / 180

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
                B[pII,:] .= 0.0

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
                B[pII,:] .= 0.0

                rhs[pII] = u[JJ] - u[KK]
            end
        end
    end

    return nothing
end

"""
    BC_LS_interior!(num, grid, iLS, A, B, rhs, BC_int, periodic_x, periodic_y)

Update levelset matrices to apply inhomogeneous Neumann boundary conditions in presence of
contact lines at the intersection of the interfaces, one with a WallNoSLip() or Navier_cl() BC 
and the other one with a FreeSurface() BC. 

Outside, the contact angle asymptotically converges to an angle of 90°. Inside, the contact
angle converges to an angle of 0° if the imposed contact angle at the contact line is
smaller than 90° and to an angle of 180° if the imposed contact angle is bigger than 90°.
"""
function BC_LS_interior!(num, grid, grid_u, grid_v, iLS, A, B, rhs, BC_int, periodic_x, periodic_y)
    @unpack x, y, nx, ny, dx, dy, ind, LS = grid

    π2 = π / 2.0
    π4 = π / 4.0

    # Distance at which a 90 deg contact angle is imposed   
    Δu = 0.5

    for i in 1:num.nLS
        if is_wall(BC_int[i])
            idx = LS[i].MIXED
            idx_mix_full = vcat(LS[i].LIQUID)

            if !isempty(intersect(LS[iLS].MIXED, idx))
                pks, _ = findminima(abs.(LS[iLS].u[idx]))
                base = get_NB_width_indices_base1(2)
                idx_ext = get_NB_width(grid, idx[pks], base)

                pks1 = idx_ext[findmin(grid.x[idx_ext])[2]]
                pkse = idx_ext[findmax(grid.x[idx_ext])[2]]

                θtmp = zeros(grid)
                JJtmp = zeros(CartesianIndex{2}, ny, nx)
                multtmp = zeros(grid)

                for II in idx
                    neighbours = Vector{CartesianIndex{2}}()
                    s = Vector{Float64}()
                    oppo = Vector{CartesianIndex{2}}()
                    perp = Vector{CartesianIndex{2}}()
                    mult = Vector{Float64}()
                    if (II[2] > 1 || periodic_x) && δx⁻(II, nx, periodic_x) in idx_mix_full
                        push!(neighbours, δx⁻(II, nx, periodic_x))
                        push!(s, π)
                        if (II[1] < ny || periodic_y)
                            push!(perp, δy⁺(II, ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁻(II, nx, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δx⁺(II, nx, periodic_x))
                    end
                    if (II[2] > 1 || periodic_x) && δx⁻(δx⁻(II, nx, periodic_x), nx, periodic_x) in idx_mix_full
                        push!(neighbours, δx⁻(δx⁻(II, nx, periodic_x), nx, periodic_x))
                        push!(s, π)
                        if (II[1] < ny || periodic_y)
                            push!(perp, δy⁺(II, ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁻(II, nx, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δx⁺(II, nx, periodic_x))
                    end
                    if (II[1] > 1 || periodic_y) && δy⁻(II, ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁻(II, ny, periodic_y))
                        push!(s, π + π2)
                        if (II[2] > 1 || periodic_x)
                            push!(perp, δx⁻(II, nx, periodic_x))
                            push!(mult, 1.0)
                        else
                            push!(perp, δx⁺(II, nx, periodic_x))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁺(II, ny, periodic_y))
                    end
                    if (II[1] > 1 || periodic_y) && δy⁻(δy⁻(II, ny, periodic_y), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁻(δy⁻(II, ny, periodic_y), ny, periodic_y))
                        push!(s, π + π2)
                        if (II[2] > 1 || periodic_x)
                            push!(perp, δx⁻(II, nx, periodic_x))
                            push!(mult, 1.0)
                        else
                            push!(perp, δx⁺(II, nx, periodic_x))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁺(II, ny, periodic_y))
                    end
                    if (II[2] < nx || periodic_x) && δx⁺(II, nx, periodic_x) in idx_mix_full
                        push!(neighbours, δx⁺(II, nx, periodic_x))
                        push!(s, 0.0)
                        if (II[1] > 1 || periodic_y)
                            push!(perp, δy⁻(II, ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁺(II, nx, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δx⁻(II, nx, periodic_x))
                    end
                    if (II[2] < nx || periodic_x) && δx⁺(δx⁺(II, nx, periodic_x), nx, periodic_x) in idx_mix_full
                        push!(neighbours, δx⁺(δx⁺(II, nx, periodic_x), nx, periodic_x))
                        push!(s, 0.0)
                        if (II[1] > 1 || periodic_y)
                            push!(perp, δy⁻(II, ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁺(II, nx, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δx⁻(II, nx, periodic_x))
                    end
                    if (II[1] < ny || periodic_y) && δy⁺(II, ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁺(II, ny, periodic_y))
                        push!(s, π2)
                        if (II[2] < nx || periodic_x)
                            push!(perp, δx⁺(II, nx, periodic_x))
                            push!(mult, 1.0)
                        else
                            push!(perp, δx⁻(II, nx, periodic_x))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁻(II, ny, periodic_y))
                    end
                    if (II[1] < ny || periodic_y) && δy⁺(δy⁺(II, ny, periodic_y), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁺(δy⁺(II, ny, periodic_y), ny, periodic_y))
                        push!(s, π2)
                        if (II[2] < nx || periodic_x)
                            push!(perp, δx⁺(II, nx, periodic_x))
                            push!(mult, 1.0)
                        else
                            push!(perp, δx⁻(II, nx, periodic_x))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁻(II, ny, periodic_y))
                    end
                    if ((II[2] > 1 || periodic_x) && (II[1] > 1 || periodic_y)) && δy⁻(δx⁻(II, nx, periodic_x), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁻(δx⁻(II, nx, periodic_x), ny, periodic_y))
                        push!(s, π + π4)
                        if ((II[2] > 1 || periodic_x) && (II[1] < ny || periodic_y))
                            push!(perp, δy⁺(δx⁻(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁻(δx⁺(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁺(δx⁺(II, nx, periodic_x), ny, periodic_y))
                    end
                    if ((II[2] < nx || periodic_x) && (II[1] > 1 || periodic_y)) && δy⁻(δx⁺(II, nx, periodic_x), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁻(δx⁺(II, nx, periodic_x), ny, periodic_y))
                        push!(s, π + π2 + π4)
                        if ((II[2] > 1 || periodic_x) && (II[1] > 1 || periodic_y))
                            push!(perp, δy⁻(δx⁻(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁺(δx⁺(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁺(δx⁻(II, nx, periodic_x), ny, periodic_y))
                    end
                    if ((II[2] < nx || periodic_x) && (II[1] < ny || periodic_y)) && δy⁺(δx⁺(II, nx, periodic_x), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁺(δx⁺(II, nx, periodic_x), ny, periodic_y))
                        push!(s, π4)
                        if ((II[2] < nx || periodic_x) && (II[1] > 1 || periodic_y))
                            push!(perp, δy⁻(δx⁺(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁺(δx⁻(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁻(δx⁻(II, nx, periodic_x), ny, periodic_y))
                    end
                    if ((II[2] > 1 || periodic_x) && (II[1] < ny || periodic_y)) && δy⁺(δx⁻(II, nx, periodic_x), ny, periodic_y) in idx_mix_full
                        push!(neighbours, δy⁺(δx⁻(II, nx, periodic_x), ny, periodic_y))
                        push!(s, π2 + π4)
                        if ((II[2] < nx || periodic_x) && (II[1] < ny || periodic_y))
                            push!(perp, δy⁺(δx⁺(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, 1.0)
                        else
                            push!(perp, δy⁻(δx⁻(II, nx, periodic_x), ny, periodic_y))
                            push!(mult, -1.0)
                        end
                        push!(oppo, δy⁻(δx⁺(II, nx, periodic_x), ny, periodic_y))
                    end

                    dist = zeros(size(neighbours))
                    for (j,neigh) in enumerate(neighbours)
                        dist[j] = distance(x[neigh], y[neigh], LS[i].α[II], x[II], y[II])
                    end
                    idmin = findmin(dist)[2]
                    
                    JJ = neighbours[idmin]
                    JJtmp[II] = JJ

                    KK = perp[idmin]
                    multtmp[II] = mult[idmin] * sign(LS[iLS].u[KK] - LS[iLS].u[II])

                    d2 = sqrt((x[pks1] - x[pkse])^2.0 + (y[pks1] - y[pkse])^2.0) / 2.0
                    
                    #Imposing contact angle for wall BC
                    printstyled(color=:magenta, @sprintf "\n Imposing contact angle: %.2e" BC_int[i].θe)
                    
                    if num.contact_angle == 1
                        advancing_receding = contact_angle_advancing_receding(grid_u, grid_v, i, II)
                        print("\n advancing/receding II ",II, " ",advancing_receding)
                        if advancing_receding>0
                            theta_adv_reced = BC_int[i].θadv
                        else
                            theta_adv_reced = BC_int[i].θrec
                        end

                        if theta_adv_reced < π2 
                            newθ = atan(tan(theta_adv_reced) * (1.0 - LS[iLS].u[II] / d2))
                        else
                            newθ = π - atan(tan(π - theta_adv_reced) * (1.0 - LS[iLS].u[II] / d2))
                        end
    
                        θtmp[II] = theta_adv_reced
                    else
                        if BC_int[i].θe < π2 
                            newθ = atan(tan(BC_int[i].θe) * (1.0 - LS[iLS].u[II] / d2))
                        else
                            newθ = π - atan(tan(π - BC_int[i].θe) * (1.0 - LS[iLS].u[II] / d2))
                        end
    
                        θtmp[II] = BC_int[i].θe
                    end


                end

                βtmp = zeros(grid)
                cl = intersect(idx, LS[iLS].MIXED)

                for II in idx
                    pII = lexicographic(II, ny)
                    pointII = Point(x[II], y[II])

                    βv = Float64[]
                    KKv = CartesianIndex{2}[]

                    ## to be changed for the normal drops!!!
                    tmpmult = 1.0

                    # if II[2] < nx || periodic_x
                    #     KK = δx⁺(II, nx, periodic_x)
                    #     push!(KKv, KK)
                    #     γ = atan(y[KK] - y[II], x[KK] - x[II])
                    #     # β = θtmp[II] + multtmp[II] * (LS[i].α[II] - γ)
                    #     if x[KK] > 0.0
                    #         tmpmult = 1.0
                    #     end
                    #     β = θtmp[II] + tmpmult * (LS[i].α[II] - γ)
                    #     push!(βv, β)
                    # end
                    if II[1] < ny || periodic_y
                        KK = δy⁺(II, ny, periodic_y)
                        push!(KKv, KK)
                        γ = atan(y[KK] - y[II], x[KK] - x[II])
                        # β = θtmp[II] + multtmp[II] * (LS[i].α[II] - γ)
                        if x[KK] > 0.0
                            tmpmult = -1.0
                        end
                        β = θtmp[II] + tmpmult * (LS[i].α[II] - γ)
                        push!(βv, β)
                    end
                    if II[1] > 1 || periodic_y
                        KK = δy⁻(II, ny, periodic_y)
                        push!(KKv, KK)
                        γ = atan(y[KK] - y[II], x[KK] - x[II])
                        # β = θtmp[II] + multtmp[II] * (LS[i].α[II] - γ)
                        if x[KK] > 0.0
                            tmpmult = -1.0
                        end
                        β = θtmp[II] + tmpmult * (LS[i].α[II] - γ)
                        push!(βv, β)
                    end
                    # if II[2] > 1 || periodic_x
                    #     KK = δx⁻(II, nx, periodic_x)
                    #     push!(KKv, KK)
                    #     γ = atan(y[KK] - y[II], x[KK] - x[II])
                    #     # β = θtmp[II] + multtmp[II] * (LS[i].α[II] - γ)
                    #     if x[KK] > 0.0
                    #         tmpmult = 1.0
                    #     end
                    #     β = θtmp[II] + tmpmult * (LS[i].α[II] - γ)
                    #     push!(βv, β)
                    # end

                    # Ensure that all angles are positive
                    for j in eachindex(βv)
                        if βv[j] < 0.0
                            βv[j] += 2π
                        end
                    end

                    αtmp = LS[i].α[II]
                    if LS[i].α[II] < 0.0
                        αtmp += 2π
                    end
                    # y-direction
                    if αtmp >= 0.0 && αtmp < π
                        KK = KKv[1]
                        β = βv[1]
                    else
                        KK = KKv[2]
                        β = βv[2]
                    end
                    # x-direction
                    # if αtmp >= π / 2.0 && αtmp < 3π / 2.0
                    #     KK = KKv[2]
                    #     β = βv[2]
                    # else
                    #     KK = KKv[1]
                    #     β = βv[1]
                    # end
                    # xy-direction
                    # if αtmp >= π / 4.0 && αtmp < 3π / 4.0
                    #     KK = KKv[2]
                    #     β = βv[2]
                    # elseif αtmp >= 3π / 4.0 && αtmp < 5π / 4.0
                    #     KK = KKv[4]
                    #     β = βv[4]
                    # elseif αtmp >= 5π / 4.0 && αtmp < 7π / 4.0
                    #     KK = KKv[3]
                    #     β = βv[3]
                    # else
                    #     KK = KKv[1]
                    #     β = βv[1]
                    # end

                    pKK = lexicographic(KK, ny)
                    pointKK = Point(x[KK], y[KK])
                    d = distance(pointII, pointKK)

                    # Ensure that the angle doesn't change too much
                    Δβ = 1.0 * π / 180
                    old = LS[iLS].u[II] - LS[iLS].u[KK]
                    # Levelset difference between two consecutive points might be bigger
                    # than the distance between them if it's not reinitialized often enough
                    if abs(old) > d
                        old = sign(old) * d
                    end
                    βold = acos(old / d)

                    if abs(β - 2π - βold) < abs(β - βold)
                        β -= 2π
                    end
                    # Compute new contact angle
                    if abs(β - βold) > Δβ
                        βnew = βold + sign(β - βold) * Δβ
                        θtmp[II] -= β - βnew
                    else
                        βnew = β
                    end

                    dist = 1e10
                    id = CartesianIndex(1,1)
                    for JJ in cl
                        pointJJ = Point(x[JJ], y[JJ])
                        _d = distance(pointJJ, pointII)
                        if _d < dist
                            dist = _d
                            id = JJ
                        end
                    end

                    A[pII,:] .= 0.0
                    A[pII,pII] = 1.0
                    A[pII,pKK] = -1.0
                    B[pII,:] .= 0.0

                    if abs(LS[iLS].u[II]) < Δu
                        rhs[pII] = d * cos(βnew) * bell_function2(LS[iLS].u[II], Δu)
                    else
                        rhs[pII] = 0.0
                    end

                    βtmp[II] = βnew
                end

                @inbounds @threads for II in LS[i].SOLID
                    pII = lexicographic(II, ny)
                    pointII = Point(x[II], y[II])

                    # Find closest interfacial cell
                    dist_mixed = zeros(size(idx))
                    for (k, mix) in enumerate(idx)
                        pointk = Point(x[mix], y[mix])
                        pointkd = Point(dx[mix], dy[mix])
                        midpoint = pointk + pointkd * LS[i].mid_point0[mix]
                        dist_mixed[k] = distance(midpoint, pointII)
                    end
                    idα = idx[findmin(dist_mixed)[2]]

                    xm = ind.all_indices[II[1],1:II[2]-1]
                    ym = ind.all_indices[1:II[1]-1,II[2]]
                    xp = ind.all_indices[II[1],II[2]+1:end]
                    yp = ind.all_indices[II[1]+1:end,II[2]]

                    mixed_xm = intersect(xm, idx)
                    mixed_ym = intersect(ym, idx)
                    mixed_xp = intersect(xp, idx)
                    mixed_yp = intersect(yp, idx)

                    dist_yp = 1e10
                    id_yp = CartesianIndex(1,1)
                    for JJ in mixed_yp
                        pointJJ = Point(x[JJ], y[JJ])
                        pointJJd = Point(dx[JJ], dy[JJ])
                        midpoint = pointJJ + pointJJd * LS[i].mid_point0[JJ]
                        d = distance(midpoint, pointII)
                        if d < dist_yp
                            dist_yp = d
                            id_yp = JJ
                        end
                    end

                    dist_ym = 1e10
                    id_ym = CartesianIndex(1,1)
                    for JJ in mixed_ym
                        pointJJ = Point(x[JJ], y[JJ])
                        pointJJd = Point(dx[JJ], dy[JJ])
                        midpoint = pointJJ + pointJJd * LS[i].mid_point0[JJ]
                        d = distance(midpoint, pointII)
                        if d < dist_ym
                            dist_ym = d
                            id_ym = JJ
                        end
                    end

                    dist_xp = 1e10
                    id_xp = CartesianIndex(1,1)
                    for JJ in mixed_xp
                        pointJJ = Point(x[JJ], y[JJ])
                        pointJJd = Point(dx[JJ], dy[JJ])
                        midpoint = pointJJ + pointJJd * LS[i].mid_point0[JJ]
                        d = distance(midpoint, pointII)
                        if d < dist_xp
                            dist_xp = d
                            id_xp = JJ
                        end
                    end

                    dist_xm = 1e10
                    id_xm = CartesianIndex(1,1)
                    for JJ in mixed_xm
                        pointJJ = Point(x[JJ], y[JJ])
                        pointJJd = Point(dx[JJ], dy[JJ])
                        midpoint = pointJJ + pointJJd * LS[i].mid_point0[JJ]
                        d = distance(midpoint, pointII)
                        if d < dist_xm
                            dist_xm = d
                            id_xm = JJ
                        end
                    end

                    if dist_yp < dist_ym
                        KK = id_yp
                    else
                        KK = id_ym
                    end
                    # if dist_xp < dist_xm
                    #     KK = id_xp
                    # else
                    #     KK = id_xm
                    # end
                    # xy-direction
                    # αtmp = LS[i].α[idα]
                    # if LS[i].α[idα] < 0.0
                    #     αtmp += 2π
                    # end
                    # if αtmp >= π / 4.0 && αtmp < 3π / 4.0
                    #     KK = id_yp
                    # elseif αtmp >= 3π / 4.0 && αtmp < 5π / 4.0
                    #     KK = id_xm
                    # elseif αtmp >= 5π / 4.0 && αtmp < 7π / 4.0
                    #     KK = id_ym
                    # else
                    #     KK = id_xp
                    # end

                    tmpLL1 = KK - II
                    tmpLL2 = CartesianIndex(sign(tmpLL1[1]), sign(tmpLL1[2]))
                    LL = II + tmpLL2
                    pLL = lexicographic(LL, ny)
                    pointLL = Point(x[LL], y[LL])

                    ## to be changed for the normal drops!!!
                    # tmpmult = 1.0
                    # if x[KK] > 0.0
                    #     tmpmult = -1.0
                    # end

                    # γ = atan(y[LL] - y[II], x[LL] - x[II])
                    # β = θtmp[KK] + multtmp[KK] * (LS[i].α[KK] - γ)
                    # β = θtmp[KK] + tmpmult * (LS[i].α[KK] - γ)

                    β = βtmp[KK]

                    # if β < π2
                    #     β = atan(tan(β) * (1.0 - LS[i].u[II] / (0.5 * abs(minimum(LS[i].u)))))
                    # end

                    d = distance(pointII, pointLL)

                    A[pII,:] .= 0.0
                    A[pII,pII] = 1.0
                    A[pII,pLL] = -1.0
                    B[pII,:] .= 0.0

                    if abs(LS[iLS].u[II]) < Δu
                        rhs[pII] = d * cos(β) * bell_function2(LS[iLS].u[II], Δu)
                    else
                        rhs[pII] = 0.0
                    end
                end
            end
        end
    end

    return nothing
end


"""
    locate_contact_line!(grid)

Intersection between mixed cells and borders of the domain.

The contact lines inside the domain (intersection of two levelsets) are
computed inside `crossing_2levelsets!` in `cutcell.jl`.
"""
function locate_contact_line!(num, grid, iLS, cl, MIXED, BC_int)
    @unpack LS, ind = grid

    intersect!(cl, [])
    # Contact lines at domain borders
    boundaries_idx = [ind.b_left, ind.b_bottom, ind.b_right, ind.b_top]
    for idx in boundaries_idx
        append!(cl, intersect(MIXED, idx[1]))
    end

    return nothing
end

"""
    extend_contact_line(grid, n)

Extend every contact line point by `n` points parallel to the boundary so that the
inhomogeneous Neumann boundary conditions are also applied on them.
"""
function extend_contact_line!(grid, cl, n_ext)
    @unpack nx, ny, ind = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

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

function extend_contact_line!(grid, LS)
    base = get_NB_width_indices_base(3)
    LS.cl = Base.union(LS.cl, get_NB_width(grid, LS.cl, base))

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

Compute dynamic contact angle to be applied when using the generalized Navier BC.
"""
function dynamic_contact_angle(grid)
    @unpack ind, κ = grid
    @unpack b_left, b_bottom, b_right, b_top = ind

    θapp = zeros(grid)

    boundaries = [b_left, b_bottom, b_right, b_top]

    θd = θapp .+ 1.5 .* κ .+ sqrt(1 + hy) ./ sin(θapp)
   
    return θd
end