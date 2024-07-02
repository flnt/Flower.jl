function indices_extension(grid, LS, inside_ext, periodic_x, periodic_y)
    _MIXED_L = collect(intersect(Set(LS.MIXED), Set(findall(LS.geoL.emptied))))
    _MIXED_S = collect(intersect(Set(LS.MIXED), Set(findall(LS.geoS.emptied))))
    _MIXED = vcat(_MIXED_L, _MIXED_S)
    indices_ext1 = vcat(LS.SOLID, _MIXED, LS.LIQUID)

    if periodic_x && periodic_y
        indices_ext = collect(intersect(Set(indices_ext1), Set(vcat(
            vec(inside_ext), grid.ind.b_left[1][2:end-1], grid.ind.b_bottom[1][2:end-1],
            grid.ind.b_right[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))))
    elseif !periodic_x && periodic_y
        indices_ext = collect(intersect(Set(indices_ext1), Set(vcat(
            vec(inside_ext), grid.ind.b_bottom[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))))
    elseif periodic_x && !periodic_y
        indices_ext = collect(intersect(Set(indices_ext1), Set(vcat(
            vec(inside_ext), grid.ind.b_left[1][2:end-1], grid.ind.b_right[1][2:end-1]
        ))))
    else
        indices_ext = collect(intersect(Set(indices_ext1), Set(inside_ext)))
    end

    left_ext = collect(intersect(Set(grid.ind.b_left[1][2:end-1]), Set(indices_ext1)))
    bottom_ext = collect(intersect(Set(grid.ind.b_bottom[1][2:end-1]), Set(indices_ext1)))
    right_ext = collect(intersect(Set(grid.ind.b_right[1][2:end-1]), Set(indices_ext1)))
    top_ext = collect(intersect(Set(grid.ind.b_top[1][2:end-1]), Set(indices_ext1)))

    return indices_ext, left_ext, bottom_ext, right_ext, top_ext
end

function update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, empty = true)
    if num.nLS > 1
        for iLS in 1:num.nLS
            update_ls_data(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].u, grid.LS[iLS].κ, BC_int, BC_int[iLS], periodic_x, periodic_y, false, empty)
            grid.LS[iLS].geoL.cap0 .= grid.LS[iLS].geoL.cap
            grid.LS[iLS].mid_point0 .= grid.LS[iLS].mid_point
        end
        combine_levelsets!(num, grid)
        NB_indices = update_ls_data(num, grid, grid_u, grid_v, num._nLS, grid.LS[end].u, grid.LS[end].κ, BC_int, DummyBC(), periodic_x, periodic_y, true, empty)
        grid.LS[end].geoL.cap0 .= grid.LS[end].geoL.cap
        grid.LS[end].mid_point0 .= grid.LS[end].mid_point
        crossing_2levelsets!(num, grid, grid.LS[1], grid.LS[2], BC_int)
        crossing_2levelsets!(num, grid_u, grid_u.LS[1], grid_u.LS[2], BC_int)
        crossing_2levelsets!(num, grid_v, grid_v.LS[1], grid_v.LS[2], BC_int)
        # for iLS in 1:num.nLS
        #     extend_contact_line!(grid, grid.LS[iLS])
        # end

        for iLS in 1:num.nLS
            postprocess_grids2!(grid, grid.LS[iLS], grid_u, grid_u.LS[iLS], grid_v, grid_v.LS[iLS], periodic_x, periodic_y, false)
        end
        postprocess_grids2!(grid, grid.LS[end], grid_u, grid_u.LS[end], grid_v, grid_v.LS[end], periodic_x, periodic_y, true)
    else
        NB_indices = update_ls_data(num, grid, grid_u, grid_v, 1, grid.LS[1].u, grid.LS[1].κ, BC_int, BC_int[1], periodic_x, periodic_y, true, empty)
        grid.LS[1].geoL.cap0 .= grid.LS[1].geoL.cap
        grid.LS[1].mid_point0 .= grid.LS[1].mid_point

        printstyled(color=:red, @sprintf "\n volume:\n")
        II = CartesianIndex(68, 1) #(id_y, id_x)
        print(grid.LS[1].geoL.dcap[II,8:11])
        print(grid.LS[1].geoL.cap[II,8:11])


        postprocess_grids2!(grid, grid.LS[1], grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y, true)
    end

    return NB_indices
end

function update_ls_data(num, grid, grid_u, grid_v, iLS, u, κ, BC_int, bc_int, periodic_x, periodic_y, neighbours, empty = true)
    NB_indices = update_ls_data_grid(num, grid, grid.LS[iLS], u, κ, periodic_x, periodic_y)

    interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.LS[iLS].u, grid_v.LS[iLS].u)

    _ = update_ls_data_grid(num, grid_u, grid_u.LS[iLS], grid_u.LS[iLS].u, grid_u.LS[iLS].κ, periodic_x, periodic_y)
    _ = update_ls_data_grid(num, grid_v, grid_v.LS[iLS], grid_v.LS[iLS].u, grid_v.LS[iLS].κ, periodic_x, periodic_y)

    postprocess_grids1!(num, grid, grid.LS[iLS], grid_u, grid_u.LS[iLS], grid_v, grid_v.LS[iLS], periodic_x, periodic_y, neighbours, empty, BC_int)

    for i in 1:num.nLS
        if is_wall(BC_int[i])
            idx_solid = Base.union(grid.LS[i].SOLID, findall(grid.LS[i].geoL.emptied))
            # @inbounds κ[idx_solid] .= 0.0
            @inbounds κ[grid.LS[i].SOLID] .= 0.0
        end
    end

    i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
    field_extension!(grid, u, κ, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

    if is_fs(bc_int)
        locate_contact_line!(num, grid, iLS, grid.LS[iLS].cl, grid.LS[iLS].MIXED, BC_int)
        locate_contact_line!(num, grid_u, iLS, grid_u.LS[iLS].cl, grid_u.LS[iLS].MIXED, BC_int)
        locate_contact_line!(num, grid_v, iLS, grid_v.LS[iLS].cl, grid_v.LS[iLS].MIXED, BC_int)

        # Apply inhomogeneous BC to more than one grid point by extending the contact line
        # extend_contact_line!(grid, grid.LS[iLS])
        # extend_contact_line!(grid_u, grid_u.LS[iLS].cl, num.n_ext_cl)
        # extend_contact_line!(grid_v, grid_v.LS[iLS].cl, num.n_ext_cl)
    end

    return NB_indices
end

function update_ls_data_grid(num, grid, LS, u, κ, periodic_x, periodic_y)
    LS.α .= NaN
    LS.faces .= 0.0
    LS.mid_point .= [Point(0.0, 0.0)]

    marching_squares!(grid, LS, u, periodic_x, periodic_y)

    LS.MIXED, LS.SOLID, LS.LIQUID = get_cells_indices(LS.iso, grid.ind.all_indices)
    NB_indices_base = get_NB_width_indices_base(num.NB)
    NB_indices = get_NB_width(grid, LS.MIXED, NB_indices_base)

    get_interface_location!(grid, LS, periodic_x, periodic_y)

    LS.geoL.emptied .= false
    LS.geoS.emptied .= false

    LS.geoL.double_emptied .= false
    LS.geoS.double_emptied .= false

    κ .= 0.0
    get_curvature(num, grid, LS.geoL, u, κ, LS.MIXED, periodic_x, periodic_y)

    return NB_indices
end

function update_stefan_velocity(num, grid, iLS, u, TS, TL, periodic_x, periodic_y, λ, Vmean)
    Stefan_velocity!(num, grid, grid.LS[iLS], grid.V, TS, TL, grid.LS[iLS].MIXED, periodic_x, periodic_y)
    grid.V[grid.LS[iLS].MIXED] .*= 1. ./ λ
    if Vmean
        a = mean(grid.V[grid.LS[iLS].MIXED])
        grid.V[grid.LS[iLS].MIXED] .= a
    end

    i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
    field_extension!(grid, u, grid.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)
end

function update_free_surface_velocity(num, grid_u, grid_v, iLS, uD, vD, periodic_x, periodic_y)
    # grid_u.V .= reshape(veci(uD,grid_u,iLS+1), grid_u)
    # grid_v.V .= reshape(veci(vD,grid_v,iLS+1), grid_v)
    grid_u.V .= reshape(vec1(uD,grid_u), grid_u)
    grid_v.V .= reshape(vec1(vD,grid_v), grid_v)

    # i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext = indices_extension(grid_u, grid_u.LS[iLS], grid_u.ind.inside, periodic_x, periodic_y)
    # i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext = indices_extension(grid_v, grid_v.LS[iLS], grid_v.ind.inside, periodic_x, periodic_y)

    # field_extension!(grid_u, grid_u.LS[iLS].u, grid_u.V, i_u_ext, l_u_ext, b_u_ext, r_u_ext, t_u_ext, num.NB, periodic_x, periodic_y)
    # field_extension!(grid_v, grid_v.LS[iLS].u, grid_v.V, i_v_ext, l_v_ext, b_v_ext, r_v_ext, t_v_ext, num.NB, periodic_x, periodic_y)
end