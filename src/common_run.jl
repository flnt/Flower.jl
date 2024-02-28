function indices_extension(grid, LS, inside_ext, periodic_x, periodic_y)
    _MIXED_L = intersect(findall(LS.geoL.emptied), LS.MIXED)
    _MIXED_S = intersect(findall(LS.geoS.emptied), LS.MIXED)
    _MIXED = vcat(_MIXED_L, _MIXED_S)
    indices_ext1 = vcat(LS.SOLID, _MIXED, LS.LIQUID)

    if periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(inside_ext), grid.ind.b_left[1][2:end-1], grid.ind.b_bottom[1][2:end-1],
            grid.ind.b_right[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif !periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(inside_ext), grid.ind.b_bottom[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif periodic_x && !periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(inside_ext), grid.ind.b_left[1][2:end-1], grid.ind.b_right[1][2:end-1]
        ))
    else
        indices_ext = intersect(indices_ext1, vec(inside_ext))
    end

    left_ext = intersect(indices_ext1, grid.ind.b_left[1][2:end-1])
    bottom_ext = intersect(indices_ext1, grid.ind.b_bottom[1][2:end-1])
    right_ext = intersect(indices_ext1, grid.ind.b_right[1][2:end-1])
    top_ext = intersect(indices_ext1, grid.ind.b_top[1][2:end-1])

    return indices_ext, left_ext, bottom_ext, right_ext, top_ext
end

function update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, empty = true)
    if num.nLS > 1
        for iLS in 1:num.nLS
            update_ls_data(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].u, grid.LS[iLS].κ, BC_int, BC_int[iLS], periodic_x, periodic_y, false, empty)
        end
        combine_levelsets!(num, grid)
        NB_indices = update_ls_data(num, grid, grid_u, grid_v, num._nLS, grid.LS[end].u, grid.LS[end].κ, BC_int, DummyBC(), periodic_x, periodic_y, true, empty)
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
        postprocess_grids2!(grid, grid.LS[end], grid_u, grid_u.LS[end], grid_v, grid_v.LS[end], periodic_x, periodic_y, true)
    end

    return NB_indices
end

function update_ls_data(num, grid, grid_u, grid_v, iLS, u, κ, BC_int, bc_int, periodic_x, periodic_y, neighbours, empty = true)
    NB_indices = update_ls_data_grid(num, grid, grid.LS[iLS], u, κ, periodic_x, periodic_y)

    interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.LS[iLS].u, grid_v.LS[iLS].u)

    _ = update_ls_data_grid(num, grid_u, grid_u.LS[iLS], grid_u.LS[iLS].u, grid_u.LS[iLS].κ, periodic_x, periodic_y)
    _ = update_ls_data_grid(num, grid_v, grid_v.LS[iLS], grid_v.LS[iLS].u, grid_v.LS[iLS].κ, periodic_x, periodic_y)

    postprocess_grids1!(grid, grid.LS[iLS], grid_u, grid_u.LS[iLS], grid_v, grid_v.LS[iLS], periodic_x, periodic_y, num.ϵ, neighbours, empty)

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

function adjoint_projection_fs(num, grid, grid_u, grid_v,
    J_u, J_v, adj, adj_ph, ph,
    RlsFS_ucorr, RlsFS_vcorr,
    Au, Bu, Av, Bv, Aϕ,
    opC_p, opC_u, opC_v, BC_p,
    current_i, last_it, periodic_x, periodic_y,
    ns_advection, free_surface)
    @unpack Re, τ = num
    @unpack nx, ny = grid

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b0_p = ones(grid)
    _b1_p = zeros(grid)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p, periodic_x, periodic_y)
    b0_p = Diagonal(vec(_b0_p))

    # Adjoint projection step
    @views adj_ph.uD[current_i,:] .= J_u(ph.u, current_i, grid_u)
    @views adj_ph.vD[current_i,:] .= J_v(ph.v, current_i, grid_v)
    if !last_it
        adj_ph.uD[current_i,:] .+= transpose(Bu) * adj_ph.ucorrD[current_i+1,:]
        adj_ph.vD[current_i,:] .+= transpose(Bv) * adj_ph.vcorrD[current_i+1,:]
    end

    # Adjoint pressure update step
    if !free_surface && !last_it
        GxT = τ .* transpose([opC_u.AxT * opC_u.Rx opC_u.Gx])
        GyT = τ .* transpose([opC_v.AyT * opC_u.Ry opC_v.Gy])

        @views adj_ph.pD[current_i,:] .+= adj_ph.pD[current_i+1,:]
        adj_ph.pD[current_i,:] .-= GxT * veci(adj_ph.ucorrD[current_i+1,:],grid_u,1)
        adj_ph.pD[current_i,:] .-= GyT * veci(adj_ph.vcorrD[current_i+1,:],grid_v,1)
    end

    # Adjoint Poisson equation
    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))
    GxT = τ .* transpose([iMu * opC_u.AxT * opC_u.Rx iMu * opC_u.Gx])
    GyT = τ .* transpose([iMv * opC_v.AyT * opC_u.Ry iMv * opC_v.Gy])

    rhs_ϕ = f2zeros(grid)
    @views rhs_ϕ .+= adj_ph.pD[current_i,:]
    @views rhs_ϕ .-= GxT * veci(adj_ph.uD[current_i,:],grid_u,1) .+ GyT * veci(adj_ph.vD[current_i,:],grid_v,1)

    blocks = DDM.decompose(transpose(Aϕ), grid.domdec, grid.domdec)
    # @views @mytime _, ch = bicgstabl!(adj_ph.ϕD[current_i,:], transpose(Aϕ), rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    adj_ph.ϕD[current_i,:] .= transpose(Aϕ) \ rhs_ϕ
    # println(ch)

    # Adjoint prediction step
    Du = iτ .* transpose([opC_p.AxT opC_p.Gx])
    Dv = iτ .* transpose([opC_p.AyT opC_p.Gy])

    # Du and Dv have to be modified updated into an equivalent
    # matrix for the adjoint system when periodic BCs are imposed.
    if periodic_x
        for i = 1:grid_u.ny
            II = CartesianIndex(i,1)
            JJ = CartesianIndex(i,grid_u.nx)
            pII_1 = lexicographic(II, grid.ny)
            pII_2 = lexicographic(II, grid_u.ny)
            pJJ = lexicographic(JJ, grid_u.ny)
            Du[pJJ,pII_1] = Du[pII_2,pII_1]
            Du[pII_2,pII_1] = 0.
        end
    end
    if periodic_y
        @inbounds @threads for i = 1:grid_v.nx
            II = CartesianIndex(1,i)
            JJ = CartesianIndex(grid_v.ny,i)
            pII_1 = lexicographic(II, grid.ny)
            pII_2 = lexicographic(II, grid_v.ny)
            pJJ = lexicographic(JJ, grid_v.ny)
            @inbounds Dv[pJJ,pII_1] = Dv[pII_2,pII_1]
            @inbounds Dv[pII_2,pII_1] = 0.0
        end
    end

    rhs_u = f2zeros(grid_u)
    rhs_u .+= adj_ph.uD[current_i,:]
    rhs_u .+= Du * veci(adj_ph.pD[current_i,:],grid,1)

    rhs_v = f2zeros(grid_v)
    rhs_v .+= adj_ph.vD[current_i,:]
    rhs_v .+= Dv * veci(adj_ph.pD[current_i,:],grid,1)

    if free_surface
        Smat = strain_rate(opC_u, opC_v)
        Su = transpose(iRe .* b0_p * [Smat[1,1] Smat[1,2]])
        Sv = transpose(iRe .* b0_p * [Smat[2,1] Smat[2,2]])

        rhs_u .+= Su * veci(adj_ph.pD[current_i,:],grid,2)
        if !last_it
            @views rhs_u .-= transpose(RlsFS_ucorr) * adj.u[current_i+1,:]
        end
        
        rhs_v .+= Sv * veci(adj_ph.pD[current_i,:],grid,2)
        if !last_it
            @views rhs_v .-= transpose(RlsFS_vcorr) * adj.u[current_i+1,:]
        end
    end

    blocks = DDM.decompose(transpose(Au), grid_u.domdec, grid_u.domdec)
    # @views @mytime _, ch = bicgstabl!(adj_ph.ucorrD[current_i,:], transpose(Au), rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    adj_ph.ucorrD[current_i,:] .= transpose(Au) \ rhs_u
    # println(ch)

    blocks = DDM.decompose(transpose(Av), grid_v.domdec, grid_v.domdec)
    # @views @mytime _, ch = bicgstabl!(adj_ph.vcorrD[current_i,:], transpose(Av), rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    adj_ph.vcorrD[current_i,:] .= transpose(Av) \ rhs_v
    # println(ch)

    return nothing
end