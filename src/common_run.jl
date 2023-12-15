function update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y, empty = true)
    grid.α .= NaN
    grid_u.α .= NaN
    grid_v.α .= NaN
    grid.faces .= 0.
    grid_u.faces .= 0.
    grid_v.faces .= 0.
    grid.mid_point .= [Point(0.0, 0.0)]
    grid_u.mid_point .= [Point(0.0, 0.0)]
    grid_v.mid_point .= [Point(0.0, 0.0)]

    marching_squares!(num, grid, u, periodic_x, periodic_y)
    interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
    marching_squares!(num, grid_u, grid_u.u, periodic_x, periodic_y)
    marching_squares!(num, grid_v, grid_v.u, periodic_x, periodic_y)

    NB_indices_base = get_NB_width_indices_base(num.NB)

    grid.ind.MIXED_ext, grid.ind.SOLID_ext, grid.ind.LIQUID_ext = get_cells_indices(grid.iso, grid.ind.all_indices, grid.nx, grid.ny, periodic_x, periodic_y)
    grid_u.ind.MIXED_ext, grid_u.ind.SOLID_ext, grid_u.ind.LIQUID_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
    grid_v.ind.MIXED_ext, grid_v.ind.SOLID_ext, grid_v.ind.LIQUID_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
    grid.ind.MIXED, grid.ind.SOLID, grid.ind.LIQUID = get_cells_indices(grid.iso, grid.ind.all_indices)
    grid_u.ind.MIXED, grid_u.ind.SOLID, grid_u.ind.LIQUID = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
    grid_v.ind.MIXED, grid_v.ind.SOLID, grid_v.ind.LIQUID = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

    NB_indices = get_NB_width(grid, grid.ind.MIXED, NB_indices_base)

    get_interface_location!(grid, grid.ind.MIXED, periodic_x, periodic_y)
    get_interface_location!(grid_u, grid_u.ind.MIXED, periodic_x, periodic_y)
    get_interface_location!(grid_v, grid_v.ind.MIXED, periodic_x, periodic_y)

    grid.geoL.emptied .= false
    grid.geoS.emptied .= false
    grid_u.geoL.emptied .= false
    grid_u.geoS.emptied .= false
    grid_v.geoL.emptied .= false
    grid_v.geoS.emptied .= false

    κ .= 0.0
    grid_u.κ .= 0.0
    grid_v.κ .= 0.0
    get_curvature(num, grid, u, κ, grid.ind.MIXED, periodic_x, periodic_y)
    get_curvature(num, grid_u, grid_u.u, grid_u.κ, grid_u.ind.MIXED, periodic_x, periodic_y)
    get_curvature(num, grid_v, grid_v.u, grid_v.κ, grid_v.ind.MIXED, periodic_x, periodic_y)
    postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, num.ϵ, empty)

    _MIXED_L = intersect(findall(grid.geoL.emptied), grid.ind.MIXED)
    _MIXED_S = intersect(findall(grid.geoS.emptied), grid.ind.MIXED)
    _MIXED = vcat(_MIXED_L, _MIXED_S)
    indices_ext1 = vcat(grid.ind.SOLID, _MIXED, grid.ind.LIQUID)

    if periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_left[1][2:end-1], grid.ind.b_bottom[1][2:end-1],
            grid.ind.b_right[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif !periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_bottom[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif periodic_x && !periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_left[1][2:end-1], grid.ind.b_right[1][2:end-1]
        ))
    else
        indices_ext = intersect(indices_ext1, vec(grid.ind.inside))
    end
    left_ext = intersect(indices_ext1, grid.ind.b_left[1][2:end-1])
    bottom_ext = intersect(indices_ext1, grid.ind.b_bottom[1][2:end-1])
    right_ext = intersect(indices_ext1, grid.ind.b_right[1][2:end-1])
    top_ext = intersect(indices_ext1, grid.ind.b_top[1][2:end-1])
    field_extension!(grid, u, κ, indices_ext, left_ext, bottom_ext, right_ext, top_ext, num.NB, periodic_x, periodic_y)

    locate_contact_line!(grid)
    locate_contact_line!(grid_u)
    locate_contact_line!(grid_v)

    # Apply inhomogeneous BC to more than one grid point by extending the contact line
    extend_contact_line!(grid_u, num.n_ext_cl)
    extend_contact_line!(grid_v, num.n_ext_cl)

    return NB_indices
end

function update_stefan_velocity(num, grid, u, TS, TL, periodic_x, periodic_y, λ, Vmean)
    Stefan_velocity!(num, grid, grid.V, TS, TL, grid.ind.MIXED, periodic_x, periodic_y)
    grid.V[grid.ind.MIXED] .*= 1. ./ λ
    if Vmean
        a = mean(grid.V[grid.ind.MIXED])
        grid.V[grid.ind.MIXED] .= a
    end
    
    _MIXED_L = intersect(findall(grid.geoL.emptied), grid.ind.MIXED)
    _MIXED_S = intersect(findall(grid.geoS.emptied), grid.ind.MIXED)
    _MIXED = vcat(_MIXED_L, _MIXED_S)
    indices_ext1 = vcat(grid.ind.SOLID, _MIXED, grid.ind.LIQUID)

    if periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_left[1][2:end-1], grid.ind.b_bottom[1][2:end-1],
            grid.ind.b_right[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif !periodic_x && periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_bottom[1][2:end-1], grid.ind.b_top[1][2:end-1]
        ))
    elseif periodic_x && !periodic_y
        indices_ext = intersect(indices_ext1, vcat(
            vec(grid.ind.inside), grid.ind.b_left[1][2:end-1], grid.ind.b_right[1][2:end-1]
        ))
    else
        indices_ext = intersect(indices_ext1, vec(grid.ind.inside))
    end
    left_ext = intersect(indices_ext1, grid.ind.b_left[1][2:end-1])
    bottom_ext = intersect(indices_ext1, grid.ind.b_bottom[1][2:end-1])
    right_ext = intersect(indices_ext1, grid.ind.b_right[1][2:end-1])
    top_ext = intersect(indices_ext1, grid.ind.b_top[1][2:end-1])
    field_extension!(grid, u, grid.V, indices_ext, left_ext, bottom_ext, right_ext, top_ext, num.NB, periodic_x, periodic_y)
end

function update_free_surface_velocity(num, grid_u, grid_v, uD, vD, periodic_x, periodic_y)
    grid_u.V .= reshape(veci(uD,grid_u,2), (grid_u.ny, grid_u.nx))
    grid_v.V .= reshape(veci(vD,grid_v,2), (grid_v.ny, grid_v.nx))

    _MIXED_u_L = intersect(findall(grid_u.geoL.emptied), grid_u.ind.MIXED)
    _MIXED_u_S = intersect(findall(grid_u.geoS.emptied), grid_u.ind.MIXED)
    _MIXED_u = vcat(_MIXED_u_L, _MIXED_u_S)
    indices_u_ext1 = vcat(grid_u.ind.SOLID, _MIXED_u, grid_u.ind.LIQUID)

    if periodic_x && periodic_y
        indices_u_ext = intersect(indices_u_ext1, vcat(
            vec(grid_u.ind.inside), grid_u.ind.b_left[1][2:end-1], grid_u.ind.b_bottom[1][2:end-1],
            grid_u.ind.b_right[1][2:end-1], grid_u.ind.b_top[1][2:end-1]
        ))
    elseif !periodic_x && periodic_y
        indices_u_ext = intersect(indices_u_ext1, vcat(
            vec(grid_u.ind.inside), grid_u.ind.b_bottom[1][2:end-1], grid_u.ind.b_top[1][2:end-1]
        ))
    elseif periodic_x && !periodic_y
        indices_u_ext = intersect(indices_u_ext1, vcat(
            vec(grid_u.ind.inside), grid_u.ind.b_left[1][2:end-1], grid_u.ind.b_right[1][2:end-1]
        ))
    else
        indices_u_ext = intersect(indices_u_ext1, vec(grid_u.ind.inside))
    end
    left_u_ext = intersect(indices_u_ext1, grid_u.ind.b_left[1][2:end-1])
    bottom_u_ext = intersect(indices_u_ext1, grid_u.ind.b_bottom[1][2:end-1])
    right_u_ext = intersect(indices_u_ext1, grid_u.ind.b_right[1][2:end-1])
    top_u_ext = intersect(indices_u_ext1, grid_u.ind.b_top[1][2:end-1])

    _MIXED_v_L = intersect(findall(grid_v.geoL.emptied), grid_v.ind.MIXED)
    _MIXED_v_S = intersect(findall(grid_v.geoS.emptied), grid_v.ind.MIXED)
    _MIXED_v = vcat(_MIXED_v_L, _MIXED_v_S)
    indices_v_ext1 = vcat(grid_v.ind.SOLID, _MIXED_v, grid_v.ind.LIQUID)

    if periodic_x && periodic_y
        indices_v_ext = intersect(indices_v_ext1, vcat(
            vec(grid_v.ind.inside), grid_v.ind.b_left[1][2:end-1], grid_v.ind.b_bottom[1][2:end-1],
            grid_v.ind.b_right[1][2:end-1], grid_v.ind.b_top[1][2:end-1]
        ))
    elseif !periodic_x && periodic_y
        indices_v_ext = intersect(indices_v_ext1, vcat(
            vec(grid_v.ind.inside), grid_v.ind.b_bottom[1][2:end-1], grid_v.ind.b_top[1][2:end-1]
        ))
    elseif periodic_x && !periodic_y
        indices_v_ext = intersect(indices_v_ext1, vcat(
            vec(grid_v.ind.inside), grid_v.ind.b_left[1][2:end-1], grid_v.ind.b_right[1][2:end-1]
        ))
    else
        indices_v_ext = intersect(indices_v_ext1, vec(grid_v.ind.inside))
    end
    left_v_ext = intersect(indices_v_ext1, grid_v.ind.b_left[1][2:end-1])
    bottom_v_ext = intersect(indices_v_ext1, grid_v.ind.b_bottom[1][2:end-1])
    right_v_ext = intersect(indices_v_ext1, grid_v.ind.b_right[1][2:end-1])
    top_v_ext = intersect(indices_v_ext1, grid_v.ind.b_top[1][2:end-1])

    field_extension!(grid_u, grid_u.u, grid_u.V, indices_u_ext, left_u_ext, bottom_u_ext, right_u_ext, top_u_ext, num.NB, periodic_x, periodic_y)
    field_extension!(grid_v, grid_v.u, grid_v.V, indices_v_ext, left_v_ext, bottom_v_ext, right_v_ext, top_v_ext, num.NB, periodic_x, periodic_y)
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