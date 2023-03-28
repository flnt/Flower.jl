function update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)
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
    get_interface_location_borders!(grid_u, grid_u.u, periodic_x, periodic_y)
    get_interface_location_borders!(grid_v, grid_v.u, periodic_x, periodic_y)

    grid.geoL.emptied .= false
    grid.geoS.emptied .= false
    grid_u.geoL.emptied .= false
    grid_u.geoS.emptied .= false
    grid_v.geoL.emptied .= false
    grid_v.geoS.emptied .= false

    get_curvature(num, grid, u, κ, grid.ind.MIXED, periodic_x, periodic_y)
    postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, num.ϵ)

    _MIXED_L_ext = intersect(findall(grid.geoL.emptied), grid.ind.MIXED_ext)
    _MIXED_S_ext = intersect(findall(grid.geoS.emptied), grid.ind.MIXED_ext)
    _MIXED_ext = vcat(_MIXED_L_ext, _MIXED_S_ext)
    indices_ext = vcat(grid.ind.SOLID_ext, _MIXED_ext, grid.ind.LIQUID_ext)
    field_extension!(grid, u, κ, indices_ext, num.NB, periodic_x, periodic_y)

    return NB_indices
end

function update_stefan_velocity(num, grid, u, TS, TL, periodic_x, periodic_y, λ, Vmean)
    Stefan_velocity!(num, grid, grid.V, TS, TL, grid.ind.MIXED, periodic_x, periodic_y)
    grid.V[grid.ind.MIXED] .*= 1. ./ λ
    if Vmean
        a = mean(grid.V[grid.ind.MIXED])
        grid.V[grid.ind.MIXED] .= a
    end
    _MIXED_L_ext = intersect(findall(grid.geoL.emptied), grid.ind.MIXED_ext)
    _MIXED_S_ext = intersect(findall(grid.geoS.emptied), grid.ind.MIXED_ext)
    _MIXED_ext = vcat(_MIXED_L_ext, _MIXED_S_ext)
    indices_ext = vcat(grid.ind.SOLID_ext, _MIXED_ext, grid.ind.LIQUID_ext)
    velocity_extension!(grid, u, grid.V, indices_ext, num.NB, periodic_x, periodic_y)
end

function update_free_surface_velocity(num, grid_u, grid_v, uD, vD, periodic_x, periodic_y)
    grid_u.V .= reshape(veci(uD,grid_u,2), (grid_u.ny, grid_u.nx))
    grid_v.V .= reshape(veci(vD,grid_v,2), (grid_v.ny, grid_v.nx))

    _MIXED_L_u_ext = intersect(findall(grid_u.geoL.emptied),
                                    grid_u.ind.MIXED_ext)
    _MIXED_S_u_ext = intersect(findall(grid_u.geoS.emptied),
                                    grid_u.ind.MIXED_ext)
    _MIXED_u_ext = vcat(_MIXED_L_u_ext, _MIXED_S_u_ext)
    indices_u_ext = vcat(grid_u.ind.SOLID_ext, _MIXED_u_ext, grid_u.ind.LIQUID_ext)

    _MIXED_L_v_ext = intersect(findall(grid_v.geoL.emptied),
                                    grid_v.ind.MIXED_ext)
    _MIXED_S_v_ext = intersect(findall(grid_v.geoS.emptied),
                                    grid_v.ind.MIXED_ext)
    _MIXED_v_ext = vcat(_MIXED_L_v_ext, _MIXED_S_v_ext)
    indices_v_ext = vcat(grid_v.ind.SOLID_ext, _MIXED_v_ext, grid_v.ind.LIQUID_ext)

    velocity_extension!(grid_u, grid_u.u, grid_u.V, indices_u_ext, num.NB, periodic_x, periodic_y)
    velocity_extension!(grid_v, grid_v.u, grid_v.V, indices_v_ext, num.NB, periodic_x, periodic_y)
end

function adjoint_projection_fs(num, grid, grid_u, grid_v,
    adj, adj_ph, RlsFS_ucorr, RlsFS_vcorr,
    Au, Bu, Av, Bv, Aϕ,
    opC_p, opC_u, opC_v, BC_p, current_i, last_it)
    @unpack Re, τ = num
    @unpack nx, ny = grid

    iRe = 1.0 / Re
    iτ = 1.0 / τ

    a0_p = zeros(grid)
    _a1_p = zeros(grid)
    _b0_p = ones(grid)
    _b1_p = zeros(grid)
    set_borders!(grid, a0_p, _a1_p, _b0_p, _b1_p, BC_p)
    b0_p = Diagonal(vec(_b0_p))

    # Adjoint projection step
    @views adj_ph.u[current_i,:] .= 0
    @views adj_ph.v[current_i,:] .= 0
    if !last_it
        adj_ph.u[current_i,:] .+= veci(transpose(Bu) * adj_ph.ucorrD[current_i+1,:], grid_u, 1)
        adj_ph.v[current_i,:] .+= veci(transpose(Bv) * adj_ph.vcorrD[current_i+1,:], grid_v, 1)
    end

    # Adjoint Poisson equation
    iMu = Diagonal(1 ./ (opC_u.M.diag .+ eps(0.01)))
    iMv = Diagonal(1 ./ (opC_v.M.diag .+ eps(0.01)))
    GxT = τ .* transpose([iMu * opC_u.AxT * opC_u.Rx iMu * opC_u.Gx])
    GyT = τ .* transpose([iMv * opC_v.AyT * opC_u.Ry iMv * opC_v.Gy])

    rhs_ϕ = zeros(2*ny*nx)
    @views rhs_ϕ .-= GxT * adj_ph.u[current_i,:] .+ GyT * adj_ph.v[current_i,:]

    blocks = DDM.decompose(transpose(Aϕ), grid.domdec, grid.domdec)
    @views @mytime _, ch = bicgstabl!(adj_ph.pD[current_i,:], transpose(Aϕ), rhs_ϕ, Pl=ras(blocks,grid.pou), log=true)
    println(ch)

    # Adjoint prediction step
    Du = transpose(iτ .* [opC_p.AxT opC_p.Gx])
    Dv = transpose(iτ .* [opC_p.AyT opC_p.Gy])

    Smat = strain_rate(opC_u, opC_v)
    Su = transpose(iRe .* b0_p * [Smat[1,1] Smat[1,2]])
    Sv = transpose(iRe .* b0_p * [Smat[2,1] Smat[2,2]])

    rhs_u = zeros(2*grid_u.ny*grid_u.nx)
    veci(rhs_u,grid_u,1) .+= adj_ph.u[current_i,:]
    rhs_u .+= Du * veci(adj_ph.pD[current_i,:],grid,1)
    rhs_u .+= Su * veci(adj_ph.pD[current_i,:],grid,2)
    if !last_it
        @views rhs_u .-= transpose(RlsFS_ucorr) * adj.u[current_i+1,:]
    end
    
    rhs_v = zeros(2*grid_v.ny*grid_v.nx)
    veci(rhs_v,grid_v,1) .+= adj_ph.v[current_i,:]
    rhs_v .+= Dv * veci(adj_ph.pD[current_i,:],grid,1)
    rhs_v .+= Sv * veci(adj_ph.pD[current_i,:],grid,2)
    if !last_it
        @views rhs_v .-= transpose(RlsFS_vcorr) * adj.u[current_i+1,:]
    end

    blocks = DDM.decompose(transpose(Au), grid_u.domdec, grid_u.domdec)
    @views @mytime _, ch = bicgstabl!(adj_ph.ucorrD[current_i,:], transpose(Au), rhs_u, Pl=ras(blocks,grid_u.pou), log=true)
    println(ch)

    blocks = DDM.decompose(transpose(Av), grid_v.domdec, grid_v.domdec)
    @views @mytime _, ch = bicgstabl!(adj_ph.vcorrD[current_i,:], transpose(Av), rhs_v, Pl=ras(blocks,grid_v.pou), log=true)
    println(ch)

    return nothing
end