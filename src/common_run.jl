function update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)
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

    grid.ind.MIXED_ext, grid.ind.SOLID_ext, grid.ind.LIQUID_ext = get_cells_indices(grid.iso, grid.ind.all_indices, grid.nx, grid.ny, periodic_x, periodic_y)
    grid_u.ind.MIXED_ext, grid_u.ind.SOLID_ext, grid_u.ind.LIQUID_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
    grid_v.ind.MIXED_ext, grid_v.ind.SOLID_ext, grid_v.ind.LIQUID_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
    grid.ind.MIXED, grid.ind.SOLID, grid.ind.LIQUID = get_cells_indices(grid.iso, grid.ind.all_indices)
    grid_u.ind.MIXED, grid_u.ind.SOLID, grid_u.ind.LIQUID = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
    grid_v.ind.MIXED, grid_v.ind.SOLID, grid_v.ind.LIQUID = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

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

    get_curvature(num, grid, u, grid.ind.MIXED, periodic_x, periodic_y)
    postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, num.ϵ)

    _MIXED_L_ext = intersect(findall(grid.geoL.emptied), grid.ind.MIXED_ext)
    _MIXED_S_ext = intersect(findall(grid.geoS.emptied), grid.ind.MIXED_ext)
    _MIXED_ext = vcat(_MIXED_L_ext, _MIXED_S_ext)
    indices_ext = vcat(grid.ind.SOLID_ext, _MIXED_ext, grid.ind.LIQUID_ext)
    field_extension!(grid, u, grid.κ, indices_ext, num.NB, periodic_x, periodic_y)
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

