function R_qi(num, grid, grid_u, grid_v, um1,
    TD0_S, TD1_S, A_S, B_S, opC_TS, BC_TS,
    TD0_L, TD1_L, A_L, B_L, opC_TL, BC_TL,
    u0, u1, LSA, LSB,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, V, iso, faces, geoS, geoL = grid
    uj = copy(um1)

    TS = reshape(veci(TD1_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD1_L, grid, 1), (ny, nx))

    # res_TS = spdiagm(2ny*nx, ny*nx, 0 => zeros(ny*nx))
    # res_TS.nzval .= 0.
    res_TS = zeros(2ny*nx, ny*nx)

    derA_S = copy(A_S)
    derA_S.nzval .= 0.
    derB_S = copy(B_S)
    derB_S.nzval .= 0.

    # res_TL = spdiagm(2ny*nx, ny*nx, 0 => zeros(ny*nx))
    # res_TL.nzval .= 0.
    res_TL = zeros(2ny*nx, ny*nx)

    derA_L = copy(A_L)
    derA_L.nzval .= 0.
    derB_L = copy(B_L)
    derB_L.nzval .= 0.

    # res_u = copy(LSA)
    # res_u.nzval .= 0.
    res_u = zeros(ny*nx, ny*nx)

    derLSA = copy(LSA)
    derLSA.nzval .= 0.
    derLSB = copy(LSB)
    derLSB.nzval .= 0.
    LSAj = copy(LSA)
    LSBj = copy(LSB)

    rows_T = rowvals(A_S)
    rows_u = rowvals(LSA)
    @inbounds for JJ in ind.all_indices
        j = lexicographic(JJ, ny)
        uj[JJ] += ϵ_adj

        # Compute capacities
        grid.α .= NaN
        grid_u.α .= NaN
        grid_v.α .= NaN
        faces .= 0.
        grid_u.faces .= 0.
        grid_v.faces .= 0.
        grid.mid_point .= [Point(0.0, 0.0)]
        grid_u.mid_point .= [Point(0.0, 0.0)]
        grid_v.mid_point .= [Point(0.0, 0.0)]
        
        marching_squares!(num, grid, uj, periodic_x, periodic_y)
        interpolate_scalar!(grid, grid_u, grid_v, uj, grid_u.u, grid_v.u)
        marching_squares!(num, grid_u, grid_u.u, periodic_x, periodic_y)
        marching_squares!(num, grid_v, grid_v.u, periodic_x, periodic_y)

        MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.all_indices, nx, ny, periodic_x, periodic_y)
        MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

        get_interface_location!(grid, MIXED, periodic_x, periodic_y)
        get_interface_location!(grid_u, MIXED_u, periodic_x, periodic_y)
        get_interface_location!(grid_v, MIXED_v, periodic_x, periodic_y)
        get_interface_location_borders!(grid_u, grid_u.u, periodic_x, periodic_y)
        get_interface_location_borders!(grid_v, grid_v.u, periodic_x, periodic_y)

        geoL.emptied .= false
        geoS.emptied .= false
        grid_u.geoL.emptied .= false
        grid_u.geoS.emptied .= false
        grid_v.geoL.emptied .= false
        grid_v.geoS.emptied .= false

        get_curvature(num, grid, uj, MIXED, periodic_x, periodic_y)
        postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, ϵ)
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        field_extension!(grid, uj, grid.κ, indices_vel_ext, NB, periodic_x, periodic_y)

        Stefan_velocity!(num, grid, V, TS, TL, MIXED, periodic_x, periodic_y)
        V[MIXED] .*= 1. ./ λ
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        velocity_extension!(grid, uj, V, indices_vel_ext, NB, periodic_x, periodic_y)

        # get perturbed matrices
        Aj_S, Bj_S, _ = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, MIXED, geoS.projection,
                                periodic_x, periodic_y)
        derA_S .= (Aj_S .- A_S) ./ ϵ_adj
        derB_S .= (Bj_S .- B_S) ./ ϵ_adj
        Rj_S = derA_S * TD1_S .- derB_S * TD0_S

        # if JJ[2] == 4
        #     @show (JJ, j)
        #     println(Rj_S[j-3:j+7])
        #     # println(A_S[j,:])
        #     # println()
        #     # println(Aj_S[j,:])
        #     # println()
        #     # println(derA_S[j,:])
        #     println()
        # end

        Aj_L, Bj_L, _ = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, MIXED, geoL.projection,
                                periodic_x, periodic_y)
        derA_L .= (Aj_L .- A_L) ./ ϵ_adj
        derB_L .= (Bj_L .- B_L) ./ ϵ_adj
        Rj_L = derA_L * TD1_L .- derB_L * TD0_L

        IIOE(grid, LSAj, LSBj, uj, V, CFL_sc, periodic_x, periodic_y)
        derLSA .= (LSAj .- LSA) ./ ϵ_adj
        derLSB .= (LSBj .- LSB) ./ ϵ_adj
        utmp = zeros(ny, nx)
        utmp[JJ] = 1.0
        Rj_u = derLSA * vec(u1) .- derLSB * vec(u0) .- LSB * vec(utmp)

        # for i in nzrange(A_S, j)
        #     @inbounds row = rows_T[i]
        #     @inbounds res_TS[row,j] = Rj_S[row]
        #     @inbounds res_TL[row,j] = Rj_L[row]
        # end
        # for i in nzrange(LSA, j)
        #     @inbounds row = rows_u[i]
        #     @inbounds res_u[row,j] = Rj_u[row]
        # end
        @inbounds res_TS[:,j] .= Rj_S
        @inbounds res_TL[:,j] .= Rj_L
        @inbounds res_u[:,j] .= Rj_u

        uj .= um1
    end

    return transpose(res_TS), transpose(res_TL), transpose(res_u)
end

function R_qi1(num, grid, grid_u, grid_v, 
    TD_S, TD_L,
    u0, u1, LSA, LSB,
    CFL_sc, periodic_x, periodic_y, ϵ_adj, λ)

    @unpack ϵ, NB = num
    @unpack nx, ny, ind, u, V, faces, iso, geoS, geoL = grid

    TS = reshape(veci(TD_S, grid, 1), (ny, nx))
    TL = reshape(veci(TD_L, grid, 1), (ny, nx))
    TSj = copy(TS)
    TLj = copy(TL)

    # res_TS = spdiagm(ny*nx, 2ny*nx, 0 => zeros(ny*nx))
    # res_TS.nzval .= 0.
    # res_TL = spdiagm(ny*nx, 2ny*nx, 0 => zeros(ny*nx))
    # res_TL.nzval .= 0.
    res_TS = zeros(ny*nx, 2ny*nx)
    res_TL = zeros(ny*nx, 2ny*nx)
    
    derLSA_S = copy(LSA)
    derLSA_S.nzval .= 0.
    derLSB_S = copy(LSB)
    derLSB_S.nzval .= 0.

    derLSA_L = copy(LSA)
    derLSA_L.nzval .= 0.
    derLSB_L = copy(LSB)
    derLSB_L.nzval .= 0.

    LSAj = copy(LSA)
    LSBj = copy(LSB)

    rows = rowvals(LSA)
    @inbounds for j = 1:ny*nx
        TSj[j] += ϵ_adj
        TLj[j] += ϵ_adj

        # Compute capacities
        grid.α .= NaN
        grid_u.α .= NaN
        grid_v.α .= NaN
        faces .= 0.
        grid_u.faces .= 0.
        grid_v.faces .= 0.
        grid.mid_point .= [Point(0.0, 0.0)]
        grid_u.mid_point .= [Point(0.0, 0.0)]
        grid_v.mid_point .= [Point(0.0, 0.0)]
        
        marching_squares!(num, grid, u, periodic_x, periodic_y)
        interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
        marching_squares!(num, grid_u, grid_u.u, periodic_x, periodic_y)
        marching_squares!(num, grid_v, grid_v.u, periodic_x, periodic_y)

        MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.all_indices, nx, ny, periodic_x, periodic_y)
        MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

        get_interface_location!(grid, MIXED, periodic_x, periodic_y)
        get_interface_location!(grid_u, MIXED_u, periodic_x, periodic_y)
        get_interface_location!(grid_v, MIXED_v, periodic_x, periodic_y)
        get_interface_location_borders!(grid_u, grid_u.u, periodic_x, periodic_y)
        get_interface_location_borders!(grid_v, grid_v.u, periodic_x, periodic_y)

        geoL.emptied .= false
        geoS.emptied .= false
        grid_u.geoL.emptied .= false
        grid_u.geoS.emptied .= false
        grid_v.geoL.emptied .= false
        grid_v.geoS.emptied .= false

        get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
        postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, ϵ)
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        field_extension!(grid, u, grid.κ, indices_vel_ext, NB, periodic_x, periodic_y)

        # get perturbed matrices
        Stefan_velocity!(num, grid, V, TSj, TL, MIXED, periodic_x, periodic_y)
        V[MIXED] .*= 1. ./ λ
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        velocity_extension!(grid, u, V, indices_vel_ext, NB, periodic_x, periodic_y)

        IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
        derLSA_S .= (LSAj .- LSA) ./ ϵ_adj
        derLSB_S .= (LSBj .- LSB) ./ ϵ_adj

        Rj_S = derLSA_S * vec(u1) .- derLSB_S * vec(u0)

        Stefan_velocity!(num, grid, V, TS, TLj, MIXED, periodic_x, periodic_y)
        V[MIXED] .*= 1. ./ λ
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        velocity_extension!(grid, u, V, indices_vel_ext, NB, periodic_x, periodic_y)

        IIOE(grid, LSAj, LSBj, u, V, CFL_sc, periodic_x, periodic_y)
        derLSA_L .= (LSAj .- LSA) ./ ϵ_adj
        derLSB_L .= (LSBj .- LSB) ./ ϵ_adj

        Rj_L = derLSA_L * vec(u1) .- derLSB_L * vec(u0)

        # for i in nzrange(LSA, j)
        #     @inbounds row = rows[i]
        #     @inbounds res_TS[row,j] = Rj_S[row]
        #     @inbounds res_TL[row,j] = Rj_L[row]
        # end
        @inbounds res_TS[:,j] .= Rj_S
        @inbounds res_TL[:,j] .= Rj_L

        TSj .= TS
        TLj .= TL
    end

    return transpose(res_TS), transpose(res_TL)
end