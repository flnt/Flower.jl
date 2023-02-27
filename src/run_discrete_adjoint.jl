function run_backward_discrete(num, grid, grid_u, grid_v,
    fwd, adj, phS, phL,
    opC_TS, opC_TL,
    J_TS, J_TL, J_u;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_TL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    stefan = false,
    heat = false,
    levelset = true,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    verbose = false,
    show_every = 100,
    Ra = 0,
    λ = 1,
    ϵ_adj = 1e-8
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack usave, uusave, uvsave, TSsave, TLsave, TDSsave, TDLsave = fwd

    local LSA; local LSB;
    local LSAm1; local LSBm1;
    local AS; local AL;
    local BS; local BL;
    local ASm1; local ALm1;
    local BSm1; local BLm1;

    local MIXED; local SOLID; local LIQUID;
    local MIXED_vel_ext; local SOLID_vel_ext; local LIQUID_vel_ext;
    local MIXED_u; local SOLID_u; local LIQUID_u;
    local MIXED_v; local SOLID_v; local LIQUID_v;

    if periodic_x
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
    end

    if periodic_y
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
    end

    CFL_sc = τ / Δ^2

    current_i = max_iterations + 1
    u .= fwd.usave[end-1,:,:]

    @views phS.T .= fwd.TSsave[end,:,:]
    @views phL.T .= fwd.TLsave[end,:,:]
    @views phS.TD .= fwd.TDSsave[end,:]
    @views phL.TD .= fwd.TDLsave[end,:]

    if levelset
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

        get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
        postprocess_grids!(grid, grid_u, grid_v, periodic_x, periodic_y, ϵ)
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        field_extension!(grid, u, κ, indices_vel_ext, NB, periodic_x, periodic_y)
    end

    if stefan
        Stefan_velocity!(num, grid, V, phS.T, phL.T, MIXED, periodic_x, periodic_y)
        V[MIXED] .*= 1. ./ λ
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        velocity_extension!(grid, u, V, indices_vel_ext, NB, periodic_x, periodic_y)
    end

    if levelset
        IIOE(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
        LSAm1 = copy(LSA)
        LSBm1 = copy(LSB)
    end

    if heat
        if heat_solid_phase
            AS, BS, _ = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, MIXED, geoS.projection,
                                    periodic_x, periodic_y)
            ASm1 = copy(AS)
            BSm1 = copy(BS)
        end
        if heat_liquid_phase
            AL, BL, _ = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, MIXED, geoL.projection,
                                    periodic_x, periodic_y)
            ALm1 = copy(AL)
            BLm1 = copy(BL)
        end
    end

    @views R3_TS_T, R3_TL_T = R_qi1(num, grid, grid_u, grid_v, 
        phS.TD, phL.TD,
        fwd.usave[current_i-1,:,:], fwd.usave[current_i,:,:], LSA, LSB,
        CFL_sc, periodic_x, periodic_y, ϵ_adj, λ)
    
    if levelset
        rhs = J_u(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
        @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
        println(ch)
    end

    if heat
        if heat_solid_phase
            rhs = J_TS(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
            rhs .-= R3_TS_T * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AS), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.TDS[current_i,:], transpose(AS), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
        if heat_liquid_phase
            rhs = J_TL(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
            rhs .-= R3_TL_T * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AL), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.TDL[current_i,:], transpose(AL), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
    end

    current_i -= 1

    while current_i > 1
        @views u .= fwd.usave[current_i-1,:,:]
        @views phS.T .= fwd.TSsave[current_i,:,:]
        @views phL.T .= fwd.TLsave[current_i,:,:]
        @views phS.TD .= fwd.TDSsave[current_i,:]
        @views phL.TD .= fwd.TDLsave[current_i,:]
    
        if levelset
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
            field_extension!(grid, u, κ, indices_vel_ext, NB, periodic_x, periodic_y)
        end

        if stefan   
            Stefan_velocity!(num, grid, V, phS.T, phL.T, MIXED, periodic_x, periodic_y)
            V[MIXED] .*= 1. ./ λ
            _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
            _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
            _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
            indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
            velocity_extension!(grid, u, V, indices_vel_ext, NB, periodic_x, periodic_y)
        end

        if levelset
            IIOE(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
        end

        if heat
            if heat_solid_phase
                AS, BS, _ = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, MIXED, geoS.projection,
                                    periodic_x, periodic_y)
            end
            if heat_liquid_phase
                AL, BL, _ = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, MIXED, geoL.projection,
                                    periodic_x, periodic_y)
            end
        end

        @views R1_u_T, R2_u_T, R3_u_T = R_qi(num, grid, grid_u, grid_v, fwd.usave[current_i,:,:],
            fwd.TDSsave[current_i,:], fwd.TDSsave[current_i+1,:], ASm1, BSm1, opC_TS, BC_TS,
            fwd.TDLsave[current_i,:], fwd.TDLsave[current_i+1,:], ALm1, BLm1, opC_TL, BC_TL,
            fwd.usave[current_i,:,:], fwd.usave[current_i+1,:,:], LSAm1, LSBm1,
            CFL_sc, periodic_x, periodic_y, ϵ_adj, λ)

        @views R3_TS_T, R3_TL_T = R_qi1(num, grid, grid_u, grid_v, 
            phS.TD, phL.TD,
            fwd.usave[current_i-1,:,:], fwd.usave[current_i,:,:], LSA, LSB,
            CFL_sc, periodic_x, periodic_y, ϵ_adj, λ)

        
        if levelset
            rhs = J_u(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
            rhs .-= R1_u_T * adj.TDS[current_i+1,:]
            rhs .-= R2_u_T * adj.TDL[current_i+1,:]
            rhs .-= R3_u_T * adj.u[current_i+1,:]
            @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
            println(ch)
            LSAm1 .= LSA
            LSBm1 .= LSB
        end
    
        if heat
            if heat_solid_phase
                rhs = J_TS(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
                rhs .-= R3_TS_T * adj.u[current_i,:]
                rhs .+= transpose(BSm1) * adj.TDS[current_i+1,:]
                kill_dead_cells!(veci(rhs, grid, 1), grid, geoS)
                kill_dead_cells!(veci(rhs, grid, 2), grid, geoS)
                ASm1 .= AS
                BSm1 .= BS
                @mytime blocks = DDM.decompose(transpose(AS), grid.domdec, grid.domdec)
    
                @views @mytime (_, ch) = bicgstabl!(adj.TDS[current_i,:], transpose(AS), rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
            end
            if heat_liquid_phase
                rhs = J_TL(num, grid, fwd.TDSsave, fwd.TDLsave, u, current_i-1)
                rhs .-= R3_TL_T * adj.u[current_i,:]
                rhs .+= transpose(BLm1) * adj.TDL[current_i+1,:]
                kill_dead_cells!(veci(rhs, grid, 1), grid, geoL)
                kill_dead_cells!(veci(rhs, grid, 2), grid, geoL)
                ALm1 .= AL
                BLm1 .= BL
                @mytime blocks = DDM.decompose(transpose(AL), grid.domdec, grid.domdec)
    
                @views @mytime (_, ch) = bicgstabl!(adj.TDL[current_i,:], transpose(AL), rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    if (heat && length(MIXED) != 0)
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                        normTL = norm(adj.TDL[current_i,:])
                        normTS = norm(adj.TDS[current_i,:])
                        normu = norm(adj.u[current_i,:])
                        print("$(@sprintf("norm(ψ) %.6e", normu))\t$(@sprintf("norm(ΘS) %.6e", normTS))\t$(@sprintf("norm(ΘL) %.6e", normTL))\n")
                    end
                catch
                    @show (MIXED)
                end
            end
        end

        current_i -= 1
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if heat
                print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
                print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            end
            print("\n\n")
        catch
            @show (length(MIXED))
        end
    end
    
    if levelset
        return MIXED, MIXED_u, MIXED_v, SOLID, LIQUID
    else
        return MIXED
    end
end