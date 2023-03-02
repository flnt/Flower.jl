function run_backward_discrete(num, grid, grid_u, grid_v,
    fwd, fwdS, fwdL, adj, adj_der, phS, phL,
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
    Vmean = false,
    verbose = false,
    show_every = 100,
    Ra = 0,
    λ = 1,
    ϵ_adj = 1e-8
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack MIXED, SOLID, LIQUID = ind
    @unpack RheatS_u, RheatL_u, RlsS_u, RlsS_TS, RlsS_TL = adj_der

    local LSA; local LSB;
    local LSAm1; local LSBm1;
    local AS; local AL;
    local BS; local BL;
    local ASm1; local ALm1;
    local BSm1; local BLm1;

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
    u .= fwd.u[end-1,:,:]

    @views phS.T .= fwdS.T[end,:,:]
    @views phL.T .= fwdL.T[end,:,:]
    @views phS.TD .= fwdS.TD[end,:]
    @views phL.TD .= fwdL.TD[end,:]

    if levelset
        update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)
    end

    if stefan
        update_stefan_velocity(num, grid, u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
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

    @views Rheat_T(num, grid, grid_u, grid_v, adj_der,
        phS.TD, phL.TD,
        fwd.u[current_i-1,:,:], fwd.u[current_i,:,:], LSA, LSB,
        CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean)
    
    if levelset
        rhs = J_u(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
        @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
        println(ch)
    end

    if heat
        if heat_solid_phase
            rhs = J_TS(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
            rhs .-= transpose(RlsS_TS) * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AS), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.phS.TD[current_i,:], transpose(AS), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
        if heat_liquid_phase
            rhs = J_TL(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
            rhs .-= transpose(RlsS_TL) * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AL), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.phL.TD[current_i,:], transpose(AL), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
    end

    current_i -= 1

    while current_i > 1
        @views u .= fwd.u[current_i-1,:,:]
        @views phS.T .= fwdS.T[current_i,:,:]
        @views phL.T .= fwdL.T[current_i,:,:]
        @views phS.TD .= fwdS.TD[current_i,:]
        @views phL.TD .= fwdL.TD[current_i,:]
    
        if levelset
            update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)
        end

        if stefan   
            update_stefan_velocity(num, grid, u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
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

        @views Rheat_u(num, grid, grid_u, grid_v, adj_der, fwd.u[current_i,:,:],
            fwdS.TD[current_i,:], fwdS.TD[current_i+1,:], ASm1, BSm1, opC_TS, BC_TS,
            fwdL.TD[current_i,:], fwdL.TD[current_i+1,:], ALm1, BLm1, opC_TL, BC_TL,
            fwd.u[current_i,:,:], fwd.u[current_i+1,:,:], LSAm1, LSBm1,
            CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean)

        @views Rheat_T(num, grid, grid_u, grid_v, adj_der,
            phS.TD, phL.TD,
            fwd.u[current_i-1,:,:], fwd.u[current_i,:,:], LSA, LSB,
            CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean)

        
        if levelset
            rhs = J_u(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
            rhs .-= transpose(RheatS_u) * adj.phS.TD[current_i+1,:]
            rhs .-= transpose(RheatL_u) * adj.phL.TD[current_i+1,:]
            rhs .-= transpose(RlsS_u) * adj.u[current_i+1,:]
            @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
            println(ch)
            LSAm1 .= LSA
            LSBm1 .= LSB
        end
    
        if heat
            if heat_solid_phase
                rhs = J_TS(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
                rhs .-= transpose(RlsS_TS) * adj.u[current_i,:]
                rhs .+= transpose(BSm1) * adj.phS.TD[current_i+1,:]
                kill_dead_cells!(veci(rhs, grid, 1), grid, geoS)
                kill_dead_cells!(veci(rhs, grid, 2), grid, geoS)
                ASm1 .= AS
                BSm1 .= BS
                @mytime blocks = DDM.decompose(transpose(AS), grid.domdec, grid.domdec)
    
                @views @mytime (_, ch) = bicgstabl!(adj.phS.TD[current_i,:], transpose(AS), rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
            end
            if heat_liquid_phase
                rhs = J_TL(num, grid, fwdS.TD, fwdL.TD, u, current_i-1)
                rhs .-= transpose(RlsS_TL) * adj.u[current_i,:]
                rhs .+= transpose(BLm1) * adj.phL.TD[current_i+1,:]
                kill_dead_cells!(veci(rhs, grid, 1), grid, geoL)
                kill_dead_cells!(veci(rhs, grid, 2), grid, geoL)
                ALm1 .= AL
                BLm1 .= BL
                @mytime blocks = DDM.decompose(transpose(AL), grid.domdec, grid.domdec)
    
                @views @mytime (_, ch) = bicgstabl!(adj.phL.TD[current_i,:], transpose(AL), rhs, Pl=ras(blocks,grid.pou), log=true)
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
                        normTS = norm(adj.phS.TD[current_i,:])
                        normTL = norm(adj.phL.TD[current_i,:])
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
        return MIXED, SOLID, LIQUID
    else
        return MIXED
    end
end