function run_backward_discrete(num, grid, grid_u, grid_v,
    fwd, fwdS, fwdL, adj, adj_der, phS, phL,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL;
    J_TS = x -> zero(phS.TD),
    J_TL = x -> zero(phL.TD),
    J_u = x -> fzeros(grid),
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
    BC_pS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_pL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_uS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_uL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_vS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()
    ),
    BC_vL = Boundaries(
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
    advection = false,
    heat = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    heat_convection = false,
    navier_stokes = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
    ns_advection = false,
    free_surface = false,
    Vmean = false,
    levelset = true,
    verbose = false,
    show_every = 100,
    Ra = 0,
    λ = 1,
    ϵ_adj = 1e-8
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re,
            max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack MIXED, SOLID, LIQUID = ind
    @unpack RheatS_ls, RheatL_ls, RlsS_ls, RlsS_TS, RlsS_TL,
            RuS_ls, RuL_ls, RvS_ls, RvL_ls, RpS_ls, RpL_ls,
            RucorrS_ls0, RucorrS_ls1, RvcorrS_ls0, RvcorrS_ls1,
            RucorrL_ls0, RucorrL_ls1, RvcorrL_ls0, RvcorrL_ls1,
            RlsFS_ls, RlsFS_ucorrS, RlsFS_vcorrS, RlsFS_ucorrL, RlsFS_vcorrL = adj_der

    local NB_indices;

    local LSA; local LSB;
    local LSAm1; local LSBm1;
    local AS; local AL;
    local BS; local BL;
    local ASm1; local ALm1;
    local BSm1; local BLm1;
    local AuS; local BuS;
    local AvS; local BvS;
    local AϕS;
    local AuL; local BuL;
    local AvL; local BvL;
    local AϕL;

    local BuSm1; local BvSm1;
    local BuLm1; local BvLm1;

    tmp_su = star2(grid_u)
    tmp_sv = star2(grid_v)
    BuS = [tmp_su tmp_su;
           tmp_su tmp_su]
    BvS = [tmp_sv tmp_sv;
           tmp_sv tmp_sv]
    BuL = [tmp_su tmp_su;
           tmp_su tmp_su]
    BvL = [tmp_sv tmp_sv;
           tmp_sv tmp_sv]

    BuSm1 = [tmp_su tmp_su;
             tmp_su tmp_su]
    BvSm1 = [tmp_sv tmp_sv;
             tmp_sv tmp_sv]
    BuLm1 = [tmp_su tmp_su;
             tmp_su tmp_su]
    BvLm1 = [tmp_sv tmp_sv;
             tmp_sv tmp_sv]

    local Mum1_S; local Mvm1_S
    local Mum1_L; local Mvm1_L

    θ_out = zeros(grid, 4)
    utmp = copy(u)

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

    iRe = 1.0 / Re
    CFL_sc = τ / Δ^2

    current_i = max_iterations + 1
    u0 = fwd.u[end-1,:,:]
    u1 = fwd.u[end,:,:]

    @views phS.T .= fwdS.T[end,:,:]
    @views phL.T .= fwdL.T[end,:,:]
    @views phS.TD .= fwdS.TD[end,:]
    @views phL.TD .= fwdL.TD[end,:]

    @views phS.uD .= fwdS.ucorrD[end-1,:]
    @views phL.uD .= fwdL.ucorrD[end-1,:]
    @views phS.vD .= fwdS.vcorrD[end-1,:]
    @views phL.vD .= fwdL.vcorrD[end-1,:]

    @views phS.ϕD .= fwdS.pD[end,:]
    @views phL.ϕD .= fwdL.pD[end,:]

    if levelset
        NB_indices = update_ls_data(num, grid, grid_u, grid_v, u0, κ, periodic_x, periodic_y)
        
        if advection && stefan
            update_stefan_velocity(num, grid, u0, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
            IIOE(grid, LSA, LSB, u0, V, CFL_sc, periodic_x, periodic_y)
        elseif advection && navier_stokes
            if free_surface && ns_liquid_phase
                update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)
            end
            if free_surface && ns_solid_phase
                update_free_surface_velocity(num, grid_u, grid_v, phS.uD, phS.vD, periodic_x, periodic_y)
            end
    
            level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            utmp .= reshape(gmres(LSA,(LSB*vec(u0))), (ny,nx))
            S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
        end

        LSAm1 = copy(LSA)
        LSBm1 = copy(LSB)
    end

    if heat
        if heat_solid_phase
            AS, BS, _ = set_heat!(dir, num, grid, opC_TS, geoS, phS, θd, BC_TS, MIXED, geoS.projection,
                                opS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                periodic_x, periodic_y, heat_convection)
            ASm1 = copy(AS)
            BSm1 = copy(BS)
        end
        if heat_liquid_phase
            AL, BL, _ = set_heat!(dir, num, grid, opC_TL, geoL, phL, θd, BC_TL, MIXED, geoL.projection,
                                opL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                periodic_x, periodic_y, heat_convection)
            ALm1 = copy(AL)
            BLm1 = copy(BL)
        end
    end

    if advection && stefan
        @views Rheat_q1(num, grid, grid_u, grid_v, adj_der,
            phS.TD, phL.TD,
            u0, u1, LSA, LSB,
            CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean,
            heat_solid_phase, heat_liquid_phase)
    end
    
    if levelset && advection && stefan
        rhs = J_u(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
        @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
        println(ch)
    end

    if heat
        if heat_solid_phase
            rhs = J_TS(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
            rhs .-= transpose(RlsS_TS) * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AS), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.phS.TD[current_i,:], transpose(AS), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
        if heat_liquid_phase
            rhs = J_TL(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
            rhs .-= transpose(RlsS_TL) * adj.u[current_i,:]
            @mytime blocks = DDM.decompose(transpose(AL), grid.domdec, grid.domdec)

            @views @mytime (_, ch) = bicgstabl!(adj.phL.TD[current_i,:], transpose(AL), rhs, Pl=ras(blocks,grid.pou), log=true)
            println(ch)
        end
    end

    if levelset && navier_stokes
        if ns_solid_phase
            _ = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                opC_pS, opC_uS, opC_vS,
                periodic_x, periodic_y, true)
        end
        if ns_liquid_phase
            _ = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                opC_pL, opC_uL, opC_vL,
                periodic_x, periodic_y, true)
        end

        Mum1_S = copy(opC_uS.M)
        Mvm1_S = copy(opC_vS.M)
        Mum1_L = copy(opC_uL.M)
        Mvm1_L = copy(opC_vL.M)

        NB_indices = update_ls_data(num, grid, grid_u, grid_v, u1, κ, periodic_x, periodic_y)
    end

    if navier_stokes && free_surface
        if ns_solid_phase
            AuS, BuS, _, AvS, BvS, _, AϕS, _ = set_navier_stokes(num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                                opC_pS, opC_uS, opC_vS, BC_pS, BC_uS, BC_vS,
                                                                Mum1_S, Mvm1_S, iRe,
                                                                opS, phS,
                                                                periodic_x, periodic_y, ns_advection)
            BuSm1 .= BuS
            BvSm1 .= BvS
        end
        if ns_liquid_phase
            AuL, BuL, _, AvL, BvL, _, AϕL, _ = set_navier_stokes(num, grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                                opC_pL, opC_uL, opC_vL, BC_pL, BC_uL, BC_vL,
                                                                Mum1_L, Mvm1_L, iRe,
                                                                opL, phL,
                                                                periodic_x, periodic_y, ns_advection)
            BuLm1 .= BuL
            BvLm1 .= BvL
        end
    end

    if levelset && navier_stokes
        Rproj_q1(num, grid, grid_u, grid_v, adj_der,
                fwdS.pD[current_i,:], opC_pS, BC_pS,
                fwdL.pD[current_i,:], opC_pL, BC_pL,
                fwdS.ucorrD[current_i-1,:], opC_uS, BC_uS,
                fwdL.ucorrD[current_i-1,:], opC_uL, BC_uL,
                fwdS.vcorrD[current_i-1,:], opC_vS, BC_vS,
                fwdS.vcorrD[current_i-1,:], opC_vL, BC_vL,
                fwdS.ucorrD[current_i,:], fwdL.ucorrD[current_i,:],
                fwdS.vcorrD[current_i,:], fwdL.vcorrD[current_i,:],
                Mum1_S, Mum1_L, Mvm1_S, Mvm1_L,
                u1, periodic_x, periodic_y, ϵ_adj,
                opS, phS, opL, phL,
                ns_solid_phase, ns_liquid_phase, ns_advection)
    end

    if navier_stokes
        if ns_solid_phase
            adjoint_projection_fs(num, grid, grid_u, grid_v,
                                  adj, adj.phS, RlsFS_ucorrS, RlsFS_vcorrS,
                                  AuS, BuS, AvS, BvS, AϕS,
                                  opC_pS, opC_uS, opC_vS, BC_pS,
                                  current_i, true, periodic_x, periodic_y)
        end
        if ns_liquid_phase
            adjoint_projection_fs(num, grid, grid_u, grid_v, 
                                  adj, adj.phL, RlsFS_ucorrL, RlsFS_vcorrL,
                                  AuL, BuL, AvL, BvL, AϕL,
                                  opC_pL, opC_uL, opC_vL, BC_pL,
                                  current_i, true, periodic_x, periodic_y)
        end
    end

    if levelset && free_surface
        rhs = J_u(num, grid, fwd.u[current_i,:,:], current_i-1)

        rhs .-= transpose(RuS_ls) * adj.phS.u[current_i,:]
        rhs .-= transpose(RuL_ls) * adj.phL.u[current_i,:]
        rhs .-= transpose(RvS_ls) * adj.phS.v[current_i,:]
        rhs .-= transpose(RvL_ls) * adj.phL.v[current_i,:]
        
        rhs .-= transpose(RpS_ls) * adj.phS.pD[current_i,:]
        rhs .-= transpose(RpL_ls) * adj.phL.pD[current_i,:]

        rhs .-= transpose(RucorrS_ls1) * adj.phS.ucorrD[current_i,:]
        rhs .-= transpose(RucorrL_ls1) * adj.phL.ucorrD[current_i,:]
        rhs .-= transpose(RvcorrS_ls1) * adj.phS.vcorrD[current_i,:]
        rhs .-= transpose(RvcorrL_ls1) * adj.phL.vcorrD[current_i,:]

        @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
        println(ch)
    end

    if verbose
        if current_i%show_every == 0
            try
                printstyled(color=:green, @sprintf "Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                if (heat && length(MIXED) != 0)
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                    normTS = norm(adj.phS.TD[current_i,:])
                    normTL = norm(adj.phL.TD[current_i,:])
                    normu = norm(adj.u[current_i,:])
                    print("$(@sprintf("norm(ψ) %.6e", normu))\t$(@sprintf("norm(ΘS) %.6e", normTS))\t$(@sprintf("norm(ΘL) %.6e", normTL))\n")
                end
                if navier_stokes
                    normΨ = norm(adj.u[current_i,:])
                    if ns_solid_phase
                        normu = norm(adj.phS.u[current_i,:])
                        normv = norm(adj.phS.v[current_i,:])
                        normp = norm(adj.phS.pD[current_i,:])
                        normucorr = norm(adj.phS.ucorrD[current_i,:])
                        normvcorr = norm(adj.phS.vcorrD[current_i,:])
                        print("$(@sprintf("norm(uS) %.6e", normu))\t$(@sprintf("norm(vS) %.6e", normv))\t$(@sprintf("norm(pS) %.6e", normp))\t$(@sprintf("norm(ucorrS) %.6e", normucorr))\t$(@sprintf("norm(vcorrS) %.6e", normvcorr))\t$(@sprintf("norm(ψ) %.6e", normΨ))\n")
                    end
                    if ns_liquid_phase
                        normu = norm(adj.phL.u[current_i,:])
                        normv = norm(adj.phL.v[current_i,:])
                        normp = norm(adj.phL.pD[current_i,:])
                        normucorr = norm(adj.phL.ucorrD[current_i,:])
                        normvcorr = norm(adj.phL.vcorrD[current_i,:])
                        print("$(@sprintf("norm(uL) %.6e", normu))\t$(@sprintf("norm(vL) %.6e", normv))\t$(@sprintf("norm(pL) %.6e", normp))\t$(@sprintf("norm(ucorrL) %.6e", normucorr))\t$(@sprintf("norm(vcorrL) %.6e", normvcorr))\t$(@sprintf("norm(ψ) %.6e", normΨ))\n")
                    end
                end
                println()
            catch
                @show (MIXED)
            end
        end
    end

    current_i -= 1

    while current_i > 1
        @views u0 .= fwd.u[current_i-1,:,:]
        @views u1 .= fwd.u[current_i,:,:]

        @views phS.T .= fwdS.T[current_i,:,:]
        @views phL.T .= fwdL.T[current_i,:,:]
        @views phS.TD .= fwdS.TD[current_i,:]
        @views phL.TD .= fwdL.TD[current_i,:]

        @views phS.uD .= fwdS.ucorrD[current_i-1,:]
        @views phL.uD .= fwdL.ucorrD[current_i-1,:]
        @views phS.vD .= fwdS.vcorrD[current_i-1,:]
        @views phL.vD .= fwdL.vcorrD[current_i-1,:]

        @views phS.ϕD .= fwdS.pD[current_i,:]
        @views phL.ϕD .= fwdL.pD[current_i,:]

        if levelset
            NB_indices = update_ls_data(num, grid, grid_u, grid_v, u0, κ, periodic_x, periodic_y)
            
            if advection && stefan
                update_stefan_velocity(num, grid, u0, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
                IIOE(grid, LSA, LSB, u0, V, CFL_sc, periodic_x, periodic_y)
            elseif advection && navier_stokes
                if free_surface && ns_liquid_phase
                    update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)
                end
                if free_surface && ns_solid_phase
                    update_free_surface_velocity(num, grid_u, grid_v, phS.uD, phS.vD, periodic_x, periodic_y)
                end
        
                level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
                utmp .= reshape(gmres(LSA,(LSB*vec(u0))), (ny,nx))
                S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u0, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            end
        end

        if heat
            if heat_solid_phase
                AS, BS, _ = set_heat!(dir, num, grid, opC_TS, geoS, phS, θd, BC_TS, MIXED, geoS.projection,
                                    opS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                    periodic_x, periodic_y, heat_convection)
            end
            if heat_liquid_phase
                AL, BL, _ = set_heat!(dir, num, grid, opC_TL, geoL, phL, θd, BC_TL, MIXED, geoL.projection,
                                    opL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                    periodic_x, periodic_y, heat_convection)
            end
        end

        if advection && stefan
            tmpχ_S = opC_TS.χ
            tmpχ_L = opC_TL.χ
            @views Rheat_q0(num, grid, grid_u, grid_v, adj_der,
                phS.TD, fwdS.TD[current_i+1,:], ASm1, BSm1, opC_TS, BC_TS,
                phL.TD, fwdL.TD[current_i+1,:], ALm1, BLm1, opC_TL, BC_TL,
                u1, fwd.u[current_i+1,:,:], LSAm1, LSBm1, tmpχ_S, tmpχ_L,
                CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean,
                opS, phS, opL, phL,
                heat_solid_phase, heat_liquid_phase, heat_convection)

            @views Rheat_q1(num, grid, grid_u, grid_v, adj_der,
                phS.TD, phL.TD,
                u0, u1, LSA, LSB,
                CFL_sc, periodic_x, periodic_y, ϵ_adj, λ, Vmean,
                heat_solid_phase, heat_liquid_phase)
        end
        
        if levelset && advection && stefan
            rhs = J_u(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
            rhs .-= transpose(RheatS_ls) * adj.phS.TD[current_i+1,:]
            rhs .-= transpose(RheatL_ls) * adj.phL.TD[current_i+1,:]
            rhs .-= transpose(RlsS_ls) * adj.u[current_i+1,:]
            @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
            println(ch)
            LSAm1 .= LSA
            LSBm1 .= LSB
        end
    
        if heat
            if heat_solid_phase
                rhs = J_TS(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
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
                rhs = J_TL(num, grid, fwdS.TD, fwdL.TD, u0, current_i-1)
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

        if levelset && navier_stokes
            if ns_solid_phase
                _ = set_laplacians!(grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                    opC_pS, opC_uS, opC_vS,
                    periodic_x, periodic_y, true)
                Mum1_S = opC_uS.M
                Mvm1_S = opC_vS.M
            end
            if ns_liquid_phase
                _ = set_laplacians!(grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                    opC_pL, opC_uL, opC_vL,
                    periodic_x, periodic_y, true)
                Mum1_L = opC_uL.M
                Mvm1_L = opC_vL.M
            end
            NB_indices = update_ls_data(num, grid, grid_u, grid_v, u1, κ, periodic_x, periodic_y)
        end

        if navier_stokes && free_surface
            if ns_solid_phase
                AuS, BuS, _, AvS, BvS, _, AϕS, _ = set_navier_stokes(num, grid, grid.geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                                    opC_pS, opC_uS, opC_vS, BC_pS, BC_uS, BC_vS,
                                                                    Mum1_S, Mvm1_S, iRe,
                                                                    opS, phS,
                                                                    periodic_x, periodic_y, ns_advection)
            end
            if ns_liquid_phase
                AuL, BuL, _, AvL, BvL, _, AϕL, _ = set_navier_stokes(num, grid, grid.geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                                    opC_pL, opC_uL, opC_vL, BC_pL, BC_uL, BC_vL,
                                                                    Mum1_L, Mvm1_L, iRe,
                                                                    opL, phL,
                                                                    periodic_x, periodic_y, ns_advection)
            end
        end

        if levelset && navier_stokes
            Rproj_q0(num, grid, grid_u, grid_v, adj_der,
                    opC_pS, BC_pS, BC_pL,
                    fwdS.ucorrD[current_i,:], opC_uS, BC_uS,
                    fwdL.ucorrD[current_i,:], BC_uL,
                    fwdS.vcorrD[current_i,:], opC_vS, BC_vS,
                    fwdL.vcorrD[current_i,:], BC_vL,
                    # BuSm1, BuLm1, BvSm1, BvLm1,
                    BuS, BuL, BvS, BvL,
                    opC_uS.M, opC_uL.M, opC_vS.M, opC_vL.M,
                    u1, fwd.u[current_i+1,:,:], LSAm1, LSBm1,
                    periodic_x, periodic_y, ϵ_adj,
                    ns_solid_phase, ns_liquid_phase)

            Rproj_q1(num, grid, grid_u, grid_v, adj_der,
                    phS.ϕD, opC_pS, BC_pS,
                    phL.ϕD, opC_pL, BC_pL,
                    phS.uD, opC_uS, BC_uS,
                    phL.uD , opC_uL, BC_uL,
                    phS.vD, opC_vS, BC_vS,
                    phL.vD, opC_vL, BC_vL,
                    fwdS.ucorrD[current_i,:], fwdL.ucorrD[current_i,:],
                    fwdS.vcorrD[current_i,:], fwdL.vcorrD[current_i,:],
                    Mum1_S, Mum1_L, Mvm1_S, Mvm1_L,
                    u1, periodic_x, periodic_y, ϵ_adj,
                    opS, phS, opL, phL,
                    ns_solid_phase, ns_liquid_phase, ns_advection)
        end

        if navier_stokes
            if ns_solid_phase
                adjoint_projection_fs(num, grid, grid_u, grid_v,
                                      adj, adj.phS, RlsFS_ucorrS, RlsFS_vcorrS,
                                      AuS, BuS, AvS, BvS, AϕS,
                                      opC_pS, opC_uS, opC_vS, BC_pS,
                                      current_i, false, periodic_x, periodic_y)

                BuSm1 .= BuS
                BvSm1 .= BvS
            end
            if ns_liquid_phase
                adjoint_projection_fs(num, grid, grid_u, grid_v, 
                                      adj, adj.phL, RlsFS_ucorrL, RlsFS_vcorrL,
                                      AuL, BuL, AvL, BvL, AϕL,
                                      opC_pL, opC_uL, opC_vL, BC_pL,
                                      current_i, false, periodic_x, periodic_y)
                            
                BuLm1 .= BuL
                BvLm1 .= BvL
            end
        end

        if levelset && free_surface
            rhs = J_u(num, grid, fwd.u[current_i,:,:], current_i-1)
            rhs .-= transpose(RuS_ls) * adj.phS.u[current_i,:]
            rhs .-= transpose(RuL_ls) * adj.phL.u[current_i,:]
            rhs .-= transpose(RvS_ls) * adj.phS.v[current_i,:]
            rhs .-= transpose(RvL_ls) * adj.phL.v[current_i,:]
            
            rhs .-= transpose(RpS_ls) * adj.phS.pD[current_i,:]
            rhs .-= transpose(RpL_ls) * adj.phL.pD[current_i,:]
    
            rhs .-= transpose(RucorrS_ls1) * adj.phS.ucorrD[current_i,:]
            rhs .-= transpose(RucorrL_ls1) * adj.phL.ucorrD[current_i,:]
            rhs .-= transpose(RvcorrS_ls1) * adj.phS.vcorrD[current_i,:]
            rhs .-= transpose(RvcorrL_ls1) * adj.phL.vcorrD[current_i,:]

            rhs .-= transpose(RlsFS_ls) * adj.u[current_i+1,:]
            rhs .-= transpose(RucorrS_ls0) * adj.phS.ucorrD[current_i+1,:]
            rhs .-= transpose(RucorrL_ls0) * adj.phL.ucorrD[current_i+1,:]
            rhs .-= transpose(RvcorrS_ls0) * adj.phS.vcorrD[current_i+1,:]
            rhs .-= transpose(RvcorrL_ls0) * adj.phL.vcorrD[current_i+1,:]
    
            @views @mytime (_, ch) = gmres!(adj.u[current_i,:], transpose(LSA), rhs, log=true)
            println(ch)
            LSAm1 .= LSA
            LSBm1 .= LSB
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    if (heat && length(MIXED) != 0)
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                        normTS = norm(adj.phS.TD[current_i,:])
                        normTL = norm(adj.phL.TD[current_i,:])
                        normu = norm(adj.u[current_i,:])
                        print("$(@sprintf("norm(ψ) %.6e", normu))\t$(@sprintf("norm(ΘS) %.6e", normTS))\t$(@sprintf("norm(ΘL) %.6e", normTL))\n")
                    end
                    if navier_stokes
                        normΨ = norm(adj.u[current_i,:])
                        if ns_solid_phase
                            normu = norm(adj.phS.u[current_i,:])
                            normv = norm(adj.phS.v[current_i,:])
                            normp = norm(adj.phS.pD[current_i,:])
                            normucorr = norm(adj.phS.ucorrD[current_i,:])
                            normvcorr = norm(adj.phS.vcorrD[current_i,:])
                            print("$(@sprintf("norm(uS) %.6e", normu))\t$(@sprintf("norm(vS) %.6e", normv))\t$(@sprintf("norm(pS) %.6e", normp))\t$(@sprintf("norm(ucorrS) %.6e", normucorr))\t$(@sprintf("norm(vcorrS) %.6e", normvcorr))\t$(@sprintf("norm(ψ) %.6e", normΨ))\n")
                        end
                        if ns_liquid_phase
                            normu = norm(adj.phL.u[current_i,:])
                            normv = norm(adj.phL.v[current_i,:])
                            normp = norm(adj.phL.pD[current_i,:])
                            normucorr = norm(adj.phL.ucorrD[current_i,:])
                            normvcorr = norm(adj.phL.vcorrD[current_i,:])
                            print("$(@sprintf("norm(uL) %.6e", normu))\t$(@sprintf("norm(vL) %.6e", normv))\t$(@sprintf("norm(pL) %.6e", normp))\t$(@sprintf("norm(ucorrL) %.6e", normucorr))\t$(@sprintf("norm(vcorrL) %.6e", normvcorr))\t$(@sprintf("norm(ψ) %.6e", normΨ))\n")
                        end
                    end
                    println()
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