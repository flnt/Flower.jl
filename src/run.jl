function run_forward(
    num, grid, grid_u, grid_v, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(),
    BC_TL = Boundaries(),
    BC_pS = Boundaries(),
    BC_pL = Boundaries(),
    BC_uS = Boundaries(),
    BC_uL = Boundaries(),
    BC_vS = Boundaries(),
    BC_vL = Boundaries(),
    BC_u = Boundaries(),
    time_scheme = CN,
    ls_scheme = weno5,
    stefan = false,
    advection = false,
    auto_reinit = false,
    heat = false,
    heat_convection = false,
    ns_advection = false,
    navier_stokes = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
    free_surface = false,
    hill = false,
    Vmean = false,
    levelset = true,
    speed = 0,
    analytical = false,
    verbose = false,
    show_every = 100,
    save_length = false,
    save_radius = false,
    adaptative_t = false,
    breakup = false,
    Ra = 0,
    λ = 1,
    )

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
    end

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations,
            current_i, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso = num
    @unpack opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL = op
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid

    iRe = 1.0 / Re

    local NB_indices;

    local Cum1S = fzeros(grid_u)
    local Cum1L = fzeros(grid_u)
    local Cvm1S = fzeros(grid_v)
    local Cvm1L = fzeros(grid_v)

    local Mm1_S
    local Mum1_S
    local Mvm1_S

    local Mm1_L
    local Mum1_L
    local Mvm1_L

    θ_out = zeros(grid, 4)
    utmp = copy(u)

    if periodic_x
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
    else
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
    end

    if periodic_y
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
    else
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
    end

    if navier_stokes
        if free_surface
            bc_type_u = neu
            bc_type_p = fs
        else
            bc_type_u = dir
            bc_type_p = neu
        end
    end

    if levelset
        NB_indices = update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)

        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, ind.MIXED)
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(ind.MIXED))
            for i in eachindex(ind.MIXED)
                a[i] = geoL.projection[ind.MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        ind.MIXED = [CartesianIndex(-1,-1)]
        grid_u.ind.MIXED = [CartesianIndex(-1,-1)]
        grid_v.ind.MIXED = [CartesianIndex(-1,-1)]
    end

    CFL_sc = τ / Δ^2
    IIOE_normal!(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
    rhs_LS = fzeros(grid)

    if save_length
        fwd.length[1] = arc_length2(geoS.projection, ind.MIXED)
        fwd.κ[1,:,:] .= κ
    end

    kill_dead_cells!(phS.T, grid, geoS)
    kill_dead_cells!(phL.T, grid, geoL)
    
    @views fwd.u[1,:,:] .= u
    @views fwd.ux[1,:,:] .= grid_u.u
    @views fwd.uy[1,:,:] .= grid_v.u
    @views fwd.T[1,:,:] .= phL.T.*geoL.cap[:,:,5] .+ phS.T[:,:].*geoS.cap[:,:,5]
    @views fwdL.T[1,:,:] .= phL.T
    @views fwdS.T[1,:,:] .= phS.T
    @views fwdS.p[1,:,:] .= phS.p
    @views fwdL.p[1,:,:] .= phL.p
    @views fwdS.u[1,:,:] .= phS.u
    @views fwdS.v[1,:,:] .= phS.v
    @views fwdL.u[1,:,:] .= phL.u
    @views fwdL.v[1,:,:] .= phL.v

    veci(phS.TD,grid,1) .= vec(phS.T)
    veci(phS.TD,grid,2) .= θd
    init_borders!(vec3(phS.TD,grid), BC_TS, θd)
    @views fwdS.TD[1,:] .= phS.TD

    veci(phL.TD,grid,1) .= vec(phL.T)
    veci(phL.TD,grid,2) .= θd
    init_borders!(vec3(phL.TD,grid), BC_TL, θd)
    @views fwdL.TD[1,:] .= phL.TD

    veci(phS.uD,grid_u,1) .= vec(phS.u)
    veci(phS.uD,grid_u,2) .= num.uD
    vec3(phS.uD,grid_u) .= num.u_inf
    veci(phL.uD,grid_u,1) .= vec(phL.u)
    veci(phL.uD,grid_u,2) .= num.uD
    vec3(phL.uD,grid_u) .= num.u_inf
    veci(phS.ucorrD,grid_u,1) .= vec(phS.u)
    veci(phS.ucorrD,grid_u,2) .= num.uD
    vec3(phS.ucorrD,grid_u) .= num.u_inf
    veci(phL.ucorrD,grid_u,1) .= vec(phL.u)
    veci(phL.ucorrD,grid_u,2) .= num.uD
    vec3(phL.ucorrD,grid_u) .= num.u_inf
    @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

    veci(phS.vD,grid_v,1) .= vec(phS.v)
    veci(phS.vD,grid_v,2) .= num.vD
    vec3(phS.vD,grid_v) .= num.v_inf
    veci(phL.vD,grid_v,1) .= vec(phL.v)
    veci(phL.vD,grid_v,2) .= num.vD
    vec3(phL.vD,grid_v) .= num.v_inf
    veci(phS.vcorrD,grid_v,1) .= vec(phS.v)
    veci(phS.vcorrD,grid_v,2) .= num.vD
    vec3(phS.vcorrD,grid_v) .= num.v_inf
    veci(phL.vcorrD,grid_v,1) .= vec(phL.v)
    veci(phL.vcorrD,grid_v,2) .= num.vD
    vec3(phL.vcorrD,grid_v) .= num.v_inf
    @views fwdS.vcorrD[1,:,:] .= phS.vcorrD
    @views fwdL.vcorrD[1,:,:] .= phL.vcorrD

    @views fwdL.pD[1,:] .= phL.pD
    @views fwdS.pD[1,:] .= phS.pD

    if is_FE(time_scheme) || is_CN(time_scheme)
        update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y, false)

        if navier_stokes || heat
            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                opC_pS, opC_uS, opC_vS, periodic_x, periodic_y
            )
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                opC_pL, opC_uL, opC_vL, periodic_x, periodic_y
            )
        end
            
        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    V0S = volume(grid.geoS)
    V0L = volume(grid.geoL)

    current_t = 0.

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TS, geoS, phS, θd, BC_TS, ind.MIXED, geoS.projection,
                                        opS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                        periodic_x, periodic_y, heat_convection)
                mul!(rhs, B, phS.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                bicgstabl!(phS.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                # @mytime (_, ch) = bicgstabl!(phS.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                # println(ch)
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TL, geoL, phL, θd, BC_TL, ind.MIXED, geoL.projection,
                                        opL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                        periodic_x, periodic_y, heat_convection)
                mul!(rhs, B, phL.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                bicgstabl!(phL.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                # @mytime (_, ch) = bicgstabl!(phL.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                # println(ch)
                phL.T .= reshape(veci(phL.TD,grid,1), grid)
            end
        end

        if advection && stefan
            update_stefan_velocity(num, grid, u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
        end

        if advection && free_surface
            update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)
        end

        if verbose && adaptative_t
            println("τ = $τ")
        end

        if advection
            CFL_sc = τ / Δ^2
            if stefan
                IIOE_normal!(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
                u .= reshape(gmres(LSA, (LSB * vec(u))), grid)
                # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
            elseif free_surface
                rhs_LS .= 0.0
                IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, τ, periodic_x, periodic_y)
                BC_LS!(grid, LSA, LSB, rhs_LS, BC_u, num.n_ext_cl)
                utmp .= reshape(gmres(LSA, (LSB * vec(u))) .+ rhs_LS, grid)

                rhs_LS .= 0.0
                S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u, θ_out, τ, periodic_x, periodic_y)
                BC_LS!(grid, LSA, LSB, rhs_LS, BC_u, num.n_ext_cl)
                u .= reshape(gmres(LSA, (LSB * vec(u))) .+ rhs_LS, grid)
            else
                @error ("Set either stefan or free_surface to true in order to advect the levelset")
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*τ);
            elseif nb_reinit > 0
                if auto_reinit
                    if rg(grid, periodic_x, periodic_y) >= δreinit
                        RK2_reinit!(ls_scheme, grid, ind, u, nb_reinit, periodic_x, periodic_y, BC_u)
                    end
                else
                    if current_i%num.reinit_every == 0
                        RK2_reinit!(ls_scheme, grid, ind, u, nb_reinit, periodic_x, periodic_y, BC_u)
                    end
                end
            end
            # numerical breakup
            if free_surface && breakup
                count = breakup(u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
                if count > 0
                    RK2_reinit!(ls_scheme, grid, ind, u, nb_reinit, periodic_x, periodic_y, BC_u)
                end
            end
        end

        if verbose
            if (current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                if heat && length(ind.MIXED) != 0
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[ind.MIXED]) findmax(V[ind.MIXED])[1] findmin(V[ind.MIXED])[1])
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[ind.MIXED]) findmax(κ[ind.MIXED])[1] findmin(κ[ind.MIXED])[1])
                elseif advection && (free_surface || stefan) && length(ind.MIXED) != 0
                    V_mean = mean([mean(grid_u.V[ind.MIXED]), mean(grid_v.V[ind.MIXED])])
                    V_max = max(findmax(grid_u.V[ind.MIXED])[1], findmax(grid_v.V[ind.MIXED])[1])
                    V_min = min(findmin(grid_u.V[ind.MIXED])[1], findmin(grid_v.V[ind.MIXED])[1])
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[ind.MIXED]) findmax(κ[ind.MIXED])[1] findmin(κ[ind.MIXED])[1])
                end
                if navier_stokes
                    if ns_solid_phase
                        normuS = norm(phS.u)
                        normvS = norm(phS.v)
                        normpS = norm(phS.p.*τ)
                        print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                    end
                    if ns_liquid_phase
                        normuL = norm(phL.u)
                        normvL = norm(phL.v)
                        normpL = norm(phL.p.*τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                end
            end
        end


        if levelset && (advection || current_i<2)
            NB_indices = update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)

            grid.geoL.fresh .= false
            grid.geoS.fresh .= false
            grid_u.geoL.fresh .= false
            grid_u.geoS.fresh .= false
            grid_v.geoL.fresh .= false
            grid_v.geoS.fresh .= false

            get_fresh_cells!(grid, grid.geoS, Mm1_S, grid.ind.all_indices)
            get_fresh_cells!(grid, grid.geoL, Mm1_L, grid.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.geoS, Mum1_S, grid_u.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.geoL, Mum1_L, grid_u.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.geoS, Mvm1_S, grid_v.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.geoL, Mvm1_L, grid_v.ind.all_indices)

            FRESH_L_u = findall(grid_u.geoL.fresh)
            FRESH_S_u = findall(grid_u.geoS.fresh)
            FRESH_L_v = findall(grid_v.geoL.fresh)
            FRESH_S_v = findall(grid_v.geoS.fresh)

            if navier_stokes
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,1), veci(phS.uD,grid_u,1), grid_u.geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,1), veci(phS.vD,grid_v,1), grid_v.geoS.projection, FRESH_S_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,2), veci(phS.uD,grid_u,1), grid_u.geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,2), veci(phS.vD,grid_v,1), grid_v.geoS.projection, FRESH_S_v, periodic_x, periodic_y)

                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1), grid_u.geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1), grid_v.geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1), grid_u.geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1), grid_v.geoL.projection, FRESH_L_v, periodic_x, periodic_y)
            end

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, ind.MIXED)
                end
                if hill
                    a = zeros(length(ind.MIXED))
                    for i in eachindex(ind.MIXED)
                        a[i] = geoL.projection[ind.MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    fwd.length[snap] = arc_length2(geoS.projection, ind.MIXED)
                    fwd.κ[snap,:,:] .= κ
                end
            end
        end

        if navier_stokes
            if !advection
                no_slip_condition!(grid, grid_u, grid_v)
                grid_u.V .= imfilter(grid_u.V, Kernel.gaussian(2))
                grid_v.V .= imfilter(grid_v.V, Kernel.gaussian(2))
                # grid_u.V .= Δ / (1 * τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, bc_type_u, bc_type_p,
                    num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, phS,
                    BC_uS, BC_vS, BC_pS,
                    opC_pS, opC_uS, opC_vS, opS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, current_i
                )
            end
            if ns_liquid_phase
                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, bc_type_u, bc_type_p,
                    num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                    BC_uL, BC_vL, BC_pL,
                    opC_pL, opC_uL, opC_vL, opL,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, advection, current_i
                )
                # linear_advection!(
                #     num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                #     BC_uL, BC_vL, opL
                # )
            end
        end

        current_t += τ
        if iszero(current_i%save_every) || current_i==max_iterations
            snap = current_i÷save_every+1
            if current_i==max_iterations
                snap = size(fwd.T,1)
            end
            fwd.t[snap] = current_t
            @views fwd.V[snap,:,:] .= V
            @views fwd.u[snap,:,:] .= u
            @views fwd.ux[snap,:,:] .= grid_u.u
            @views fwd.uy[snap,:,:] .= grid_v.u

            if heat_solid_phase && heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T.*geoL.cap[:,:,5] .+ phS.T.*geoS.cap[:,:,5]
                @views fwdL.T[snap,:,:] .= phL.T
                @views fwdL.TD[snap,:] .= phL.TD
                @views fwdS.T[snap,:,:] .= phS.T
                @views fwdS.TD[snap,:] .= phS.TD
            elseif heat_solid_phase
                @views fwd.T[snap,:,:] .= phS.T
                @views fwdS.T[snap,:,:] .= phS.T
                @views fwdS.TD[snap,:] .= phS.TD
            elseif heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T
                @views fwdL.T[snap,:,:] .= phL.T
                @views fwdL.TD[snap,:] .= phL.TD
            end

            if ns_solid_phase && ns_liquid_phase
                @views fwdS.p[snap,:,:] .= phS.p
                @views fwdL.p[snap,:,:] .= phL.p
                @views fwdL.pD[snap,:] .= phL.pD
                @views fwdS.pD[snap,:] .= phS.pD
                @views fwdS.ϕ[snap,:,:] .= phS.ϕ
                @views fwdL.ϕ[snap,:,:] .= phL.ϕ
                @views fwdS.u[snap,:,:] .= phS.u
                @views fwdL.u[snap,:,:] .= phL.u
                @views fwdS.v[snap,:,:] .= phS.v
                @views fwdL.v[snap,:,:] .= phL.v
                @views fwdS.ucorr[snap,:,:] .= phS.ucorr
                @views fwdL.ucorr[snap,:,:] .= phL.ucorr
                @views fwdS.vcorr[snap,:,:] .= phS.vcorr
                @views fwdL.vcorr[snap,:,:] .= phL.vcorr
                @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
                @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
                @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
                @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, fwdL; step=snap)
            elseif ns_solid_phase
                @views fwdS.p[snap,:,:] .= phS.p
                @views fwdS.pD[snap,:] .= phS.pD
                @views fwdS.ϕ[snap,:,:] .= phS.ϕ
                @views fwdS.u[snap,:,:] .= phS.u
                @views fwdS.v[snap,:,:] .= phS.v
                @views fwdS.ucorr[snap,:,:] .= phS.ucorr
                @views fwdS.vcorr[snap,:,:] .= phS.vcorr
                @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
                @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
            elseif ns_liquid_phase
                @views fwdL.p[snap,:,:] .= phL.p
                @views fwdL.pD[snap,:] .= phL.pD
                @views fwdL.ϕ[snap,:,:] .= phL.ϕ
                @views fwdL.u[snap,:,:] .= phL.u
                @views fwdL.v[snap,:,:] .= phL.v
                @views fwdL.ucorr[snap,:,:] .= phL.ucorr
                @views fwdL.vcorr[snap,:,:] .= phL.vcorr
                @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
                @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, fwdL; step=snap)
            end
            if advection
                fwdS.Vratio[snap] = volume(grid.geoS) / V0S
                fwdL.Vratio[snap] = volume(grid.geoL) / V0L
            end
        end

        if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
            norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
            println(@sprintf "\n CRASHED after %d iterations \n" current_i)
            return ind.MIXED, ind.SOLID, ind.LIQUID
        end

        current_i += 1

        if adaptative_t
            τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))
        end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if heat || advection
                print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[ind.MIXED]) findmax(V[ind.MIXED])[1] findmin(V[ind.MIXED])[1] std(V[ind.MIXED]))
                print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[ind.MIXED]) findmax(κ[ind.MIXED])[1] findmin(κ[ind.MIXED])[1] std(κ[ind.MIXED]))
            end
            if navier_stokes
                if ns_solid_phase
                    normuS = norm(phS.u)
                    normvS = norm(phS.v)
                    normpS = norm(phS.p.*τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                end
            end
            print("\n\n")
        catch
            @show (length(ind.MIXED))
        end
    end

    if levelset
        if save_radius || hill
            return ind.MIXED, ind.SOLID, ind.LIQUID, radius
        end
        return ind.MIXED, ind.SOLID, ind.LIQUID
    else
        return ind.MIXED
    end
end

function run_backward(num, grid, opS, opL, fwd, adj;
    periodic_x = false,
    periodic_y = false,
    BC_TS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_TL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    stefan = false,
    advection = false,
    heat = false,
    heat_convection = false,
    liquid_phase = false,
    solid_phase = false,
    hill = false,
    Vmean = false,
    levelset = true,
    speed = 0,
    analytical = false,
    verbose = false,
    show_every = 100,
    save_length = false,
    save_radius = false
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, max_iterations, current_i, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack nx, ny, ind, faces, geoS, geoL, mid_point = grid
    @unpack all_indices, inside, b_left, b_bottom, b_right, b_top = ind
    @unpack usave, TSsave, TLsave, Tsave, Vsave, κsave = fwd
    @unpack iso, u, TS, TL, DTS, DTL, κ, V = adj

    local MIXED; local SOLID; local LIQUID;
    local WAS_SOLID; local WAS_LIQUID;
    local NB_indices_base; local NB_indices;
    local FRESH_L; local FRESH_S;

    if periodic_x
        BC_u.left.ind = idx.periodicL;
        BC_u.right.ind = idx.periodicR;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = b_left;
        BC_u.right.ind = b_right;
    end

    if periodic_y
        BC_u.bottom.ind = idx.periodicB;
        BC_u.top.ind = idx.periodicT;
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = b_bottom;
        BC_u.top.ind = b_top;
    end

    current_i = max_iterations + 1

    if levelset
        marching_squares!(num, grid, u, periodic_x, periodic_y)

        NB_indices_base = get_NB_width_indices_base(NB)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
        NB_indices = get_NB_width(MIXED, NB_indices_base)

        get_iterface_location!(grid, MIXED)
        get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
    end

    HS = zeros(ny, nx)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HS[II] = distance(mid_point[II], geoS.centroid[II], dx[II], dy[II])
    end
    
    HL = zeros(ny, nx)
    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
        HL[II] = distance(mid_point[II], geoL.centroid[II], dx[II], dy[II])
    end

    DTS .= θd
    DTL .= θd
    bcS = similar(DTS)
    bcL = similar(DTL)
    apply_curvature(bcT, DTS, num, grid, all_indices)
    apply_curvature(bcT, DTL, num, grid, all_indices)
    if aniso
        apply_anisotropy(bcS, DTS, MIXED, num, grid, geoS.projection)
        apply_anisotropy(bcL, DTL, MIXED, num, grid, geoS.projection)
    end
    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)
    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

    laplacian!(dir, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])

    while current_i > 1

        if heat
            opL.CUTT .= zeros(nx*ny)
            opS.CUTT .= zeros(nx*ny)

            try
                if solid_phase
                    HS .= 0.
                    for II in vcat(ind.b_left[1], ind.b_bottom[1], ind.b_right[1], ind.b_top[1])
                        HS[II] = distance(mid_point[II], geoS.centroid[II], dx[II], dy[II])
                    end

                    DTS .= θd
                    apply_curvature(bcT, DTS, num, grid, all_indices)
                    if aniso
                        apply_anisotropy(bcS, DTS, MIXED, num, grid, geoS.projection)
                    end
                    bcSx, bcSy = set_bc_bnds(dir, bcS, HS, BC_TS)

                    laplacian!(dir, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoS, opS)
                    TS .= reshape(gmres(opS.A,(opS.B*vec(TS) + 2.0*τ*opS.CUTT)), (ny,nx))
                end
                if liquid_phase
                    HL .= 0.
                    for II in vcat(b_left[1], b_bottom[1], b_right[1], b_top[1])
                        HL[II] = distance(mid_point[II], geoL.centroid[II], dx[II], dy[II])
                    end

                    DTL .= θd
                    apply_curvature(bcT, DTL, num, grid, all_indices)
                    if aniso
                        apply_anisotropy(bcL, DTL, MIXED, num, grid, geoS.projection)
                    end
                    bcLx, bcLy = set_bc_bnds(dir, bcL, HL, BC_TL)

                    laplacian!(dir, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoL, opL)
                    TL .= reshape(gmres(opL.A,(opL.B*vec(TL) + 2.0*τ*opL.CUTT)), (ny,nx))
                end
            catch
                @error ("Unphysical temperature field, iteration $current_i")
                break
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                    print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                catch
                    @show (MIXED)
                end
            end
        end

        if levelset
            marching_squares!(num, grid, u, periodic_x, periodic_y)

            WAS_LIQUID = copy(LIQUID)
            WAS_SOLID = copy(SOLID)

            MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)
            NB_indices = Flower.get_NB_width(MIXED, NB_indices_base)

            get_iterface_location!(grid, MIXED)

            FRESH_L = intersect(MIXED, WAS_SOLID)
            FRESH_S = intersect(MIXED, WAS_LIQUID)

            init_fresh_cells!(grid, TS, geoS.projection, FRESH_S, periodic_x, periodic_y)
            init_fresh_cells!(grid, TL, geoL.projection, FRESH_L, periodic_x, periodic_y)

            get_curvature(num, grid, u, MIXED, periodic_x, periodic_y)
        end

        current_i -= 1
        κ .= κsave[current_i,:,:]
        u .= usave[current_i,:,:]
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
            print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            print("\n \n")
        catch
            @show (length(MIXED))
        end
    end
end
