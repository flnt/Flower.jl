function run_forward(num, grid, grid_u, grid_v,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, 
    phS, phL, fwd;
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
    BC_pS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_pL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_uS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_uL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_vS = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_vL = Boundaries(
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
    Ra = 0,
    λ = 1,
    )

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
    end

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack Tall, usave, uusave, uvsave, TSsave, TLsave, TDSsave, TDLsave, Tsave, psave, ϕsave, Uxsave, Uysave, Uxcorrsave, Uycorrsave, Vsave, κsave, lengthsave, tv, Cd, Cl = fwd

    iRe = 1.0 / Re

    local MIXED; local SOLID; local LIQUID;
    local MIXED_vel_ext; local SOLID_vel_ext; local LIQUID_vel_ext;
    local MIXED_u_vel_ext; local SOLID_u_vel_ext; local LIQUID_u_vel_ext;
    local MIXED_v_vel_ext; local SOLID_v_vel_ext; local LIQUID_v_vel_ext;
    local indices_u_vel_ext;
    local MIXED_u; local SOLID_u; local LIQUID_u;
    local MIXED_v; local SOLID_v; local LIQUID_v;

    local Cum1S = zeros(grid_u.nx*grid_u.ny)
    local Cum1L = zeros(grid_u.nx*grid_u.ny)
    local Cvm1S = zeros(grid_v.nx*grid_v.ny)
    local Cvm1L = zeros(grid_v.nx*grid_v.ny)

    local Lum1_S
    local bc_Lum1_S
    local Lvm1_S
    local bc_Lvm1_S
    local Mm1_S
    local Mum1_S
    local Mvm1_S

    local Lum1_L
    local bc_Lum1_L
    local Lvm1_L
    local bc_Lvm1_L
    local Mm1_L
    local Mum1_L
    local Mvm1_L

    θ_out = zeros(ny, nx, 4)
    utmp = copy(u)

    HTS = zeros(ny,nx)
    HTL = zeros(ny,nx)
    phS.DT .= phS.T
    phL.DT .= phL.T
    bcTS = copy(phS.DT)
    bcTL = copy(phL.DT)

    HuS = zeros(grid_u.ny,grid_u.nx)
    HuL = zeros(grid_u.ny,grid_u.nx)

    HvS = zeros(grid_v.ny,grid_v.nx)
    HvL = zeros(grid_v.ny,grid_v.nx)

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

    if levelset
        grid.mid_point .= [Point(0.0, 0.0)]
        grid_u.mid_point .= [Point(0.0, 0.0)]
        grid_v.mid_point .= [Point(0.0, 0.0)]
        
        marching_squares!(num, grid, u, periodic_x, periodic_y)
        interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
        uusave[1,:,:] .= grid_u.u
        uvsave[1,:,:] .= grid_v.u
        marching_squares!(num, grid_u, grid_u.u, periodic_x, periodic_y)
        marching_squares!(num, grid_v, grid_v.u, periodic_x, periodic_y)

        MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.all_indices, nx, ny, periodic_x, periodic_y)
        MIXED_u_vel_ext, SOLID_u_vel_ext, LIQUID_u_vel_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
        MIXED_v_vel_ext, SOLID_v_vel_ext, LIQUID_v_vel_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
        MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

        kill_dead_cells!(phS.T, grid, geoS)
        kill_dead_cells!(phL.T, grid, geoL)

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

        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, MIXED)
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(MIXED))
            for i in eachindex(MIXED)
                a[i] = geoL.projection[MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
        MIXED_u = [CartesianIndex(-1,-1)]
        MIXED_v = [CartesianIndex(-1,-1)]
    end

    # Heat
    set_heat!(num, grid, geoS, geoS.projection,
            opS, phS, HTS, bcTS, HuS, HvS,
            BC_TS, BC_uS, BC_vS,
            MIXED, LIQUID, heat_convection)
    
    set_heat!(num, grid, geoL, geoS.projection,
            opL, phL, HTL, bcTL, HuL, HvL,
            BC_TL, BC_uL, BC_vL,
            MIXED, SOLID, heat_convection)

    if ns_advection
        Cum1S .= opS.Cu * vec(phS.u) .+ opS.CUTCu
        Cvm1S .= opS.Cv * vec(phS.v) .+ opS.CUTCv
        Cum1L .= opL.Cu * vec(phL.u) .+ opL.CUTCu
        Cvm1L .= opL.Cv * vec(phL.v) .+ opL.CUTCv
    end

    CFL_sc = τ / Δ^2
    IIOE(LSA, LSB, u, V, ind.inside, CFL_sc, ny)

    if save_length
        lengthsave[1] = arc_length2(geoS.projection, MIXED)
        κsave[1,:,:] .= κ
    end

    usave[1,:,:] .= u
    Tsave[1,:,:] .= phL.T.*geoL.cap[:,:,5] .+ phS.T[:,:].*geoS.cap[:,:,5]
    TLsave[1,:,:] .= phL.T
    TSsave[1,:,:] .= phS.T
    psave[1,:,:] .= phL.p .+ phS.p
    Uxsave[1,:,:] .= phL.u .+ phS.u
    Uysave[1,:,:] .= phL.v .+ phS.v

    tmpDS = copy(phS.T)
    tmpDS[2:end-1,2:end-1] .= θd
    init_borders!(tmpDS, BC_TS, θd)
    tmpDL = copy(phL.T)
    tmpDL[2:end-1,2:end-1] .= θd
    init_borders!(tmpDL, BC_TL, θd)

    veci(phS.TD,grid,1) .= vec(phS.T)
    veci(phL.TD,grid,1) .= vec(phL.T)
    veci(phS.TD,grid,2) .= vec(tmpDS)
    veci(phL.TD,grid,2) .= vec(tmpDL)
    @views TDLsave[1,:] .= phL.TD
    @views TDSsave[1,:] .= phS.TD

    tmpu = ones(grid_u.ny, grid_u.nx) .* num.u_inf
    tmpu[2:end-1,2:end-1] .= 0.0
    veci(phS.uD,grid_u,1) .= vec(phS.u)
    veci(phS.uD,grid_u,2) .= vec(tmpu)
    veci(phL.uD,grid_u,1) .= vec(phL.u)
    veci(phL.uD,grid_u,2) .= vec(tmpu)

    tmpv = ones(grid_v.ny, grid_v.nx) .* num.v_inf
    tmpv[2:end-1,2:end-1] .= 0.0
    veci(phS.vD,grid_v,1) .= vec(phS.v)
    veci(phS.vD,grid_v,2) .= vec(tmpv)
    veci(phL.vD,grid_v,1) .= vec(phL.v)
    veci(phL.vD,grid_v,2) .= vec(tmpv)

    if free_surface
        Lpm1_S, bc_Lpm1_S, Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Lpm1_fs_S, bc_Lpm1_fs_S, Lum1_fs_S, bc_Lum1_fs_S, Lvm1_fs_S, bc_Lvm1_fs_S = set_laplacians!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                  opC_pS, opC_uS, opC_vS, periodic_x, periodic_y, true)
        Lpm1_L, bc_Lpm1_L, Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Lpm1_fs_L, bc_Lpm1_fs_L, Lum1_fs_L, bc_Lum1_fs_L, Lvm1_fs_L, bc_Lvm1_fs_L = set_laplacians!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                  opC_pL, opC_uL, opC_vL, periodic_x, periodic_y, true)
        
        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)
    else
        Lpm1_S, bc_Lpm1_S, Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S = set_laplacians!(grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS,
                                                  opC_pS, opC_uS, opC_vS, periodic_x, periodic_y)
        Lpm1_L, bc_Lpm1_L, Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L = set_laplacians!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                  opC_pL, opC_uL, opC_vL, periodic_x, periodic_y)
        
        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)
    end

    current_t = 0.

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            if heat_solid_phase
                veci(phS.TD,grid,1) .= vec(phS.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, MIXED, geoS.projection,
                                        periodic_x, periodic_y)
                mul!(rhs, B, phS.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                @mytime (_, ch) = bicgstabl!(phS.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
                phS.T .= reshape(veci(phS.TD,grid,1), (ny, nx))
            end
            if heat_liquid_phase
                veci(phL.TD,grid,1) .= vec(phL.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, MIXED, geoL.projection,
                                        periodic_x, periodic_y)
                mul!(rhs, B, phL.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                @mytime (_, ch) = bicgstabl!(phL.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
                phL.T .= reshape(veci(phL.TD,grid,1), (ny, nx))
            end
        end

        if stefan
            Stefan_velocity!(num, grid, V, phS.T, phL.T, MIXED, periodic_x, periodic_y)
            V[MIXED] .*= 1. ./ λ
            if Vmean
                a = mean(V[MIXED])
                V[MIXED] .= a
            end
            _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
            _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
            _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
            indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
            velocity_extension!(grid, u, V, indices_vel_ext, NB, periodic_x, periodic_y)
        end

        if free_surface
            grid_u.V .= reshape(veci(phL.uD,grid_u,2), (grid_u.ny, grid_u.nx))
            grid_v.V .= reshape(veci(phL.vD,grid_v,2), (grid_v.ny, grid_v.nx))

            _MIXED_L_u_vel_ext = intersect(findall(grid_u.geoL.emptied),
                                           MIXED_u_vel_ext)
            _MIXED_S_u_vel_ext = intersect(findall(grid_u.geoS.emptied),
                                           MIXED_u_vel_ext)
            _MIXED_u_vel_ext = vcat(_MIXED_L_u_vel_ext, _MIXED_S_u_vel_ext)
            indices_u_vel_ext = vcat(SOLID_u_vel_ext, _MIXED_u_vel_ext, LIQUID_u_vel_ext)

            _MIXED_L_v_vel_ext = intersect(findall(grid_v.geoL.emptied),
                                           MIXED_v_vel_ext)
            _MIXED_S_v_vel_ext = intersect(findall(grid_v.geoS.emptied),
                                           MIXED_v_vel_ext)
            _MIXED_v_vel_ext = vcat(_MIXED_L_v_vel_ext, _MIXED_S_v_vel_ext)
            indices_v_vel_ext = vcat(SOLID_v_vel_ext, _MIXED_v_vel_ext, LIQUID_v_vel_ext)

            velocity_extension!(grid_u, grid_u.u, grid_u.V, indices_u_vel_ext, NB, periodic_x, periodic_y)
            velocity_extension!(grid_v, grid_v.u, grid_v.V, indices_v_vel_ext, NB, periodic_x, periodic_y)
        end

        if verbose && adaptative_t
            println("τ = $τ")
        end

        if advection
            CFL_sc = τ / Δ^2
            if stefan
                IIOE(grid, LSA, LSB, u, V, CFL_sc, periodic_x, periodic_y)
                try
                    u .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                    # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
                catch
                    @error ("Inadequate level set function, iteration $current_i")
                    break
                end
            elseif free_surface
                level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, MIXED, τ, periodic_x, periodic_y)
                try
                    utmp .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                catch
                    @error ("Inadequate level set function, iteration $current_i")
                    break
                end
                S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u, θ_out, MIXED, τ, periodic_x, periodic_y)
                try
                    u .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                catch
                    @error ("Inadequate level set function, iteration $current_i")
                    break
                end
            else
                @error ("Set either stefan or free_surface to true in order to advect the levelset")
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*τ);
            elseif nb_reinit > 0
                if current_i%num.reinit_every == 0
                    FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
                end
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    if (heat && length(MIXED) != 0)
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                    elseif free_surface
                        V_mean = mean([mean(grid_u.V[MIXED]), mean(grid_v.V[MIXED])])
                        V_max = max(findmax(grid_u.V[MIXED])[1], findmax(grid_v.V[MIXED])[1])
                        V_min = min(findmin(grid_u.V[MIXED])[1], findmin(grid_v.V[MIXED])[1])
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" V_mean V_max V_min)
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
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
                catch
                    @show (MIXED)
                end
            end
        end


        if levelset && (advection || current_i<2)
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
            MIXED_u_vel_ext, SOLID_u_vel_ext, LIQUID_u_vel_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
            MIXED_v_vel_ext, SOLID_v_vel_ext, LIQUID_v_vel_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
            MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
            MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
            MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

            kill_dead_cells!(phS.T, grid, geoS)
            kill_dead_cells!(phL.T, grid, geoL)

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

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, MIXED)
                end
                if hill
                    a = zeros(length(MIXED))
                    for i in eachindex(MIXED)
                        a[i] = geoL.projection[MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    lengthsave[snap] = arc_length2(geoS.projection, MIXED)
                    κsave[snap,:,:] .= κ
                end
            end
        end

        if navier_stokes
            if !free_surface
                no_slip_condition!(grid, grid_u, grid_v)
                grid_u.V .= imfilter(grid_u.V, Kernel.gaussian(2))
                grid_v.V .= imfilter(grid_v.V, Kernel.gaussian(2))
                # grid_u.V .= Δ / (1 * τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                if free_surface
                    Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mm1_S, Mum1_S, Mvm1_S = projection_fs!(num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, phS,
                                                                                          BC_uS, BC_vS, BC_pS,
                                                                                          opC_pS, opC_uS, opC_vS,
                                                                                          Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S,
                                                                                          SOLID, MIXED, periodic_x, periodic_y, current_i)
                else
                    Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S = projection_no_slip!(num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, phS,
                                                                                               BC_uS, BC_vS, BC_pS,
                                                                                               opC_pS, opC_uS, opC_vS,
                                                                                               Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S,
                                                                                               SOLID, MIXED, periodic_x, periodic_y)
                end
                
            end
            if ns_liquid_phase
                if free_surface
                    Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mm1_L, Mum1_L, Mvm1_L = projection_fs!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                                                                                          BC_uL, BC_vL, BC_pL,
                                                                                          opC_pL, opC_uL, opC_vL,
                                                                                          Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                                                                                          LIQUID, MIXED, periodic_x, periodic_y, current_i)
                else
                   Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L = projection_no_slip!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                                                                                              BC_uL, BC_vL, BC_pL,
                                                                                              opC_pL, opC_uL, opC_vL,
                                                                                              Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                                                                                              LIQUID, MIXED, periodic_x, periodic_y)
                end
            end
        end

        current_t += τ
        if iszero(current_i%save_every) || current_i==max_iterations
            snap = current_i÷save_every+1
            if current_i==max_iterations
                snap = size(Tsave,1)
            end
            tv[snap] = current_t
            @views Vsave[snap,:,:] .= V
            @views usave[snap,:,:] .= u
            @views uusave[snap,:,:] .= grid_u.u
            @views uvsave[snap,:,:] .= grid_v.u

            if heat_solid_phase && heat_liquid_phase
                @views Tsave[snap,:,:] .= phL.T.*geoL.cap[:,:,5] .+ phS.T[:,:].*geoS.cap[:,:,5]
                @views TLsave[snap,:,:] .= phL.T
                @views TDLsave[snap,:] .= phL.TD
                @views TSsave[snap,:,:] .= phS.T
                @views TDSsave[snap,:] .= phS.TD
            elseif heat_solid_phase
                @views Tsave[snap,:,:] .= phS.T
                @views TSsave[snap,:,:] .= phS.T
                @views TDSsave[snap,:] .= phS.TD
            elseif heat_liquid_phase
                @views Tsave[snap,:,:] .= phL.T
                @views TLsave[snap,:,:] .= phL.T
                @views TDLsave[snap,:] .= phL.TD
            end

            if ns_solid_phase && ns_liquid_phase
                @views psave[snap,:,:] .= phL.p.*geoL.cap[:,:,5] .+ phS.p.*geoS.cap[:,:,5]
                @views ϕsave[snap,:,:] .= phL.ϕ.*geoL.cap[:,:,5] .+ phS.ϕ.*geoS.cap[:,:,5]
                @views Uxsave[snap,:,:] .= phL.u.*grid_u.geoL.cap[:,:,5] .+ phS.u.*geoS.capu[:,:,5]
                @views Uysave[snap,:,:] .= phL.v.*grid_v.geoL.cap[:,:,5] .+ phS.v.*geoS.capv[:,:,5]
                @views Uxcorrsave[snap,:,:] .= phL.ucorr.*grid_u.geoL.cap[:,:,5] .+ phS.ucorr.*geoS.capu[:,:,5]
                @views Uycorrsave[snap,:,:] .= phL.vcorr.*grid_v.geoL.cap[:,:,5] .+ phS.vcorr.*geoS.capv[:,:,5]
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd; step=snap)
            elseif ns_solid_phase
                @views psave[snap,:,:] .= phS.p
                @views ϕsave[snap,:,:] .= phS.ϕ
                @views Uxsave[snap,:,:] .= phS.u
                @views Uysave[snap,:,:] .= phS.v
                @views Uxcorrsave[snap,:,:] .= phS.ucorr
                @views Uycorrsave[snap,:,:] .= phS.vcorr
            elseif ns_liquid_phase
                @views psave[snap,:,:] .= phL.p
                @views ϕsave[snap,:,:] .= phL.ϕ
                @views Uxsave[snap,:,:] .= phL.u
                @views Uysave[snap,:,:] .= phL.v
                @views Uxcorrsave[snap,:,:] .= phL.ucorr
                @views Uycorrsave[snap,:,:] .= phL.vcorr
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd; step=snap)
            end
        end

        current_i += 1

        if adaptative_t
            τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))
        end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if heat || free_surface
                print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
                print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
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
            @show (length(MIXED))
        end
    end

    if levelset
        if save_radius || hill
            return MIXED, SOLID, LIQUID, radius
        end
        return MIXED, MIXED_u, MIXED_v, SOLID, LIQUID
    else
        return MIXED
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

        bcs!(faces, BC_u.left, dx[1,1])
        bcs!(faces, BC_u.right, dx[1,end])
        bcs!(faces, BC_u.bottom, dy[1,1])
        bcs!(faces, BC_u.top, dy[end,1])

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

            bcs!(faces, BC_u.left, dx[1,1])
            bcs!(faces, BC_u.right, dx[1,end])
            bcs!(faces, BC_u.bottom, dy[1,1])
            bcs!(faces, BC_u.top, dy[end,1])

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
