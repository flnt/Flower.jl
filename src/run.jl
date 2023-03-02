function run_forward(num, grid, grid_u, grid_v,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, 
    phS, phL, fwd, fwdS, fwdL;
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

    iRe = 1.0 / Re

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
        update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)

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

    # Heat
    set_heat!(num, grid, geoS, geoS.projection,
            opS, phS, HTS, bcTS, HuS, HvS,
            BC_TS, BC_uS, BC_vS,
            ind.MIXED, ind.LIQUID, heat_convection)
    
    set_heat!(num, grid, geoL, geoS.projection,
            opL, phL, HTL, bcTL, HuL, HvL,
            BC_TL, BC_uL, BC_vL,
            ind.MIXED, ind.SOLID, heat_convection)

    if ns_advection
        Cum1S .= opS.Cu * vec(phS.u) .+ opS.CUTCu
        Cvm1S .= opS.Cv * vec(phS.v) .+ opS.CUTCv
        Cum1L .= opL.Cu * vec(phL.u) .+ opL.CUTCu
        Cvm1L .= opL.Cv * vec(phL.v) .+ opL.CUTCv
    end

    CFL_sc = τ / Δ^2
    IIOE(LSA, LSB, u, V, ind.inside, CFL_sc, ny)

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

    tmp = copy(phS.T)
    tmp[2:end-1,2:end-1] .= θd
    init_borders!(tmp, BC_TS, θd)
    veci(phS.TD,grid,1) .= vec(phS.T)
    veci(phS.TD,grid,2) .= vec(tmp)
    @views fwdS.TD[1,:] .= phS.TD

    tmp = copy(phL.T)
    tmp[2:end-1,2:end-1] .= θd
    init_borders!(tmp, BC_TL, θd)    
    veci(phL.TD,grid,1) .= vec(phL.T)
    veci(phL.TD,grid,2) .= vec(tmp)
    @views fwdL.TD[1,:] .= phL.TD

    tmp = ones(grid_u.ny, grid_u.nx) .* num.u_inf
    tmp[2:end-1,2:end-1] .= 0.0
    veci(phS.uD,grid_u,1) .= vec(phS.u)
    veci(phS.uD,grid_u,2) .= vec(tmp)
    veci(phL.uD,grid_u,1) .= vec(phL.u)
    veci(phL.uD,grid_u,2) .= vec(tmp)
    @views fwdS.ucorrD[1,:,:] .= phS.uD
    @views fwdL.ucorrD[1,:,:] .= phL.uD

    tmp = ones(grid_v.ny, grid_v.nx) .* num.v_inf
    tmp[2:end-1,2:end-1] .= 0.0
    veci(phS.vD,grid_v,1) .= vec(phS.v)
    veci(phS.vD,grid_v,2) .= vec(tmp)
    veci(phL.vD,grid_v,1) .= vec(phL.v)
    veci(phL.vD,grid_v,2) .= vec(tmp)
    @views fwdS.vcorrD[1,:,:] .= phS.vD
    @views fwdL.vcorrD[1,:,:] .= phL.vD

    @views fwdL.pD[1,:] .= phL.pD
    @views fwdS.pD[1,:] .= phS.pD

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
                kill_dead_cells!(phS.T, grid, geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TS, geoS, BC_TS, ind.MIXED, geoS.projection,
                                        periodic_x, periodic_y)
                mul!(rhs, B, phS.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                @mytime (_, ch) = bicgstabl!(phS.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
                phS.T .= reshape(veci(phS.TD,grid,1), (ny, nx))
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                A_T, B, rhs = set_heat!(dir, num, grid, opC_TL, geoL, BC_TL, ind.MIXED, geoL.projection,
                                        periodic_x, periodic_y)
                mul!(rhs, B, phL.TD, 1.0, 1.0)
                @mytime blocks = DDM.decompose(A_T, grid.domdec, grid.domdec)

                @mytime (_, ch) = bicgstabl!(phL.TD, A_T, rhs, Pl=ras(blocks,grid.pou), log=true)
                println(ch)
                phL.T .= reshape(veci(phL.TD,grid,1), (ny, nx))
            end
        end

        if stefan
            update_stefan_velocity(num, grid, u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
        end

        if free_surface
            update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)
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
                level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
                try
                    utmp .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                catch
                    @error ("Inadequate level set function, iteration $current_i")
                    break
                end
                S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
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
                    if (heat && length(ind.MIXED) != 0)
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" mean(V[ind.MIXED]) findmax(V[ind.MIXED])[1] findmin(V[ind.MIXED])[1])
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[ind.MIXED]) findmax(κ[ind.MIXED])[1] findmin(κ[ind.MIXED])[1])
                    elseif free_surface
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
                catch
                    @show (ind.MIXED)
                end
            end
        end


        if levelset && (advection || current_i<2)
            update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)

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
                                                                                          ind.SOLID, ind.MIXED, periodic_x, periodic_y, current_i)
                else
                    Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S = projection_no_slip!(num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, phS,
                                                                                               BC_uS, BC_vS, BC_pS,
                                                                                               opC_pS, opC_uS, opC_vS,
                                                                                               Lum1_S, bc_Lum1_S, Lvm1_S, bc_Lvm1_S, Mum1_S, Mvm1_S,
                                                                                               ind.SOLID, ind.MIXED, periodic_x, periodic_y)
                end
                
            end
            if ns_liquid_phase
                if free_surface
                    Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mm1_L, Mum1_L, Mvm1_L = projection_fs!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                                                                                          BC_uL, BC_vL, BC_pL,
                                                                                          opC_pL, opC_uL, opC_vL,
                                                                                          Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                                                                                          ind.LIQUID, ind.MIXED, periodic_x, periodic_y, current_i)
                else
                   Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L = projection_no_slip!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                                                                                              BC_uL, BC_vL, BC_pL,
                                                                                              opC_pL, opC_uL, opC_vL,
                                                                                              Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                                                                                              ind.LIQUID, ind.MIXED, periodic_x, periodic_y)
                end
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
                @views fwdS.ucorrD[snap,:,:] .= phS.uD
                @views fwdL.ucorrD[snap,:,:] .= phL.uD
                @views fwdS.vcorrD[snap,:,:] .= phS.vD
                @views fwdL.vcorrD[snap,:,:] .= phL.vD
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, fwdL; step=snap)
            elseif ns_solid_phase
                @views fwdS.p[snap,:,:] .= phS.p
                @views fwdS.pD[snap,:] .= phS.pD
                @views fwdS.ϕ[snap,:,:] .= phS.ϕ
                @views fwdS.u[snap,:,:] .= phS.u
                @views fwdS.v[snap,:,:] .= phS.v
                @views fwdS.ucorr[snap,:,:] .= phS.ucorr
                @views fwdS.vcorr[snap,:,:] .= phS.vcorr
                @views fwdS.ucorrD[snap,:,:] .= phS.uD
                @views fwdS.vcorrD[snap,:,:] .= phS.vD
            elseif ns_liquid_phase
                @views fwdL.p[snap,:,:] .= phL.p
                @views fwdL.pD[snap,:] .= phL.pD
                @views fwdL.ϕ[snap,:,:] .= phL.ϕ
                @views fwdL.u[snap,:,:] .= phL.u
                @views fwdL.v[snap,:,:] .= phL.v
                @views fwdL.ucorr[snap,:,:] .= phL.ucorr
                @views fwdL.vcorr[snap,:,:] .= phL.vcorr
                @views fwdL.ucorrD[snap,:,:] .= phL.uD
                @views fwdL.vcorrD[snap,:,:] .= phL.vD
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, fwdL; step=snap)
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
