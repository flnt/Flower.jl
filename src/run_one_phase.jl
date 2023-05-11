function run_forward_one_phase(num, grid, grid_u, grid_v, opL, opC_pL, opC_uL, opC_vL, phL, fwd, fwdL, tracer;
    periodic_x=false,
    periodic_y=false,
    BC_pL=Boundaries(
        left=Boundary(),
        right=Boundary(),
        bottom=Boundary(),
        top=Boundary()),
    BC_uL=Boundaries(
        left=Boundary(),
        right=Boundary(),
        bottom=Boundary(),
        top=Boundary()),
    BC_vL=Boundaries(
        left=Boundary(),
        right=Boundary(),
        bottom=Boundary(),
        top=Boundary()),
    BC_u=Boundaries(
        left=Boundary(),
        right=Boundary(),
        bottom=Boundary(),
        top=Boundary()),
    advection=false, #move the level set
    ns_advection=false,
    navier_stokes=false,
    levelset=true,
    verbose=false,
    adaptative_t=false,
    show_every=100
)

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re,
    max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, θext, LSA, LSB, α = grid

    local NB_indices, Mm1_L, Mum1_L, Mvm1_L
    local Cum1L = zeros(grid_u.nx * grid_u.ny)
    local Cvm1L = zeros(grid_v.nx * grid_v.ny)

    θ_out = zeros(grid, 4)
    tmp_tracer = copy(tracer)
    LSRHS = fzeros(grid)

    if periodic_x
        BC_u.left.ind = ind.b_left
        BC_u.right.ind = ind.b_right
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = ind.b_left
        BC_u.right.ind = ind.b_right
    end

    if periodic_y
        BC_u.bottom.ind = ind.b_bottom
        BC_u.top.ind = ind.b_top
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = ind.b_bottom
        BC_u.top.ind = ind.b_top
    end

    NB_indices = update_ls_data(num, grid, grid_u, grid_v, tracer, κ, periodic_x, periodic_y)
    gp_iso_tmp = copy(grid.iso)
    gp_α_tmp = copy(grid.α)
    gpu_iso_tmp = copy(grid_u.iso)
    gpu_α_tmp = copy(grid_u.α)
    gpv_iso_tmp = copy(grid_v.iso)
    gpv_α_tmp = copy(grid_v.α)

    NB_indices = update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)
    grid.iso .= gp_iso_tmp
    grid.α .= gp_α_tmp
    grid_u.iso .= gpu_iso_tmp
    grid_u.α .= gpu_α_tmp
    grid_v.iso .= gpv_iso_tmp
    grid_v.α .= gpv_α_tmp

    if ns_advection
        Cum1L .= opL.Cu * vec(phL.u) .+ opL.CUTCu
        Cvm1L .= opL.Cv * vec(phL.v) .+ opL.CUTCv
    end

    CFL_sc = τ / Δ^2
    IIOE(LSA, LSB, u, V, ind.inside, CFL_sc, ny)

    @views fwd.u[1, :, :] .= u
    @views fwd.ux[1, :, :] .= grid_u.u
    @views fwd.uy[1, :, :] .= grid_v.u
    @views fwdL.p[1, :, :] .= phL.p
    @views fwdL.u[1, :, :] .= phL.u
    @views fwdL.v[1, :, :] .= phL.v

    tmp = ones(grid_u.ny, grid_u.nx) .* num.u_inf
    tmp[2:end-1, 2:end-1] .= 0.0
    veci(phL.uD, grid_u, 1) .= vec(phL.u)
    veci(phL.uD, grid_u, 2) .= vec(tmp)

    veci(phL.ucorrD, grid_u, 1) .= vec(phL.u)
    veci(phL.ucorrD, grid_u, 2) .= vec(tmp)
    @views fwdL.ucorrD[1, :, :] .= phL.ucorrD

    tmpv = ones(grid_v.ny, grid_v.nx) .* num.v_inf
    tmpv[2:end-1, 2:end-1] .= 0.0
    veci(phL.vD, grid_v, 1) .= vec(phL.v)
    veci(phL.vD, grid_v, 2) .= vec(tmpv)

    _, _, Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L = set_laplacians!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
        opC_pL, opC_uL, opC_vL, periodic_x, periodic_y)
    Mm1_L, Mum1_L, Mvm1_L = copy(opC_pL.M), copy(opC_uL.M), copy(opC_vL.M)

    #-----------------------------------------------------------------------------------------------------
    # LOOP STARTS HERE
    current_t = 0.0
    update_levelset_matrices(periodic_x, periodic_y, grid, LSA, LSB)

    # for i in eachindex(grid.α)
    #     if !isnan(grid.α[i])
    #         θext[i] = π * 0.50
    #     end
    # end

    while current_i < max_iterations + 1

        #update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)

        # 1 is bulk field 
        grid_u.V .= reshape(veci(phL.uD, grid_u, 1), (grid_u.ny, grid_u.nx))
        grid_v.V .= reshape(veci(phL.vD, grid_v, 1), (grid_v.ny, grid_v.nx))

        update_levelset_contact_angle(periodic_x, periodic_y, grid, tracer, LSRHS) # here only gp 

        if advection
            CFL_sc = τ / Δ^2
            # using Inflow-Implicit Outflow-Explicit (IIOE) Eq. 12 from Mikula (2014)
            # First θ=1/2 (constant coefficient; Stefan problem)
            # Second θ=min() following Eq. 19/20 from Mikula (2014) (general case)
            level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, false, false)

            try
                tmp_tracer .= reshape(gmres(LSA, (LSB * vec(tracer) .+ LSRHS); verbose=false, log=false), (ny, nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            S2IIOE!(grid, grid_u, grid_v, LSA, LSB, tmp_tracer, tracer, θ_out, ind.MIXED, τ, false, false)
            try
                tracer .= reshape(gmres(LSA, (LSB * vec(tracer) .+ LSRHS); verbose=false, log=false), (ny, nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end

            # if nb_reinit > 0
            #     if current_i%num.reinit_every == 0
            #         FE_reinit(grid, ind, tracer, nb_reinit, Boundaries(), false, false)
            #     end
            # end
            # if free_surface # numerical breakup
            #     count = breakup(u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
            #     if count > 0
            #         FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
            #     end
            # end
        end

        if verbose
            if current_i % show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i - 1) 100 * (current_i - 1) / max_iterations)
                    print(@sprintf "t = %3.2f  dt = %.6f\n" current_t τ)
                    if length(ind.MIXED) != 0
                        V_mean = mean([mean(grid_u.V[ind.MIXED]), mean(grid_v.V[ind.MIXED])])
                        V_max = max(findmax(grid_u.V[ind.MIXED])[1], findmax(grid_v.V[ind.MIXED])[1])
                        V_min = min(findmin(grid_u.V[ind.MIXED])[1], findmin(grid_v.V[ind.MIXED])[1])
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" V_mean V_max V_min)
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[ind.MIXED]) findmax(κ[ind.MIXED])[1] findmin(κ[ind.MIXED])[1])
                    end
                    if navier_stokes
                        normuL = norm(phL.u)
                        normvL = norm(phL.v)
                        normpL = norm(phL.p .* τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                catch
                    @show (ind.MIXED)
                end
            end
        end

        if levelset && (advection || current_i < 2)

            NB_indices = update_ls_data(num, grid, grid_u, grid_v, tracer, κ, periodic_x, periodic_y)
            gp_iso_tmp = copy(grid.iso)
            gp_α_tmp = copy(grid.α)
            gpu_iso_tmp = copy(grid_u.iso)
            gpu_α_tmp = copy(grid_u.α)
            gpv_iso_tmp = copy(grid_v.iso)
            gpv_α_tmp = copy(grid_v.α)

            NB_indices = update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)
            grid.iso .= gp_iso_tmp
            grid.α .= gp_α_tmp
            grid_u.iso .= gpu_iso_tmp
            grid_u.α .= gpu_α_tmp
            grid_v.iso .= gpv_iso_tmp
            grid_v.α .= gpv_α_tmp

            if iszero(current_i % save_every) || current_i == max_iterations
                snap = current_i ÷ save_every + 1
            end
        end

        #if contact_lines
        update_Young_stress(ind.MIXED, grid_u, grid_v, num)

        if navier_stokes
            no_slip_condition!(grid, grid_u, grid_v)
            Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L = projection_no_slip!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                BC_uL, BC_vL, BC_pL,
                opC_pL, opC_uL, opC_vL,
                Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                ind.LIQUID, ind.MIXED, periodic_x, periodic_y)
        end

        current_t += τ
        if iszero(current_i % save_every) || current_i == max_iterations
            snap = current_i ÷ save_every + 1
            if current_i == max_iterations
                snap = size(fwd.T, 1)
            end
            fwd.t[snap] = current_t
            @views fwd.V[snap, :, :] .= V
            @views fwd.u[snap, :, :] .= u
            @views fwd.ux[snap, :, :] .= grid_u.u
            @views fwd.uy[snap, :, :] .= grid_v.u
            @views fwdL.T[snap, :, :] .= tracer
            @views fwdL.p[snap, :, :] .= phL.p
            @views fwdL.ϕ[snap, :, :] .= phL.ϕ
            @views fwdL.u[snap, :, :] .= phL.u
            @views fwdL.v[snap, :, :] .= phL.v
            @views fwdL.ucorr[snap, :, :] .= phL.ucorr
            @views fwdL.vcorr[snap, :, :] .= phL.vcorr
        end
        current_i += 1
        if adaptative_t
            τ = min(CFL * Δ^2 * Re, CFL * Δ / max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)...))
        end
    end
    # LOOP ENDS HERE
    #-----------------------------------------------------------------------------------------------------
    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i - 1) 100 * (current_i - 1) / max_iterations)
            if navier_stokes
                normuL = norm(phL.u)
                normvL = norm(phL.v)
                normpL = norm(phL.p .* τ)
                print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
            end
            print("\n\n")
        catch
            @show (length(ind.MIXED))
        end
    end
end