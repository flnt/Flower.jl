function run_forward_one_phase(num, grid, grid_u, grid_v,
    opL, opC_pL, opC_uL, opC_vL, 
    phL, fwd, tracer;
    periodic_x = false,
    periodic_y = false,
    BC_pL = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    BC_uL = Boundaries(
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
    advection = false, #move the level set
    ns_advection = false,
    navier_stokes = false,
    levelset = true,
    verbose = false,
    adaptative_t = false,
    show_every = 100,
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid
    @unpack Tall, usave, uusave, uvsave, TSsave, TLsave, Tsave, psave, ϕsave, Uxsave, Uysave, Uxcorrsave, Uycorrsave, Vsave, κsave, lengthsave, tv, Cd, Cl = fwd

    local MIXED; local SOLID; local LIQUID;
    local MIXED_vel_ext; local SOLID_vel_ext; local LIQUID_vel_ext;
    local MIXED_u_vel_ext; local SOLID_u_vel_ext; local LIQUID_u_vel_ext;
    local MIXED_v_vel_ext; local SOLID_v_vel_ext; local LIQUID_v_vel_ext;
    local MIXED_u; local SOLID_u; local LIQUID_u;
    local MIXED_v; local SOLID_v; local LIQUID_v;

    local Cum1L = zeros(grid_u.nx*grid_u.ny)
    local Cvm1L = zeros(grid_v.nx*grid_v.ny)

    local Lum1_L
    local bc_Lum1_L
    local Lvm1_L
    local bc_Lvm1_L
    local Mm1_L
    local Mum1_L
    local Mvm1_L

    θ_out = zeros(grid, 4)
    utmp = copy(u)

    tmp_tracer = copy(tracer)

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

    usave[1,:,:] .= u
    psave[1,:,:] .= phL.p 
    Uxsave[1,:,:] .= phL.u 
    Uysave[1,:,:] .= phL.v 

    if levelset

        grid.mid_point .= [Point(0.0, 0.0)]
        grid_u.mid_point .= [Point(0.0, 0.0)]
        grid_v.mid_point .= [Point(0.0, 0.0)]
        
        marching_squares!(num, grid)
        interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
        uusave[1,:,:] .= grid_u.u
        uvsave[1,:,:] .= grid_v.u
        marching_squares!(num, grid_u)        
        marching_squares!(num, grid_v)

        MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.all_indices, nx, ny, periodic_x, periodic_y)
        MIXED_u_vel_ext, SOLID_u_vel_ext, LIQUID_u_vel_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
        MIXED_v_vel_ext, SOLID_v_vel_ext, LIQUID_v_vel_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
        MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
        MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
        MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

        get_iterface_location!(grid, MIXED)
        get_iterface_location!(grid_u, MIXED_u)
        get_iterface_location!(grid_v, MIXED_v)
        get_interface_location_borders!(grid_u, periodic_x, periodic_y)
        get_interface_location_borders!(grid_v, periodic_x, periodic_y)

        get_curvature(num, grid, MIXED, periodic_x, periodic_y)
        postprocess_grids!(grid, grid_u, grid_v, MIXED, MIXED_u, MIXED_v, periodic_x, periodic_y, advection, ϵ)
        _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
        _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
        _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
        indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
        field_extension!(grid, grid.κ, indices_vel_ext, NB, periodic_x, periodic_y)

    elseif !levelset
        MIXED = [CartesianIndex(-1,-1)]
        MIXED_u = [CartesianIndex(-1,-1)]
        MIXED_v = [CartesianIndex(-1,-1)]
    end

    if ns_advection
        Cum1L .= opL.Cu * vec(phL.u) .+ opL.CUTCu
        Cvm1L .= opL.Cv * vec(phL.v) .+ opL.CUTCv
    end

    tmpu = ones(grid_u.ny, grid_u.nx) .* num.u_inf
    tmpu[2:end-1,2:end-1] .= 0.0
    veci(phL.uD,grid_u,1) .= vec(phL.u)
    veci(phL.uD,grid_u,2) .= vec(tmpu)

    tmpv = ones(grid_v.ny, grid_v.nx) .* num.v_inf
    tmpv[2:end-1,2:end-1] .= 0.0
    veci(phL.vD,grid_v,1) .= vec(phL.v)
    veci(phL.vD,grid_v,2) .= vec(tmpv)

    _, _, Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L = set_laplacians!(grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL,
                                                opC_pL, opC_uL, opC_vL, periodic_x, periodic_y)
    
    Mm1_L = copy(opC_pL.M)
    Mum1_L = copy(opC_uL.M)
    Mvm1_L = copy(opC_vL.M)

#-----------------------------------------------------------------------------------------------------
# LOOP STARTS HERE
    current_t = 0.

    # Imposing Newman boundary condition to the levelset function
    # LSA contains the product of the mass matrix and the Laplacian for the advection scheme
    # LSAΦ^{n+1}=LSBΦ^{n}+LSC*u_target 
    # 1*Φ^{n+1}=1*Φ^{n}
    # Φ^{n+1}_{i,j}+Φ^{n+1}_{i,j+1}=0

    if ! periodic_y
        for (II,JJ) in zip(grid.ind.b_bottom[1], grid.ind.b_bottom[2]) # first and second rows
        i = lexicographic(II, grid.ny)
        j = lexicographic(JJ, grid.ny)
        LSA[i,i] = 1. # set 1st row of LSA to 1
        LSB[i,i] = 0. # set 1st row of LSB to 0
        LSA[i,j] = -1. # set 2nd row of LSA to -1
        end
    end
#-----------------------------------------------------------------------------------------------------
    # LOOP STARTS HERE
    while current_i < max_iterations + 1
 
        grid_u.V .= reshape(veci(phL.uD,grid_u,1), (grid_u.ny, grid_u.nx))
        grid_v.V .= reshape(veci(phL.vD,grid_v,1), (grid_v.ny, grid_v.nx))

        if advection
            # using Inflow-Implicit Outflow-Explicit (IIOE) Eq. 12 from Mikula (2014)
            # First θ=1/2 (constant coefficient; Stefan problem)
            # Second θ=min() following Eq. 19/20 from Mikula (2014) (general case)
            level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, MIXED, τ, false, false)
            try
                tmp_tracer .= reshape(gmres(LSA,(LSB*vec(tracer))), (ny,nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            S2IIOE!(grid, grid_u, grid_v, LSA, LSB, tmp_tracer, tracer, θ_out, MIXED, τ, false, false)
            try
                tracer .= reshape(IterativeSolvers.gmres(LSA,(LSB*vec(tracer))), (ny,nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            if nb_reinit > 0
                if current_i%num.reinit_every == 0
                    FE_reinit(grid, ind, tracer, nb_reinit, Boundaries(), false, false)
                end
            end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                    print(@sprintf "t = %3.2f  dt = %.6f\n" current_t τ)
                    if length(MIXED) != 0
                        V_mean = mean([mean(grid_u.V[MIXED]), mean(grid_v.V[MIXED])])
                        V_max = max(findmax(grid_u.V[MIXED])[1], findmax(grid_v.V[MIXED])[1])
                        V_min = min(findmin(grid_u.V[MIXED])[1], findmin(grid_v.V[MIXED])[1])
                        print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f\n" V_mean V_max V_min)
                        print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1])
                    end
                    if navier_stokes
                        normuL = norm(phL.u)
                        normvL = norm(phL.v)
                        normpL = norm(phL.p.*τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                catch
                    @show (MIXED)
                end
            end
        end

        if levelset && current_i<2
            update_ls_data(num, grid, grid_u, grid_v, u, periodic_x, periodic_y)
            grid.α .= NaN
            grid_u.α .= NaN
            grid_v.α .= NaN
            faces .= 0.
            grid_u.faces .= 0.
            grid_v.faces .= 0.
            grid.mid_point .= [Point(0.0, 0.0)]
            grid_u.mid_point .= [Point(0.0, 0.0)]
            grid_v.mid_point .= [Point(0.0, 0.0)]

            marching_squares!(num, grid)
            interpolate_scalar!(grid, grid_u, grid_v, u, grid_u.u, grid_v.u)
            marching_squares!(num, grid_u)
            marching_squares!(num, grid_v)

            MIXED_vel_ext, SOLID_vel_ext, LIQUID_vel_ext = get_cells_indices(iso, ind.all_indices, nx, ny, periodic_x, periodic_y)
            MIXED_u_vel_ext, SOLID_u_vel_ext, LIQUID_u_vel_ext = get_cells_indices(grid_u.iso, grid_u.ind.all_indices, grid_u.nx, grid_u.ny, periodic_x, periodic_y)
            MIXED_v_vel_ext, SOLID_v_vel_ext, LIQUID_v_vel_ext = get_cells_indices(grid_v.iso, grid_v.ind.all_indices, grid_v.nx, grid_v.ny, periodic_x, periodic_y)
            MIXED, SOLID, LIQUID = get_cells_indices(iso, ind.all_indices)
            MIXED_u, SOLID_u, LIQUID_u = get_cells_indices(grid_u.iso, grid_u.ind.all_indices)
            MIXED_v, SOLID_v, LIQUID_v = get_cells_indices(grid_v.iso, grid_v.ind.all_indices)

            get_iterface_location!(grid, MIXED)
            get_iterface_location!(grid_u, MIXED_u)
            get_iterface_location!(grid_v, MIXED_v)
            get_interface_location_borders!(grid_u, periodic_x, periodic_y)
            get_interface_location_borders!(grid_v, periodic_x, periodic_y)

            geoL.emptied .= false
            geoS.emptied .= false
            grid_u.geoL.emptied .= false
            grid_u.geoS.emptied .= false
            grid_v.geoL.emptied .= false
            grid_v.geoS.emptied .= false

            get_curvature(num, grid, MIXED, periodic_x, periodic_y)
            postprocess_grids!(grid, grid_u, grid_v, MIXED, MIXED_u, MIXED_v, periodic_x, periodic_y, advection, ϵ)

            _MIXED_L_vel_ext = intersect(findall(geoL.emptied), MIXED_vel_ext)
            _MIXED_S_vel_ext = intersect(findall(geoS.emptied), MIXED_vel_ext)
            _MIXED_vel_ext = vcat(_MIXED_L_vel_ext, _MIXED_S_vel_ext)
            indices_vel_ext = vcat(SOLID_vel_ext, _MIXED_vel_ext, LIQUID_vel_ext)
            field_extension!(grid, grid.κ, indices_vel_ext, NB, periodic_x, periodic_y)
        end

        if navier_stokes
            no_slip_condition!(grid, grid_u, grid_v)
            Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L = projection_no_slip!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
                                                                                              BC_uL, BC_vL, BC_pL,
                                                                                              opC_pL, opC_uL, opC_vL,
                                                                                              Lum1_L, bc_Lum1_L, Lvm1_L, bc_Lvm1_L, Mum1_L, Mvm1_L,
                                                                                              LIQUID, MIXED, periodic_x, periodic_y)
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
            @views Tsave[snap,:,:] .= tracer

            @views psave[snap,:,:] .= phL.p[:,:]
            @views ϕsave[snap,:,:] .= phL.ϕ[:,:]
            @views Uxsave[snap,:,:] .= phL.u[:,:]
            @views Uysave[snap,:,:] .= phL.v[:,:]
            @views Uxcorrsave[snap,:,:] .= phL.ucorr[:,:]
            @views Uycorrsave[snap,:,:] .= phL.vcorr[:,:]
        end
        current_i += 1
        if adaptative_t
            τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)...))
        end
    end
    # LOOP ENDS HERE
#-----------------------------------------------------------------------------------------------------

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
            if navier_stokes
                normuL = norm(phL.u)
                normvL = norm(phL.v)
                normpL = norm(phL.p.*τ)
                print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
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