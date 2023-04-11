
# function locate_index(gp, indices)
#     index_CL = zeros(Bool,size(indices)) # bolean vector all false with size of indices
#     @inbounds @threads for II in indices # loop over the indices
#         if gp.iso[II] in (1.,2.,6.,9.,13.,14.)
#             index_CL[II] = true  #then we are in a CL cell
#         end
#     end
#     return index_CL # CL contact line 
# end

function locate_index(gp, ind)
    index_CL = zeros(Bool, size(ind))
    @inbounds @threads for i in eachindex(ind)
        iso_value = gp.iso[ind[i]]
        index_CL[i] = iso_value != 0.0 && iso_value in (1., 2., 6., 9., 13., 14.) ? true : index_CL[i]
        #@show(i, ind[i], iso_value)
    end
    return index_CL
end


# # make this function general to work with any direction
function get_cut_points(gp, indices, index_CL)
    x_Cl_vec = zeros(size(indices)) # zero vector of size indices
    for i in findall(index_CL) # loop over true values of index_CL
    if indices == gp.ind.b_bottom[1]
        x = gp.cut_points[1,i] # n,n matrix
        if x[1].y == -0.5
            x_Cl_vec[i] = x[1].x
        elseif x[2].y == -0.5
            x_Cl_vec[i] = x[2].x
        end
    elseif indices == gp.ind.b_top[1]
        x = gp.cut_points[end,i]
        if x[1].y == 0.5
            x_Cl_vec[i] = x[1].x
        elseif x[2].y == 0.5
            x_Cl_vec[i] = x[2].x
        end
    elseif indices == gp.ind.b_left[1]
        x = gp.cut_points[i,1]
        if x[1].x == -0.5
            x_Cl_vec[i] = x[1].y
        elseif x[2].x == -0.5
            x_Cl_vec[i] = x[2].y
        end
    elseif indices == gp.ind.b_right[1]
        x = gp.cut_points[i,end]
        if x[1].x == 0.5
            x_Cl_vec[i] = x[1].y
        elseif x[2].x == 0.5
            x_Cl_vec[i] = x[2].y
        end
    end
end
    return x_Cl_vec
end

function compute_bell_function(gp,num,x_Cl_vec,index_CL,indices)
    bell_function = zero(x_Cl_vec)
    for i in findall(index_CL) # loop over all non-zeros
        if indices == gp.ind.b_bottom[1]
            rel_x = gp.x[1,i] .- x_Cl_vec[i]
        elseif indices == gp.ind.b_top[1]
            rel_x = gp.x[end,i] .- x_Cl_vec[i]
        elseif indices == gp.ind.b_left[1]
            rel_x = gp.y[i,1] .- x_Cl_vec[i]
        elseif indices == gp.ind.b_right[1]
            rel_x = gp.y[i,end] .- x_Cl_vec[i]
        end
        bell_function .+= (1.0-tanh^2(rel_x/num.εCA))/num.εCA
    end 
    return bell_function
end

# value to fill a0 
function compute_young_stress(gp,num,indices)
    if indices == gp.ind.b_bottom[1]
        index_CL = locate_index(gp, indices)
        x_Cl_vec = get_cut_points(gp, indices, index_CL)
        bellf = compute_bell_function(gp,num,x_Cl_vec,index_CL,indices)
        return bellf[:].*(1.0/num.Ca).*(cos.(gp.α[1,:]).-cos(num.θe*π/180))
    end
end

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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re,
            max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid


    local NB_indices;

    local Cum1L = zeros(grid_u.nx*grid_u.ny)
    local Cvm1L = zeros(grid_v.nx*grid_v.ny)

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

    if levelset
        
        ind.MIXED = [CartesianIndex(-1,-1)]
        grid_u.ind.MIXED = [CartesianIndex(-1,-1)]
        grid_v.ind.MIXED = [CartesianIndex(-1,-1)]


    if ns_advection

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
    veci(phL.uD,grid_u,1) .= vec(phL.u)
    veci(phL.uD,grid_u,2) .= vec(tmp)

    veci(phS.ucorrD,grid_u,1) .= vec(phS.u)
    veci(phS.ucorrD,grid_u,2) .= vec(tmp)
    veci(phL.ucorrD,grid_u,1) .= vec(phL.u)
    veci(phL.ucorrD,grid_u,2) .= vec(tmp)
    @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

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

        if free_surface
            update_free_surface_velocity(num, grid_u, grid_v, phL.uD, phL.vD, periodic_x, periodic_y)
        end

        if advection
            CFL_sc = τ / Δ^2
                        # using Inflow-Implicit Outflow-Explicit (IIOE) Eq. 12 from Mikula (2014)
            # First θ=1/2 (constant coefficient; Stefan problem)
            # Second θ=min() following Eq. 19/20 from Mikula (2014) (general case)
            #level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            level_update_IIOE!(grid, grid_u, grid_v, LSA, LSB, θ_out, ind.MIXED, τ, false, false)

            try
                #utmp .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                tmp_tracer .= reshape(gmres(LSA,(LSB*vec(tracer))), (ny,nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end
            #S2IIOE!(grid, grid_u, grid_v, LSA, LSB, utmp, u, θ_out, ind.MIXED, τ, periodic_x, periodic_y)
            S2IIOE!(grid, grid_u, grid_v, LSA, LSB, tmp_tracer, tracer, θ_out, MIXED, τ, false, false)
            try
                #u .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))
                tracer .= reshape(gmres(LSA,(LSB*vec(tracer))), (ny,nx))
            catch
                @error ("Inadequate level set function, iteration $current_i")
                break
            end

            if nb_reinit > 0
                if current_i%num.reinit_every == 0
                    #FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
                    FE_reinit(grid, ind, tracer, nb_reinit, Boundaries(), false, false)
                end
            end
            # # numerical breakup
            # if free_surface
            #     count = breakup(u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
            #     if count > 0
            #         FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)
            #     end
            # end
        end

        if verbose
            if current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
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
                        normpL = norm(phL.p.*τ)
                        print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    end
                catch
                    @show (ind.MIXED)
                end
            end
        end


        if levelset && (advection || current_i<2)
            NB_indices = update_ls_data(num, grid, grid_u, grid_v, u, κ, periodic_x, periodic_y)

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, ind.MIXED)
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
                    Mm1_S, Mum1_S, Mvm1_S = projection_fs!(num, grid, geoS, grid_u, grid_u.geoS, grid_v, grid_v.geoS, phS,
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
                    Mm1_L, Mum1_L, Mvm1_L = projection_fs!(num, grid, geoL, grid_u, grid_u.geoL, grid_v, grid_v.geoL, phL,
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
            @views Tsave[snap,:,:] .= tracer
            @views fwdL.p[snap,:,:] .= phL.p
            @views fwdL.ϕ[snap,:,:] .= phL.ϕ
            @views fwdL.u[snap,:,:] .= phL.u
            @views fwdL.v[snap,:,:] .= phL.v
            @views fwdL.ucorr[snap,:,:] .= phL.ucorr
            @views fwdL.vcorr[snap,:,:] .= phL.vcorr
        end
        current_i += 1
        if adaptative_t
            τ = min(CFL*Δ^2*Re, CFL*Δ/max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...))
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
            @show (length(ind.MIXED))
        end
    end