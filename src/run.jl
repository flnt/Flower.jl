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
    BC_trans_scal = Vector{Boundaries}(),
    BC_phi_eleL = BoundariesInt(),
    BC_int = [Wall()],
    time_scheme = CN,
    ls_scheme = weno5,
    auto_reinit = false,
    heat = false,
    heat_convection = false,
    heat_liquid_phase = false,
    heat_solid_phase = false,
    navier_stokes = false,
    ns_advection = false,
    ns_liquid_phase = false,
    ns_solid_phase = false,
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
    Ra = 0.0,
    λ = 1,
    electrolysis = false,
    electrolysis_convection = false,  
    electrolysis_liquid_phase = false,
    electrolysis_phase_change = false,
    electrolysis_reaction = "nothing",
    )
    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations,
            current_i, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso, nLS, _nLS,
            concentration0, diffusion_coeff, nb_transported_scalars, temperature0, i0, phi_ele0, phi_ele1, alphac,
            alphaa, Ru, Faraday = num
    @unpack opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL = op
    @unpack x, y, nx, ny, dx, dy, ind, LS, V = grid

    if length(BC_int) != nLS
        @error ("You have to specify $(nLS) boundary conditions.")
        return nothing
    end

    free_surface = false
    stefan = false
    if any(is_fs, BC_int)
        free_surface = true
    end
    if any(is_stefan, BC_int)
        stefan = true
    end

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan
        advection = true
    else
        advection = false
    end

    iRe = 1.0 / Re
    CFL_sc = τ / Δ^2

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
    utmp = copy(LS[1].u)
    rhs_LS = fzeros(grid)

    if levelset
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

        if save_radius
            n_snaps = iszero(max_iterations%save_every) ? max_iterations÷save_every+1 : max_iterations÷save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, LS[1])
        end
        if hill
            local radius = zeros(max_iterations+1)
            a = zeros(length(LS[1].MIXED))
            for i in eachindex(LS[1].MIXED)
                a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_u.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_v.LS[1].MIXED = [CartesianIndex(-1,-1)]
    end

    if save_length
        fwd.length[1] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
    end

    kill_dead_cells!(phS.T, grid, LS[end].geoS)
    kill_dead_cells!(phL.T, grid, LS[end].geoL)

    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    # TODO kill_dead_cells! for [:,:,iscal]
    if electrolysis
        for iscal=1:nb_transported_scalars
            # kill_dead_cells!(@view(phL.trans_scal[:,:,iscal]), grid, LS[end].geoL) #TODO  end or 1
            @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid, LS[end].geoL) #TODO  end or 1
        end
    end     
    ####################################################################################################

    if nb_transported_scalars>0
        printstyled(color=:green, @sprintf "\n after kill \n")
        print_electrolysis_statistics(phL)
        print("average", average!(phL.T, grid, LS[1].geoL, num))
    end

    for iLS in 1:_nLS
        @views fwd.u[iLS,1,:,:] .= LS[iLS].u
        @views fwd.ux[iLS,1,:,:] .= grid_u.LS[iLS].u
        @views fwd.uy[iLS,1,:,:] .= grid_v.LS[iLS].u
        @views fwd.κ[iLS,1,:,:] .= LS[iLS].κ
    end
    @views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]
    @views fwdL.T[1,:,:] .= phL.T
    @views fwdS.T[1,:,:] .= phS.T
    @views fwdS.p[1,:,:] .= phS.p
    @views fwdL.p[1,:,:] .= phL.p
    @views fwdS.u[1,:,:] .= phS.u
    @views fwdS.v[1,:,:] .= phS.v
    @views fwdL.u[1,:,:] .= phL.u
    @views fwdL.v[1,:,:] .= phL.v
    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    if electrolysis
        for iscal=1:nb_transported_scalars
            @views fwd.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal]
            @views fwdL.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal]
        end
        @views fwdL.phi_ele[1,:,:] .= phL.phi_ele
        @views fwdL.i_current_mag[1,:,:] .= phL.i_current_mag
    end     
    ####################################################################################################

    vec1(phS.TD,grid) .= vec(phS.T)
    vec2(phS.TD,grid) .= θd
    init_borders!(phS.TD, grid, BC_TS, θd)
    @views fwdS.TD[1,:] .= phS.TD

    vec1(phL.TD,grid) .= vec(phL.T)
    vec2(phL.TD,grid) .= θd
    init_borders!(phL.TD, grid, BC_TL, θd)
    @views fwdL.TD[1,:] .= phL.TD

    
    if electrolysis
        printstyled(color=:green, @sprintf "\n Check \n")
        print_electrolysis_statistics(phL)
    end

    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check %s %s %s %.2e %.2e\n" heat heat_convection electrolysis τ θd)

        for iscal=1:nb_transported_scalars
            @views vec1(phL.trans_scalD[:,iscal],grid) .= vec(phL.trans_scal[:,:,iscal])
            @views vec2(phL.trans_scalD[:,iscal],grid) .= concentration0[iscal]
            @views init_borders!(phL.trans_scalD[:,iscal], grid, BC_trans_scal[iscal], concentration0[iscal])
            @views fwdL.trans_scalD[1,:,iscal] .= phL.trans_scalD[:,iscal]
        end
        vec1(phL.phi_eleD,grid) .= vec(phL.phi_ele)
        vec2(phL.phi_eleD,grid) .= phi_ele0 #TODO
        init_borders!(phL.phi_eleD, grid, BC_phi_eleL, phi_ele0)
        @views fwdL.phi_eleD[1,:] .= phL.phi_eleD
    end     
    ####################################################################################################

    vec1(phS.uD,grid_u) .= vec(phS.u)
    vec2(phS.uD,grid_u) .= num.uD
    vecb(phS.uD,grid_u) .= num.u_inf
    vec1(phL.uD,grid_u) .= vec(phL.u)
    vec2(phL.uD,grid_u) .= num.uD
    vecb(phL.uD,grid_u) .= num.u_inf
    vec1(phS.ucorrD,grid_u) .= vec(phS.u)
    vec2(phS.ucorrD,grid_u) .= num.uD
    vecb(phS.ucorrD,grid_u) .= num.u_inf
    vec1(phL.ucorrD,grid_u) .= vec(phL.u)
    vec2(phL.ucorrD,grid_u) .= num.uD
    vecb(phL.ucorrD,grid_u) .= num.u_inf
    @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

    vec1(phS.vD,grid_v) .= vec(phS.v)
    vec2(phS.vD,grid_v) .= num.vD
    vecb(phS.vD,grid_v) .= num.v_inf
    vec1(phL.vD,grid_v) .= vec(phL.v)
    vec2(phL.vD,grid_v) .= num.vD
    vecb(phL.vD,grid_v) .= num.v_inf
    vec1(phS.vcorrD,grid_v) .= vec(phS.v)
    vec2(phS.vcorrD,grid_v) .= num.vD
    vecb(phS.vcorrD,grid_v) .= num.v_inf
    vec1(phL.vcorrD,grid_v) .= vec(phL.v)
    vec2(phL.vcorrD,grid_v) .= num.vD
    vecb(phL.vcorrD,grid_v) .= num.v_inf
    @views fwdS.vcorrD[1,:,:] .= phS.vcorrD
    @views fwdL.vcorrD[1,:,:] .= phL.vcorrD

    @views fwdL.pD[1,:] .= phL.pD
    @views fwdS.pD[1,:] .= phS.pD

    if is_FE(time_scheme) || is_CN(time_scheme)
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        if navier_stokes || heat || electrolysis
            geoS = [LS[iLS].geoS for iLS in 1:_nLS]
            geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
            geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS,
                opC_pS, opC_uS, opC_vS, periodic_x, periodic_y
            )

            geoL = [LS[iLS].geoL for iLS in 1:_nLS]
            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
                opC_pL, opC_uL, opC_vL, periodic_x, periodic_y
            )
        end

        Mm1_L = copy(opC_pL.M)
        Mm1_S = copy(opC_pS.M)
        Mum1_L = copy(opC_uL.M)
        Mum1_S = copy(opC_uS.M)
        Mvm1_L = copy(opC_vL.M)
        Mvm1_S = copy(opC_vS.M)

        if navier_stokes || heat || electrolysis
            AuS, BuS, _ = FE_set_momentum(
                BC_int, num, grid_u, opC_uS,
                false, false,
                iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                true
            )
            AvS, BvS, _ = FE_set_momentum(
                BC_int, num, grid_v, opC_vS,
                false, false,
                iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                true
            )
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            AϕS, _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pS, opC_uS, opC_vS,
                false, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                true
            )

            AuL, BuL, _ = FE_set_momentum(
                BC_int, num, grid_u, opC_uL,
                false, false,
                iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                true
            )
            AvL, BvL, _ = FE_set_momentum(
                BC_int, num, grid_v, opC_vL,
                false, false,
                iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                true
            )
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            AϕL, _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                false, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    V0S = volume(LS[end].geoS)
    V0L = volume(LS[end].geoL)

    current_t = 0.

    while current_i < max_iterations + 1

        if !stefan
            V .= speed*ones(ny, nx)
        end

        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                A_T, B, rhs = set_heat!(BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                                        opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                                        periodic_x, periodic_y, heat_convection)
                mul!(rhs, B, phS.TD, 1.0, 1.0)

                phS.TD .= A_T \ rhs
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                A_T, B, rhs = set_heat!(BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                                        opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                                        periodic_x, periodic_y, heat_convection)
                mul!(rhs, B, phL.TD, 1.0, 1.0)

                phL.TD .= A_T \ rhs
                phL.T .= reshape(veci(phL.TD,grid,1), grid)
            end
        end

        ####################################################################################################
        #Electrolysis
        ####################################################################################################  
        if electrolysis
            if electrolysis_liquid_phase

                for iscal=1:nb_transported_scalars
                    @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL)
                    @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

                    A_T, B, rhs = set_scalar_transport!(BC_trans_scal[iscal].int, num, grid, opC_TL, LS[1].geoL, phL, θd, BC_trans_scal[iscal],
                                                        LS[1].MIXED, LS[1].geoL.projection,
                                                        opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                                                        periodic_x, periodic_y, electrolysis_convection, diffusion_coeff[iscal])
                    mul!(rhs, B, phL.trans_scalD[:,iscal], 1.0, 1.0)

                    @views phL.trans_scalD[:,iscal] .= A_T \ rhs
                    @views phL.trans_scal[:,:,iscal] .= reshape(veci(phL.trans_scalD[:,iscal],grid,1), grid)
                end


                #TODO heat for electrolysis

                #TODO Poisson
                ####################################################################################################
                #Electrolysis: Poisson
                ####################################################################################################  

                #electroneutrality assumption


                a0_p = [] 
                for i in 1:num.nLS
                    push!(a0_p, zeros(grid))
                end

                # Constant electrical conductivity assumption
                #TODO electrical conductivity depends on concentration
                #iKOH index of KOH 
                # kappa_ele=2*Faraday^2*concentration0[iKOH]*diffusion_coeff[iKOH]/(Ru*T)
                # elec_cond=2*Faraday^2*trans_scal[iKOH]*diffusion_coeff[iKOH]/(Ru*T)

                #so after initialisation

                if heat
                    elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru.*phL.T) #phL.T
                else
                    elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru*temperature0) 
                end

                #TODO Poisson with variable coefficients
                #TODO need to iterate? since nonlinear

                #TODO no concentration prefactor

                #Update Butler-Volmer Boundary Condition with new potential 
               
                # eta = phi_ele1 .- phL.phi_ele[:,1]
                # i_current = i0*(exp(alphaa*Faraday*eta/(Ru*temperature0))-exp(-alphac*Faraday*eta/(Ru*temperature0)))

                if electrolysis_reaction == "Butler_no_concentration"

                    if heat
                        BC_phi_eleL.left.val = -butler_volmer_no_concentration.(alphaa,alphac,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,phL.T)./elec_cond[:,1]
                    else
                        BC_phi_eleL.left.val = -butler_volmer_no_concentration.(alphaa,alphac,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,temperature0)./elec_cond[:,1]
                    end    

                elseif electrolysis_reaction == "Butler_no_concentration"
                    # BC_phi_eleL.left.val = -butler_volmer_concentration.(alphaa,alphac,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,temperature0)./elec_cond
                end


                Aphi_eleL, rhs_phi_ele = set_poisson(
                    [BC_phi_eleL.int], num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                    false, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_phi_eleL,
                    true
                )

                # A, rhs = set_poisson(BC_int, num, grid, a0_p, op.opC_pL, op.opC_uL, op.opC_vL, 1, Lp, bc_Lp, bc_Lp_b, BC, true)
                # AϕL, _ = set_poisson(
                #     BC_int, num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                #     false, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                #     true
                # )

                # b = Δf.(
                #     grid.x .+ getproperty.(grid.LS[1].geoL.centroid, :x) .* grid.dx,
                #     grid.y .+ getproperty.(grid.LS[1].geoL.centroid, :y) .* grid.dy
                # )
                # b_phi_ele = zeros(
                #     grid.x .+ getproperty.(grid.LS[1].geoL.centroid, :x) .* grid.dx,
                #     grid.y .+ getproperty.(grid.LS[1].geoL.centroid, :y) .* grid.dy
                # )

                # b_phi_ele = zeros(
                #     grid.x .+ getproperty.(grid.LS[1].geoL.centroid, :x) .* grid.dx,
                #     grid.y .+ getproperty.(grid.LS[1].geoL.centroid, :y) .* grid.dy
                # )

                b_phi_ele = zeros(grid)


                veci(rhs_phi_ele,grid,1) .+= op.opC_pL.M * vec(b_phi_ele)
            
                res_phi_ele = zeros(size(rhs_phi_ele))
            
                @time @inbounds @threads for i in 1:Aphi_eleL.m
                    @inbounds Aphi_eleL[i,i] += 1e-10
                end
                
                # @time res_phi_ele .= Aphi_eleL \ rhs_phi_ele

                phL.phi_eleD .= Aphi_eleL \ rhs_phi_ele
                
                phL.phi_ele .= reshape(veci(phL.phi_eleD,grid,1), grid)

            
                # TODO compute magnitude of exchange current
                # gradient!(::Neumann, Ox, Oy, Bx, By, HNx, HNy, Divx, Divy, dcap, n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
                # TODO add post-treatment variables

                #TODO update BC concentration

                # compute_grad_phi_ele!(grid, phL, V, periodic_x, periodic_y) #TODO current
                compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, opC_pL, opC_pS) #TODO current

                # scal_magnitude


                if heat
                    elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru.*phL.T) #phL.T
                else
                    elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru*temperature0) 
                end

                phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude

                # IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)


                ####################################################################################################
                

            end
        end
        ####################################################################################################

        for iLS in 1:nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
            elseif is_fs(BC_int[iLS])
                if electrolysis_phase_change
                    update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, phL.ucorrD, phL.vcorrD, periodic_x, periodic_y, Vmean, phL.trans_scal[:,:,1],diffusion_coeff[1],concentration0[1],opC_pL)
                else
                    update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.ucorrD, phL.vcorrD, periodic_x, periodic_y)
                end
            end
            #TODO molar flux of hydrogen with Fick's law
            # $\frac{dm}{dt}=\int_S D \nabla c \cdot \mathbf{n} dS$
        end

        if verbose && adaptative_t
            println("τ = $τ")
        end

        if advection
            for (iLS, bc) in enumerate(BC_int)
                if is_stefan(bc)
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, (LS[iLS].B * vec(LS[iLS].u))), grid)
                    # u .= sqrt.((x .- current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
                elseif is_fs(bc)
                    # rhs_LS .= 0.0
                    # LS[iLS].A.nzval .= 0.0
                    # LS[iLS].B.nzval .= 0.0
                    # IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, θ_out, τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    # utmp .= reshape(gmres(LS[iLS].A, (LS[iLS].B * vec(LS[iLS].u))) .+ rhs_LS, grid)

                    # rhs_LS .= 0.0
                    # S2IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, utmp, LS[iLS].u, θ_out, τ, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, (LS[iLS].B * vec(LS[iLS].u))) .+ rhs_LS, grid)

                    # Project velocities to the normal and use advecion scheme for advection just
                    # in the normal direction
                    tmpVx = zeros(grid)
                    tmpVy = zeros(grid)
                    V .= 0.0
                    @inbounds @threads for II in grid.LS[iLS].MIXED
                        cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                        cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                        tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))

                        cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                        cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                        tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                        tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                        β = atan(tmpVy[II], tmpVx[II])
                        V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
                    end

                    i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
                    field_extension!(grid, grid.LS[iLS].u, grid.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                    rhs_LS .= 0.0
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, (LS[iLS].B * vec(LS[iLS].u))) .+ rhs_LS, grid)
                end
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*current_i*τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*current_i*τ);
            elseif nb_reinit > 0
                if auto_reinit && current_i%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)
                            println(ls_rg)
                            if ls_rg >= δreinit
                                println("yes")
                                RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                            end
                        end
                    end
                elseif current_i%num.reinit_every == 0
                    for iLS in 1:nLS
                        RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                    end
                end
            end
            # numerical breakup
            if free_surface && breakup
                count = breakup_f(LS[1].u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 1e-5)
                if count > 0
                    RK2_reinit!(ls_scheme, grid, ind, 1, LS[1].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                end
            end
        end

        if verbose
            if (current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (current_i-1) 100*(current_i-1)/max_iterations)
                if heat && length(LS[end].MIXED) != 0
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1])
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
                elseif advection && length(LS[end].MIXED) != 0
                    V_mean = mean([mean(grid_u.V[LS[1].MIXED]), mean(grid_v.V[LS[1].MIXED])])
                    V_max = max(findmax(grid_u.V[LS[1].MIXED])[1], findmax(grid_v.V[LS[1].MIXED])[1])
                    V_min = min(findmin(grid_u.V[LS[1].MIXED])[1], findmin(grid_v.V[LS[1].MIXED])[1])
                    print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1])
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
                        if electrolysis
                            print_electrolysis_statistics(phL)
                        end 
                    end
                end
            end
        end


        if levelset && (advection || current_i<2)
            NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)

            LS[end].geoL.fresh .= false
            LS[end].geoS.fresh .= false
            grid_u.LS[end].geoL.fresh .= false
            grid_u.LS[end].geoS.fresh .= false
            grid_v.LS[end].geoL.fresh .= false
            grid_v.LS[end].geoS.fresh .= false

            get_fresh_cells!(grid, grid.LS[end].geoS, Mm1_S, grid.ind.all_indices)
            get_fresh_cells!(grid, grid.LS[end].geoL, Mm1_L, grid.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoS, Mum1_S, grid_u.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoL, Mum1_L, grid_u.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoS, Mvm1_S, grid_v.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoL, Mvm1_L, grid_v.ind.all_indices)

            FRESH_L_u = findall(grid_u.LS[end].geoL.fresh)
            FRESH_S_u = findall(grid_u.LS[end].geoS.fresh)
            FRESH_L_v = findall(grid_v.LS[end].geoL.fresh)
            FRESH_S_v = findall(grid_v.LS[end].geoS.fresh)

            if navier_stokes
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,1), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,1), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phS.uD,grid_u,2), veci(phS.uD,grid_u,1),
                    grid_u.LS[end].geoS.projection, FRESH_S_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phS.vD,grid_v,2), veci(phS.vD,grid_v,1),
                    grid_v.LS[end].geoS.projection, FRESH_S_v, periodic_x, periodic_y)

                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
            end

            if iszero(current_i%save_every) || current_i==max_iterations
                snap = current_i÷save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, LS[1])
                end
                if hill
                    a = zeros(length(LS[1].MIXED))
                    for i in eachindex(LS[1].MIXED)
                        a[i] = LS[1].geoL.projection[LS[1].MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                if save_length
                    fwd.length[snap] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
                end
            end
        end

        if navier_stokes
            if !advection
                no_slip_condition!(num, grid, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
                # grid_u.V .= Δ / (1 * τ)
                # grid_v.V .= 0.0
            end

            if ns_solid_phase
                geoS = [LS[iLS].geoS for iLS in 1:_nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
                AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, AuS, BuS, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, AvS, BvS, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                    BC_uS, BC_vS, BC_pS,
                    opC_pS, opC_uS, opC_vS, opS,
                    AuS, BuS, AvS, BvS, AϕS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra
                )
            end
            if ns_liquid_phase
                geoL = [LS[iLS].geoL for iLS in 1:_nLS]
                geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
                geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, AuL, BvL, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, AvL, BvL, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    opC_pL, opC_uL, opC_vL, opL,
                    AuL, BuL, AvL, BvL, AϕL,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, advection, current_i, Ra
                )
                # if current_i == 1
                #     phL.u .= -0.5 .* grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy
                #     phL.v .= 0.5 .* grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
                #     phL.u[grid_u.LS[1].SOLID] .= 0.0
                #     phL.v[grid_v.LS[1].SOLID] .= 0.0
                # end
                # linear_advection!(
                #     num, grid, LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
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
            for iLS in 1:_nLS
                @views fwd.u[iLS,snap,:,:] .= LS[iLS].u
                @views fwd.ux[iLS,snap,:,:] .= grid_u.LS[iLS].u
                @views fwd.uy[iLS,snap,:,:] .= grid_v.LS[iLS].u
                @views fwd.κ[iLS,snap,:,:] .= LS[iLS].κ
            end

            if heat_solid_phase && heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T.*LS[end].geoS.cap[:,:,5]
            end
            if heat_solid_phase
                @views fwd.T[snap,:,:] .= phS.T
                @views fwdS.T[snap,:,:] .= phS.T
                @views fwdS.TD[snap,:] .= phS.TD
            end
            if heat_liquid_phase
                @views fwd.T[snap,:,:] .= phL.T
                @views fwdL.T[snap,:,:] .= phL.T
                @views fwdL.TD[snap,:] .= phL.TD
            end
            if electrolysis_liquid_phase #TODO

                for iscal=1:nb_transported_scalars
                    @views fwd.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal]
                    @views fwdL.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal]
                    # @views fwdL.trans_scalD[snap,:,iscal] .= phL.trans_scalD
                end
                @views fwdL.phi_ele[snap,:,:] .= phL.phi_ele

                @views fwdL.phi_eleD[snap,:] .= phL.phi_eleD

                @views fwdL.i_current_mag[snap,:,:] .= phL.i_current_mag


            end
            if ns_solid_phase
                @views fwdS.p[snap,:,:] .= phS.p
                @views fwdS.pD[snap,:] .= phS.pD
                @views fwdS.ϕ[snap,:,:] .= phS.ϕ
                @views fwdS.u[snap,:,:] .= phS.u
                @views fwdS.v[snap,:,:] .= phS.v
                @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
                @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
            end
            if ns_liquid_phase
                @views fwdL.p[snap,:,:] .= phL.p
                @views fwdL.pD[snap,:] .= phL.pD
                @views fwdL.ϕ[snap,:,:] .= phL.ϕ
                @views fwdL.u[snap,:,:] .= phL.u
                @views fwdL.v[snap,:,:] .= phL.v
                @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
                @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
                force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, phL; step=snap)
            end
            if advection
                fwdS.Vratio[snap] = volume(LS[end].geoS) / V0S
                fwdL.Vratio[snap] = volume(LS[end].geoL) / V0L
            end
        end
        # force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, phL; step=current_i+1)

        if electrolysis
        
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                any(isnan, phL.trans_scal) || any(isnan, phL.phi_ele) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8 || norm(phL.trans_scal) > 1e8 || norm(phL.phi_ele) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" current_i)
                return nothing
            end
        else
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" current_i)
                return nothing
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
            if stefan && advection
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
            end
            if free_surface && advection
                print(@sprintf "Vol_ratio = %.3f%%\n" (volume(LS[end].geoL) / V0L * 100))
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(V[LS[1].MIXED]) findmax(V[LS[1].MIXED])[1] findmin(V[LS[1].MIXED])[1] std(V[LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(LS[1].κ[LS[1].MIXED]) findmax(LS[1].κ[LS[1].MIXED])[1] findmin(LS[1].κ[LS[1].MIXED])[1] std(LS[1].κ[LS[1].MIXED]))
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
                    if electrolysis
                       print_electrolysis_statistics(phL)
                    end 
                end
            end
            print("\n\n")
        catch
            @show (length(LS[end].MIXED))
        end
    end

    if levelset && (save_radius || hill)
        return radius
    else
        return nothing
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
