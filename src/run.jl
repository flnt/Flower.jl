function run_forward(
    num, grid, grid_u, grid_v, op, phS, phL;
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
    BC_trans_scal = Vector{BoundariesInt}(),
    BC_phi_ele = BoundariesInt(),
    BC_int = [WallNoSlip()],
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
    electrolysis_solid_phase = false,
    electrolysis_phase_change_case = "Khalighi",
    electrolysis_reaction = "nothing",
    bulk_conductivity = true,
    imposed_velocity = "none",
    adapt_timestep_mode = 0,
    non_dimensionalize=1,
    mode_2d=0,
    conductivity_mode = 1,
    convection_Cdivu = 0,
    test_laplacian = false,
    )

    # print("\n before unpack \n")

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso, nLS, _nLS, nNavier,
            concentration0, diffusion_coeff, nb_transported_scalars, nb_saved_scalars, temperature0, i0, phi_ele0, phi_ele1, alpha_c,
            alpha_a, Ru, Faraday, mu1, rho1 = num
    @unpack opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL = op
    @unpack x, y, nx, ny, dx, dy, ind, LS, V = grid

    # print("\n after unpack \n")

    if num.epsilon_mode == 1 || num.epsilon_mode ==2
        num.epsilon_dist = eps(0.01) * num.Δ
        num.epsilon_vol = (eps(0.01)*num.Δ)^2
        #TODO kill dead cells
        #TODO 1e-...
    end


    if length(BC_int) != nLS
        @error ("You have to specify $(nLS) boundary conditions.")
        return nothing
    end
    crashed=false
    # num.current_i=0
    free_surface = false
    stefan = false
        navier = false
    if any(is_fs, BC_int)
        free_surface = true
    end
    if any(is_stefan, BC_int)
        stefan = true
    end
    if any(is_navier_cl, BC_int) || any(is_navier, BC_int)
        navier = true
    end

    if nNavier > 1
        @warn ("When using more than 1 Navier BC, the interfaces shouldn't cross")
    end

    if electrolysis
        electric_potential = true 
        electrolysis_advection = true
        if nb_saved_scalars<1
            @error("number of saved scalars should be superior to one")
            return
        end

    end

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan || electrolysis_phase_change_case !="none"
        advection = true
    else
        advection = false
    end

    # The count threshold shouldn't be smaller than 2
    count_limit_breakup = 6
    
    printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL num.τ)
    if adapt_timestep_mode !=0
        num.τ = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
        # print("after adapt_timestep!")
        printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e num.τ : %.2e\n" CFL num.τ num.τ)
    end

    pres_free_surfaceS = 0.0
    # pres_free_surfaceL = 0.0
    pres_free_surfaceL = num.pres0

    if electrolysis_phase_change_case == "levelset"
        jump_mass_fluxS = false 
        jump_mass_fluxL = true
    else
        jump_mass_fluxS = false 
        jump_mass_fluxL = false
    end

    mass_fluxS = 0.0
    mass_fluxL = 0.0
    

    iRe = 1.0 / Re
    CFL_sc = num.τ / Δ^2

    irho = 1.0

    if non_dimensionalize==0
        #force L=1 u=1
        Re=rho1/mu1
        iRe = 1.0/Re
        irho=1.0/rho1
        num.visc_coeff=iRe
    else 
        num.visc_coeff=iRe
    end

    printstyled(color=:green, @sprintf "\n Re : %.2e %.2e\n" Re num.visc_coeff)


    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    current_radius = 0.0
    # TODO kill_dead_cells! for [:,:,iscal]
    if electrolysis
        current_radius = num.R

        printstyled(color=:green, @sprintf "\n radius: %.2e \n" current_radius)


        p_liq= num.pres0 + mean(veci(phL.pD,grid,2)) #TODO here one bubble
        p_g=p_liq + 2 * num.σ / current_radius

        #TODO using temperature0
        if mode_2d==0
            nH2 = p_g * 4.0 / 3.0 * pi * current_radius ^ 3 / (temperature0 * num.Ru) 
        elseif mode_2d == 1 #reference thickness for a cylinder
            nH2 = p_g * pi * current_radius ^ 2 * num.ref_thickness_2d / (temperature0 * num.Ru) 
        elseif mode_2d==2 #mol/meter
            nH2=concentration0[1]* pi * current_radius ^ 2
        elseif mode_2d==3 #mol/meter half circle
            nH2=1.0/2.0*concentration0[1]* pi * current_radius ^ 2
        end
        # nH2 = 4.0/3.0 * pi * current_radius^3 * num.rho2 / num.MWH2

        printstyled(color=:green, @sprintf "\n Mole: %.2e \n" nH2)

        printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(temperature0*num.Ru))



    end     
    ####################################################################################################

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

        # At every iteration, update_all_ls_data is called twice, 
        # once inside run.jl and another one (if there's advection of the levelset) inside set_heat!. 
        # The difference between both is a flag as last argument, inside run.jl is 
        # implicitly defined as true and inside set_heat is false. 
        # If you're calling your version of set_heat several times, then you're calling the version with the flag set to false, but for the convective term it has to be set to true, so that's why
        # The flag=true, the capacities are set for the convection, the flag=false they are set for the other operators

        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
       
        # printstyled(color=:red, @sprintf "\n levelset:\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

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

    # if save_length
    #     fwd.length[1] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
    # end


    ####################################################################################################
    # Initialisation
    ####################################################################################################

    HSu = get_height(grid_u,ind,dx,dy,LS[end].geoS) #TODO which LS
    HLu = get_height(grid_u,ind,dx,dy,LS[end].geoL)

    HSv = get_height(grid_v,ind,dx,dy,LS[end].geoS) #TODO which LS
    HLv = get_height(grid_v,ind,dx,dy,LS[end].geoL)

    HS = get_height(grid,ind,dx,dy,LS[end].geoS) #TODO which LS
    HL = get_height(grid,ind,dx,dy,LS[end].geoL)

    presintfc = 0.0
    # presintfc = pres0 + p_lapl ? #TODO init pressure


    #TODO perio, intfc, ... check init_fields_2!

    init_fields_2!(phS.pD,phS.p,HS,BC_pS,grid,presintfc)
    init_fields_2!(phL.pD,phL.p,HL,BC_pL,grid,presintfc)

    init_fields_2!(phS.uD,phS.u,HSu,BC_uS,grid_u,num.uD)
    # init_fields_2!(phS.ucorrD,phS.u,HSu,BC_uS,grid_u,num.uD)

    init_fields_2!(phL.uD,phL.u,HLu,BC_uL,grid_u,num.uD)
    # init_fields_2!(phL.ucorrD,phL.u,HLu,BC_uL,grid_u,num.uD)

    init_fields_2!(phS.vD,phS.v,HSv,BC_vS,grid_v,num.vD)
    # init_fields_2!(phS.vcorrD,phS.v,HSv,BC_vS,grid_v,num.vD)

    init_fields_2!(phL.vD,phL.v,HLv,BC_vL,grid_v,num.vD)
    # init_fields_2!(phL.vcorrD,phL.v,HLv,BC_vL,grid_v,num.vD)

    print("\n ","check borders",vecb(phL.vD,grid_v)[ny-1],vecb(phL.vD,grid_v)[ny],vecb(phL.vD,grid_v)[ny+1],vecb(phL.vD,grid_v)[ny+2])
    print("\n ",vecb_L(phL.vD,grid_v))
    print("\n ",vecb_B(phL.vD,grid_v))
    print("\n ",vecb_R(phL.vD,grid_v))
    print("\n ",vecb_T(phL.vD,grid_v))

    ####################################################################################################

    init_fields_2!(phS.TD,phS.T,HS,BC_TS,grid,θd)
    init_fields_2!(phL.TD,phL.T,HL,BC_TL,grid,θd)


    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check %s %s %s %s %.2e %.2e %2i\n" heat heat_convection electrolysis electrolysis_convection num.τ θd nb_transported_scalars)

        # print("\n before init",vecb_L(phL.trans_scalD[:,1], grid))

        # print("\n vecb_L(elec_condD, grid) \n ",vecb_L(phL.trans_scalD[:,2], grid) )


        for iscal=1:nb_transported_scalars

            # printstyled(color=:green, @sprintf "\n Init scal %.2e \n" concentration0[iscal])

            #Solid phase: necessary?

            # @views vec1(phS.trans_scalD[:,iscal],grid) .= vec(phS.trans_scal[:,:,iscal])
            # @views vec2(phS.trans_scalD[:,iscal],grid) .= concentration0[iscal]
            # @views init_borders!(phS.trans_scalD[:,iscal], grid, BC_trans_scal[iscal], concentration0[iscal])

            # @views vec1(phL.trans_scalD[:,iscal],grid) .= vec(phL.trans_scal[:,:,iscal])
            # @views vec2(phL.trans_scalD[:,iscal],grid) .= concentration0[iscal]
            # @views init_borders!(phL.trans_scalD[:,iscal], grid, BC_trans_scal[iscal], concentration0[iscal])

            @views phS.trans_scal[:,:,iscal] .= concentration0[iscal]
            @views phL.trans_scal[:,:,iscal] .= concentration0[iscal]

            @views init_fields_2!(phS.trans_scalD[:,iscal],phS.trans_scal[:,:,iscal],HS,BC_trans_scal[iscal],grid,concentration0[iscal])
            @views init_fields_2!(phL.trans_scalD[:,iscal],phL.trans_scal[:,:,iscal],HL,BC_trans_scal[iscal],grid,concentration0[iscal])
            
            # @views fwd.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal].*LS[end].geoL.cap[:,:,5] .+ phS.trans_scal[:,:,iscal].*LS[end].geoS.cap[:,:,5]
            # @views fwdS.trans_scalD[1,:,iscal] .= phS.trans_scalD[:,iscal]
            # @views fwdL.trans_scalD[1,:,iscal] .= phL.trans_scalD[:,iscal]
            # @views fwdL.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal]

        end


        # print("\n vecb_L(elec_condD, grid) \n ", vecb_L(phL.trans_scalD[:,2], grid) )


        # print("\n after init",vecb_L(phL.trans_scalD[:,1], grid))


        # vec1(phL.phi_eleD,grid) .= vec(phL.phi_ele)
        # vec2(phL.phi_eleD,grid) .= phi_ele0 #TODO
        # init_borders!(phL.phi_eleD, grid, BC_phi_ele, phi_ele0)

        init_fields_2!(phS.phi_eleD,phS.phi_ele,HS,BC_phi_ele,grid,phi_ele0) #
        init_fields_2!(phL.phi_eleD,phL.phi_ele,HL,BC_phi_ele,grid,phi_ele0)

        # @views fwdS.phi_eleD[1,:] .= phS.phi_eleD
        # @views fwdL.phi_eleD[1,:] .= phL.phi_eleD
    end  
    
    if electrolysis
        printstyled(color=:green, @sprintf "\n Check init_fields_2!\n")
        print_electrolysis_statistics(nb_transported_scalars,grid,phL)

        # print("\n test left H2", vecb_L(phL.trans_scalD[:,1], grid))
        # print("\n test left potential", vecb_L(phL.phi_eleD, grid))

        # plot_python_bc(num,(grid.x[:,1] - 0.5 * grid.dx[:,1])/num.plot_xscale, grid.y[:,1]/num.plot_xscale,vecb_L(phL.trans_scalD[:,1], grid),"vecb_L_H2",num.plot_prefix,grid)
        # plot_python_bc(num,(grid.x[:,1] - 0.5 * grid.dx[:,1])/num.plot_xscale, grid.y[:,1]/num.plot_xscale,vecb_L(phL.phi_eleD, grid),"vecb_L_phi",num.plot_prefix,grid)
        
        # print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
        # "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
        # "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)


    
        
        if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
            any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
            norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8 || norm(phL.trans_scal) > 1e8 || norm(phL.phi_ele) > 1e8 ||
            any(phL.trans_scal .<0))
            println(@sprintf "\n CRASHED start \n")

            # println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
            
            print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
            "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
            "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , 
            "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 ,
             "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8,"\n any(phL.trans_scal .<0): ", any(phL.trans_scal .<0))

            crashed=true
            # return nothing
            return num.current_i

        end
    end
    ####################################################################################################

   
    # @views fwdS.ucorrD[1,:,:] .= phS.ucorrD
    # @views fwdL.ucorrD[1,:,:] .= phL.ucorrD

    # @views fwdS.vcorrD[1,:,:] .= phS.vcorrD
    # @views fwdL.vcorrD[1,:,:] .= phL.vcorrD

    # @views fwdL.pD[1,:] .= phL.pD
    # @views fwdS.pD[1,:] .= phS.pD

    # @views fwdS.TD[1,:] .= phS.TD
    # @views fwdL.TD[1,:] .= phL.TD

    # @views fwd.radius[1] = current_radius


    kill_dead_cells!(phS.T, grid, LS[end].geoS)
    kill_dead_cells!(phL.T, grid, LS[end].geoL)


    #TODO check timestep coefficients n-1 
    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    # TODO kill_dead_cells! for [:,:,iscal]
    if electrolysis
        for iscal=1:nb_transported_scalars
            # @views kill_dead_cells!(phS.trans_scal[:,:,iscal], grid, LS[end].geoS) #TODO
            # @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid, LS[end].geoL) 
            # @views kill_dead_cells_val!(phS.trans_scal[:,:,iscal], grid, LS[end].geoS) #TODO
            # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[end].geoL,concentration0[iscal]) 
            @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[end].geoL,0.0) 
            @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

        end
    end     
    ####################################################################################################

    # print("\n vecb_L(elec_condD, grid) after kill \n ", vecb_L(phL.trans_scalD[:,2], grid) )


    # if nb_transported_scalars>0
    #     printstyled(color=:green, @sprintf "\n after kill \n")
    #     print_electrolysis_statistics(nb_transported_scalars,phL)
    #     printstyled(color=:green, @sprintf "\n average T %s\n" average!(phL.T, grid, LS[1].geoL, num))
    # end

    # for iLS in 1:_nLS
        # @views fwd.u[iLS,1,:,:] .= LS[iLS].u
        # @views fwd.ux[iLS,1,:,:] .= grid_u.LS[iLS].u
        # @views fwd.uy[iLS,1,:,:] .= grid_v.LS[iLS].u
        # @views fwd.κ[iLS,1,:,:] .= LS[iLS].κ
    # end
    # @views fwd.T[1,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T[:,:].*LS[end].geoS.cap[:,:,5]
    # @views fwdL.T[1,:,:] .= phL.T
    # @views fwdS.T[1,:,:] .= phS.T
    # @views fwdS.p[1,:,:] .= phS.p
    # @views fwdL.p[1,:,:] .= phL.p
    # @views fwdS.u[1,:,:] .= phS.u
    # @views fwdS.v[1,:,:] .= phS.v
    # @views fwdL.u[1,:,:] .= phL.u
    # @views fwdL.v[1,:,:] .= phL.v
    # @views fwdL.vD[1,:] .= phL.vD

    ####################################################################################################
    #Electrolysis
    ####################################################################################################
    # if electrolysis
    #     for iscal=1:nb_transported_scalars #TODO updated later cf init_fields_2!
    #         @views fwd.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal].*LS[end].geoL.cap[:,:,5] .+ phS.trans_scal[:,:,iscal].*LS[end].geoS.cap[:,:,5]
    #         @views fwdL.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal]
    #     end
        
    #     # @views fwd.mass_flux[1,:,:] .= 0.0
    #     @views fwdL.phi_ele[1,:,:] .= phL.phi_ele
    #     @views fwdL.i_current_mag[1,:,:] .= phL.i_current_mag
    #     @views fwdS.Eu[1,:,:] .= phS.Eu
    #     @views fwdS.Ev[1,:,:] .= phS.Ev
    #     @views fwdL.Eu[1,:,:] .= phL.Eu
    #     @views fwdL.Ev[1,:,:] .= phL.Ev
    # end     
    

    



    if is_FE(time_scheme) || is_CN(time_scheme)
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        # printstyled(color=:red, @sprintf "\n levelset 2:\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

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
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb
            AϕS = spzeros(nt, nt)
            AϕL = spzeros(nt, nt)
            if electrolysis 
                Aphi_eleL = spzeros(nt, nt)
            end

            ATS = spzeros(nt, nt)
            ATL = spzeros(nt, nt)
            BTS = spzeros(nt, nt)
            BTL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny
            nt = (num.nLS + 1) * ni + nb
            AuS = spzeros(nt, nt)
            AuL = spzeros(nt, nt)
            BuS = spzeros(nt, nt)
            BuL = spzeros(nt, nt)

            ni = grid_v.nx * grid_v.ny
            nb = 2 * grid_v.nx + 2 * grid_v.ny
            nt = (num.nLS + 1) * ni + nb
            AvS = spzeros(nt, nt)
            AvL = spzeros(nt, nt)
            BvS = spzeros(nt, nt)
            BvL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny + grid_v.nx * grid_v.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny + 2 * grid_v.nx + 2 * grid_v.ny
            nt = (nLS - nNavier + 1) * ni + nNavier * grid.nx * grid.ny + nb
            AuvS = spzeros(nt, nt)
            AuvL = spzeros(nt, nt)
            BuvS = spzeros(nt, nt)
            BuvL = spzeros(nt, nt)

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uS,
                    AuS, BuS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vS,
                    AvS, BvS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pS, opC_uS, opC_vS,
                    AuvS, BuvS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pS, opC_uS, opC_vS,
                AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                ATS, BTS,
                opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )

            # if electrolysis_solid_phase
                # set_heat!(
                #     BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                #     ATS, BTS,
                #     opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                #     periodic_x, periodic_y, heat_convection, true, BC_int
                # )    

                # set_scalar_transport!(BC_trans_scal[iscal].int, num, grid, opC_TS, LS[1].geoS, phS, concentration0[iscal], BC_trans_scal[iscal],
                #                                             LS[1].MIXED, LS[1].geoS.projection,
                #                                             ATL, BTL,
                #                                             opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                #                                             periodic_x, periodic_y, electrolysis_convection, true, diffusion_coeff[iscal])
            # end 

            if test_laplacian
                exact_laplacian = test_laplacian_pressure(num,grid_v,phL,opC_pL, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L)
                return exact_laplacian
            end

            if !navier
                _ = FE_set_momentum(
                    BC_int, num, grid_u, opC_uL,
                    AuL, BuL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    true
                )
                _ = FE_set_momentum(
                    BC_int, num, grid_v, opC_vL,
                    AvL, BvL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    opC_pL, opC_uL, opC_vL,
                    AuvL, BuvL,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true
                )
            end
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            _ = set_poisson(
                BC_int, num, grid, a0_p, opC_pL, opC_uL, opC_vL,
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )

            set_heat!(
                BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                ATL, BTL,
                opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )
            # if electrolysis
            #     for iscal=1:nb_transported_scalars

            #         set_scalar_transport!(BC_trans_scal[iscal].int, num, grid, opC_TL, LS[1].geoL, phL, concentration0[iscal], BC_trans_scal[iscal],
            #                                                 LS[1].MIXED, LS[1].geoL.projection,
            #                                                 ATL,BTL,
            #                                                 opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
            #                                                 periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff[iscal])
            #     end
            # end 
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    if heat_convection || electrolysis_convection
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
        # printstyled(color=:red, @sprintf "\n levelset 3:\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])
    end

    V0S = volume(LS[end].geoS)
    V0L = volume(LS[end].geoL)



    current_t = 0.

    while num.current_i < max_iterations + 1

        ####################################################################################################
        #Adapt timestep
        # printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL num.τ)
        if adapt_timestep_mode !=0
            num.τ = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
            # print("after adapt_timestep!")
            printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e num.τ : %.2e\n" CFL num.τ num.τ)
        end
        ####################################################################################################

        printstyled(color=:red, @sprintf "\n iter: %5i\n" num.current_i)
        println(grid.LS[1].geoL.dcap[1,1,:])

        if electrolysis

            ####################################################################################################
            #PDI (IO)
            ####################################################################################################
            
            #TODO not necessary to expose everything now for ex only LS ? and the rest later

            if num.io_pdi>0

                try
                    printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                    time = current_t #Cdouble
                    nstep = num.current_i
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
            
                    Eus,Evs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.Eu, phL.Ev)
            
                    us,vs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.u,phL.v)
            
                    # print("\n before write \n ")
            
                    iLSpdi = 1 # all LS iLS = 1 # or all LS ?


                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name

                    @ccall "libpdi".PDI_multi_expose("write_data_start_loop"::Cstring,
                    "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,          
                    C_NULL::Ptr{Cvoid})::Cvoid
            
                    # print("\n after write \n ")
            
                    # @ccall "libpdi".PDI_finalize()::Cvoid
            
                    # printstyled(color=:red, @sprintf "\n PDI test end\n" )
            
                catch error
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                    print(error)
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                end

               

            end #if io_pdi

            ####################################################################################################

        end


        if !stefan
            V .= speed*ones(ny, nx)
        end

        # Solve heat equation
        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                    ATS, BTS,
                    opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                rhs = set_heat!(
                    BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                    ATL, BTL,
                    opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs
                phL.T .= reshape(veci(phL.TD,grid,1), grid)

            end
        end

        ####################################################################################################
        #Electrolysis
        ####################################################################################################  
        if electrolysis
            # if electrolysis_solid_phase

            #     for iscal=1:nb_transported_scalars
            #         @views kill_dead_cells!(phS.trans_scal[:,:,iscal], grid, LS[1].geoS)
            #         @views veci(phS.trans_scalD[:,iscal],grid,1) .= vec(phS.trans_scal[:,:,iscal])

            #         rhs = set_scalar_transport!(BC_trans_scal[iscal].int, num, grid, opC_TS, LS[1].geoS, phS, concentration0[iscal], BC_trans_scal[iscal],
            #                                             LS[1].MIXED, LS[1].geoS.projection,
            #                                             ATL, BTL,
            #                                             opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
            #                                             periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff[iscal])
            #         mul!(rhs, BTL, phS.trans_scalD[:,iscal], 1.0, 1.0)

            #         @views phS.trans_scalD[:,iscal] .= ATL \ rhs
            #         @views phS.trans_scal[:,:,iscal] .= reshape(veci(phS.trans_scalD[:,iscal],grid,1), grid)
            #     end
            # end

            if electrolysis_liquid_phase

                # printstyled(color=:green, @sprintf "\n Before transport\n")
                # print_electrolysis_statistics(nb_transported_scalars,grid,phL)

                # print("\n before res",vecb_L(phL.trans_scalD[:,1], grid))

                #TODO check BC not overwritten by different scalars
                #TODO check ls advection true when n scalars

                # print("\n vecb_L",vecb_L(phL.phi_eleD, grid))
                # print("\n phL.phi_ele[:,1]",phL.phi_ele[:,1])

                #TODO phL.phi_ele[:,1] or vecb_L(phL.phi_eleD, grid)
                # we suppose phi(x=0)=... cf Khalighi
                #but here BC
                
                if electrolysis_reaction == "Butler_no_concentration"

                    if heat
                        i_butler = butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,vecb_L(phL.phi_eleD, grid),
                        phi_ele1,Ru,phL.T)
                    else
                        i_butler = butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,vecb_L(phL.phi_eleD, grid),
                        phi_ele1,Ru,temperature0)
                            
                    end   
                end
                # print("\n i_butler",i_butler)

                #################################################

                
                # phL.saved_scal[:,:,1],phL.saved_scal[:,:,2],phL.saved_scal[:,:,3],phL.saved_scal[:,:,4]=compute_mass_flux!(num,grid, grid_u, 
                # grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,1)

                # for iscal=1:nb_saved_scalars
                #     @views fwd.saved_scal[1,:,:,iscal] .= phL.saved_scal[:,:,iscal] #varflux[iscal] #reshape(varfluxH2, grid)
                # end

                #################################################


                # print("\n vecb_L(elec_condD, grid) before res \n ", vecb_L(phL.trans_scalD[:,2], grid) )

                ####################################################################################################
                # New start scalar loop
                ####################################################################################################

                if imposed_velocity == "zero"
                    phL.u .= 0.0
                    phL.v .= 0.0
                elseif imposed_velocity == "constant"
                        
                        #Required to modify whole uD vD
                        phL.u .= 0.0
                        phL.v .= BC_vL.bottom.val    
                        phL.uD .= 0.0
                        phL.vD .= BC_vL.bottom.val    

                        phL.pD .= 0.0

                        if ((num.current_i-1)%show_every == 0) 
                            printstyled(color=:red, @sprintf "\n Imposed velocity v min %.2e max %.2e\n" minimum(phL.vD) maximum(phL.vD))
                            printstyled(color=:red, @sprintf "\n Imposed velocity u min %.2e max %.2e\n" minimum(phL.uD) maximum(phL.uD))
                        end

                elseif imposed_velocity == "Poiseuille"

                    vPoiseuille = Poiseuille_fmax.(grid_v.x,num.v_inlet,num.L0)

                    
                    #Required to modify whole uD vD
                    phL.u .= 0.0
                    phL.v .= vPoiseuille 
 
                    phL.uD .= 0.0 #Dirichlet and Neumann
                    # phL.vD .= BC_vL.bottom.val    
                    # vecb...
                    # ....
                    
                    if ((num.current_i-1)%show_every == 0) 
                        printstyled(color=:red, @sprintf "\n Imposed velocity v min %.2e max %.2e\n" minimum(phL.vD) maximum(phL.vD))
                        printstyled(color=:red, @sprintf "\n Imposed velocity u min %.2e max %.2e\n" minimum(phL.uD) maximum(phL.uD))       
                        
                        printstyled(color=:cyan, @sprintf "\n before scalar transport 0\n")

                        print_electrolysis_statistics(nb_transported_scalars,grid,phL)
                    end 
                    
                elseif imposed_velocity == "radial"
                    
                    phL.v .= 0.0
                    phL.u .= 0.0
                    phS.v .= 0.0
                    phS.u .= 0.0
                
                    for ju in 1:grid_u.ny
                        for iu in 1:grid_u.nx
                            xcell = grid_u.x[ju,iu]
                            ycell = grid_u.y[ju,iu]
                
                            vec0 = [num.xcoord, num.ycoord]
                            vec1 = [xcell, ycell]
                
                            vecr = vec1-vec0
                            normr = norm(vecr)
                            vecr .*= 1.0/normr
                            factor = 1.0/normr 
                            factor *= num.radial_vel_factor
                            # if normr>radius                              
                            #     phL.u[ju,iu] = factor * vecr[1]
                            #     phS.u[ju,iu] = factor * vecr[1]                
                            # end
                            phL.u[ju,iu] = factor * vecr[1]
                            phS.u[ju,iu] = factor * vecr[1]   
                        end
                    end
                
                    for jv in 1:grid_v.ny
                        for iv in 1:grid_v.nx    
                            xcell = grid_v.x[jv,iv]
                            ycell = grid_v.y[jv,iv]
                
                            vec0 = [num.xcoord, num.ycoord]
                            vec1 = [xcell, ycell]
                
                            vecr = vec1-vec0
                            normr = norm(vecr)
                            vecr .*= 1.0/normr
                            factor = 1.0/normr 
                            factor *= num.radial_vel_factor

                            # if normr>radius           
                            #     phL.v[jv,iv] = factor * vecr[2]
                            #     phS.v[jv,iv] = factor * vecr[2]                
                            # end
                            #print("\n i ",iv," j ",jv," vec",vecr[2]," y ",ycell," ycoord ",ycoord," v ",phL.v[jv,iv])
                            phL.v[jv,iv] = factor * vecr[2]
                            phS.v[jv,iv] = factor * vecr[2]   

                        end
                    end

                    # grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (current_radius) * ones(ny, nx)                  



                end #imposed_velocity

                


                for iscal=1:nb_transported_scalars
                    @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL,0.0) 
                    # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 
                    @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

                    if electrolysis_reaction == "Butler_no_concentration"


                        BC_trans_scal[iscal].left.val = i_butler./(2*Faraday*diffusion_coeff[iscal])

                        # printstyled(color=:red, @sprintf "\n Butler deactivated \n")
                        # if iscal==1 || iscal==2
                        #     BC_trans_scal[iscal].left.val = i_butler./(2*Faraday*diffusion_coeff[iscal])
                        # end

                        if iscal==1 || iscal==2
                            BC_trans_scal[iscal].left.val .*=-1 #H2O consummed
                        end
                        # print("\n")
                        # print("\n left BC ", BC_trans_scal[iscal].left.val)

                        # for testn in 1:ny
                        #     printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i border %.5e\n" testn ny-testn+1 vecb_L(phL.trans_scalD[:,iscal], grid)[testn])
                        # end

                    end
                end

                #TODO convection_Cdivu BC divergence
                #TODO check nb_transported_scalars>1

                if ((num.current_i-1)%show_every == 0) 
                    printstyled(color=:cyan, @sprintf "\n before scalar transport \n")
                    print_electrolysis_statistics(nb_transported_scalars,grid,phL)
                end
             
                # printstyled(color=:red, @sprintf "\n levelset: before scalar_transport\n")
                # println(grid.LS[1].geoL.dcap[1,1,:])

                if electrolysis_advection
                    update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
                end
                scalar_transport!(BC_trans_scal, num, grid, opC_TL, LS[1].geoL, phL, concentration0,
                LS[1].MIXED, LS[1].geoL.projection, opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff,convection_Cdivu)

            

                # scalar_transport_2!(BC_trans_scal, num, grid, opC_TL, LS[1].geoL, phL, concentration0,
                # LS[1].MIXED, LS[1].geoL.projection, opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                # periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff,convection_Cdivu)

                if ((num.current_i-1)%show_every == 0) 
                    printstyled(color=:cyan, @sprintf "\n after scalar transport \n")
                    print_electrolysis_statistics(nb_transported_scalars,grid,phL)
                end

                if imposed_velocity != "none" || num.debug== "scalar_testing"
                    scal_error=0.0
                    for iscal in 1:nb_transported_scalars

                        # print("\n maximum ",maximum(phL.trans_scalD[:,iscal]), )
                        # printstyled(color=:cyan, @sprintf "\n error after scalar transport max %.2e min %.2e\n" maximum(phL.trans_scalD[:,iscal]) minimum(phL.trans_scalD[:,iscal]))

                        scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-concentration0[iscal])./concentration0[iscal])
                        scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid).-concentration0[iscal])./concentration0[iscal])
                        scal_error = max(scal_error_bulk,scal_error_border,scal_error)

                    end

                    printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e CFL %.2e\n" scal_error num.v_inlet*num.dt0/grid.dx[1,1])

                    printstyled(color=:red, @sprintf "\n Poiseuille \n")

                    # Check the velocity field before the scalar transport
                    test_Poiseuille(num,phL.vD,grid_v)

                    printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

                    printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

                    printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

                    compute_grad_p!(num,grid, grid_u, grid_v, phL.pD, opC_pL, opC_uL, opC_vL)

                end

                # if imposed_velocity =="none"
                #     printstyled(color=:red, @sprintf "\n after scalar transport \n")

                #     # Check the velocity field before the scalar transport
                #     test_Poiseuille(num,phL,grid_v)
                    
                # end
                ####################################################################################################
                # New start scalar loop
                ####################################################################################################


                # ####################################################################################################
                # # Start scalar loop
                # ####################################################################################################

                # for iscal=1:nb_transported_scalars


                #     # print("before trans")

                #     # phL.saved_scal[:,:,1],phL.saved_scal[:,:,2],phL.saved_scal[:,:,3],phL.saved_scal[:,:,4]=compute_mass_flux!(num,grid, grid_u, 
                #     # grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,1)
            
                #     # for iscal=1:nb_saved_scalars
                #     #     @views fwd.saved_scal[1,:,:,iscal] .= phL.saved_scal[:,:,iscal] #varflux[iscal] #reshape(varfluxH2, grid)
                #     # end

                #     # print("\n vecb_L(elec_condD, grid) before BC \n ", vecb_L(phL.trans_scalD[:,2], grid) )

                #     # print("before trans")
                #     if electrolysis_reaction == "Butler_no_concentration"
                #         BC_trans_scal[iscal].left.val = i_butler./(2*Faraday*diffusion_coeff[iscal])

                #         # print("\n Butler",BC_trans_scal[iscal].left.val)

                #         if iscal==1 || iscal==2
                #             BC_trans_scal[iscal].left.val .*=-1 #H2O
                #         end
                #     end
                  
                 
                #     # print("\n vecb_L(elec_condD, grid) after BC \n ", vecb_L(phL.trans_scalD[:,2], grid) )


                #     # @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL)
                #     @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

                #     @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

                #     # print("\n",BC_trans_scal[iscal].int, " ",maximum(veci(phL.trans_scalD[:,iscal],grid,2)))
                #     # print("\n test chi",maximum(opC_TL.χ[1])," ",maximum(vec2(opC_TL.χ[1],grid)))

                #     # phL.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)




                #     diffusion_coeff_iscal = diffusion_coeff[iscal]

                #     # diffusion_coeff_iscal = 0.0
                #     # printstyled(color=:red, @sprintf "\n 0 diffusion \n")

                #     # print("BC top ",BC_trans_scal[iscal].top)
                #     # print("BC bottom ",BC_trans_scal[iscal].bottom)

                #     # #########################################################################"

                #     printstyled(color=:red, @sprintf "\n levelset: before set_scalar_transport!\n")
                #     println(grid.LS[1].geoL.dcap[1,1,:])

                #     rhs = set_scalar_transport!(BC_trans_scal[iscal].int, num, grid, opC_TL, LS[1].geoL, phL, concentration0[iscal], BC_trans_scal[iscal],
                #                                         LS[1].MIXED, LS[1].geoL.projection,
                #                                         ATL,BTL,
                #                                         opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                #                                         periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff_iscal)
                    
                #     ##########################################################################################################"

                #     printstyled(color=:red, @sprintf "\n levelset: after set_scalar_transport!\n")
                #     println(grid.LS[1].geoL.dcap[1,1,:])

                #     if convection_Cdivu
                #         # Duv = fzeros(grid)
                #         # Duv = fnzeros(grid,num)

                #         Duv = opC_pL.AxT * vec1(phL.uD,grid_u) .+ opC_pL.Gx_b * vecb(phL.uD,grid_u) .+
                #         opC_pL.AyT * vec1(phL.vD,grid_v) .+ opC_pL.Gy_b * vecb(phL.vD,grid_v)
                #         for iLS in 1:nLS
                #             if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS])
                #                 Duv .+= opC_pL.Gx[iLS] * veci(phL.uD,grid_u,iLS+1) .+ 
                #                         opC_pL.Gy[iLS] * veci(phL.vD,grid_v,iLS+1)
                #             end
                #         end
                    
                #         # rhs .+= Duv #.* phL.trans_scalD[:,iscal] multiplied just after

                #         # vec1(rhs_ϕ,grid) .= rho1 .* iτ .* Duv #TODO
                #         vec1(rhs,grid) .+= Duv 

                #         printstyled(color=:green, @sprintf "\n max Duv for C.div(u): %.2e\n" maximum(Duv))


                #     end
                #     ##########################################################################################################"
                    
                #     @views mul!(rhs, BTL, phL.trans_scalD[:,iscal], 1.0, 1.0) #TODO @views not necessary ?


                #     # if nb_saved_scalars>1
                #     #     # phL.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
                #     #     phL.saved_scal[:,:,2]=reshape(veci(rhs,grid,1), grid)
                #     # end

                #     if nb_saved_scalars>4
                #         # phL.saved_scal[:,:,5]=reshape(opC_TL.χ[1].diag,grid)
                #         phL.saved_scal[:,:,5]=reshape(veci(rhs,grid,1), grid)

                #         if nb_saved_scalars>5
                #             phL.saved_scal[:,:,6]=reshape(opC_TL.χ[1].diag,grid)
                #         end
                #     end


                #     # print("\n test left before A/r L", vecb_L(phL.trans_scalD[:,iscal], grid))
                #     # print("\n test left before A/r T", vecb_T(phL.trans_scalD[:,iscal], grid))


                #     @views phL.trans_scalD[:,iscal] .= ATL \ rhs



                #     # print("\n test left after A/r L", vecb_L(phL.trans_scalD[:,iscal], grid))
                #     # print("\n test left after A/r T", vecb_T(phL.trans_scalD[:,iscal], grid))

                #     # @views phL.trans_scal[:,:,iscal] .= reshape(veci(phL.trans_scalD[:,iscal],grid,2), grid)


                #     # printstyled(color=:green, @sprintf "\n average c %s\n" average!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL, num))


                #     @views phL.trans_scal[:,:,iscal] .= reshape(veci(phL.trans_scalD[:,iscal],grid,1), grid)

                #     printstyled(color=:cyan, @sprintf "\n after resol\n")
                #     print("\n",phL.trans_scal[1,:,iscal])
                #     print("\n",reshape(veci(rhs,grid,1), grid)[1,:])
                #     print("\n",reshape(veci(rhs,grid,1), grid)[2,:])
                #     print("\n",reshape(veci(rhs,grid,1), grid)[3,:])


                    


                #     nonzero = veci(phL.trans_scalD[:,iscal],grid,2)[abs.(veci(phL.trans_scalD[:,iscal],grid,2)) .> 0.0]
                #     # print("nonzero\n")
                #     # print(nonzero)
                #     print("\n mean ",mean(nonzero))

                #     if iscal!=2 #H2O consummed at the electrode, would need to make distinction to make sure the decrease in H2O is physical or not
                      
                      
                #         @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, LS[1].geoL,concentration0[iscal]) 

                #         # @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

                #         if any(phL.trans_scal[:,:,iscal].<concentration0[iscal]*(1-num.concentration_check_factor))
                #             print("iscal ",iscal)
                #             printstyled(color=:red, @sprintf "\n concentration: %.2e %.2e \n" minimum(phL.trans_scal[:,:,iscal]) concentration0[iscal]*(1-num.concentration_check_factor))
                #             printstyled(color=:red, @sprintf "\n concentration drop: %.2e%% \n" (minimum(phL.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]*100)
                #             @error("concentration too low")
                #             return num.current_i
                #         end
                #     end

                #     # print("\n vecb_L(elec_condD, grid) after res 0 \n ", vecb_L(phL.trans_scalD[:,2], grid) )


                #     # print("all\n")
                #     # print(phL.trans_scalD[:,iscal])

                #     # print("\n test", phL.trans_scal[:,1,iscal])
                #     # print("\n test", phL.trans_scal[1,:,iscal])

                #     # print("\n test left", vecb_L(phL.trans_scalD[:,iscal], grid))
                #     # print("\n test right", vecb_R(phL.trans_scalD[:,iscal], grid))
                #     print("\n test bottom", vecb_B(phL.trans_scalD[:,iscal], grid))
                #     # print("\n test top", vecb_T(phL.trans_scalD[:,iscal], grid))

                #     print("\n test v bottom", maximum(vecb_B(phL.vD, grid)))
                #     print("\n test v top", maximum(vecb_T(phL.vD, grid)))
                #     print("\n test v 1 ", maximum(phL.v[1,:]))
                #     print("\n test v 2 ", maximum(phL.v[1,:]))

                #     print("\n test concentration ", minimum(phL.trans_scal[:,:,iscal]), " ", maximum(phL.trans_scal[:,:,iscal]))


                #     # print("\n test diff", phL.v[1,:] .- vecb_T(phL.vD, grid))



                #     # for jplot in 1:ny
                #     #     for iplot in 1:nx

                #             # II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                #     #         pII = lexicographic(II, grid.ny)

                #     #         if phL.trans_scalD[pII,iscal] < 0.0
                #     #         # if phL.trans_scalD[pII,iscal] < concentration0[iscal]

                #     #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e %.2e %.2e %.2e \n" iplot jplot grid.x[iplot]/num.plot_xscale grid.y[jplot]/num.plot_xscale phL.trans_scalD[pII,iscal] rhs[pII])
                        
                #     #         end
                #     #     end
                #     # end


                #     print("end scal",iscal)

                #     printstyled(color=:red, @sprintf "\n levelset: end scal set_scalar_transport!\n")
                #     println(grid.LS[1].geoL.dcap[1,1,:])

                # end
                # ####################################################################################################
                # # End scalar loop
                # ####################################################################################################




                # print("\n vecb_L(elec_condD, grid) after res 1 \n ", vecb_L(phL.trans_scalD[:,2], grid) )


                # print("\n after res",vecb_L(phL.trans_scalD[:,1], grid))


                # printstyled(color=:green, @sprintf "\n After transport\n")
                # print_electrolysis_statistics(nb_transported_scalars,grid,phL)


                #TODO heat for electrolysis

                ####################################################################################################
                #Electrolysis: Poisson
                ####################################################################################################  
                if electric_potential
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

                    # #print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")

                    if electrolysis && nb_transported_scalars>1 && conductivity_mode != 0
                        if heat 
                            elec_cond  = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru.*phL.T) 
                            elec_condD = 2*Faraday^2 .*phL.trans_scalD[:,2].*diffusion_coeff[2]./(Ru.*phL.TD)
                        else
                            elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru*temperature0) 
                            elec_condD = 2*Faraday^2 .*phL.trans_scalD[:,2].*diffusion_coeff[2]./(Ru.*temperature0)

                        end
                    else 
                        elec_cond = ones(grid)
                        elec_condD = fnones(grid, num)
                        printstyled(color=:green, @sprintf "\n conductivity one")
                        # b_phi_ele = zeros(grid)

                    end 


                    # for iLS in 1:nLS
                    #     kill_dead_cells_val!(veci(elec_condD,grid,iLS+1), grid, LS[iLS].geoL,1.0)
                    # end

                    # [LS[iLS].geoS for iLS in 1:_nLS]

                    # for iLS in 1:nLS
                    #     kill_dead_cells_val!(veci(elec_condD,grid,iLS+1), grid, LS[1].geoL,1.0)
                    # end

                    # kill_dead_cells_val_wall!(vecb(elec_condD,grid), grid, LS[1].geoL,1.0)



            
                    
                    # iLS=1
                    # kill_dead_cells_val!(veci(elec_condD,grid,iLS+1), grid, LS[1].geoL,1.0)

                    # kill_dead_cells_val!(elec_condD, grid, LS[1].geoL,1.0) #overwrite null conductivity in solid cells for BC 1/cond

                    # print("\n \n elec_condD",vecb_L(elec_condD, grid))


                    #TODO Poisson with variable coefficients
                    #TODO need to iterate? since nonlinear

                    #TODO no concentration prefactor

                    #Update Butler-Volmer Boundary Condition with new potential 
                
                    # eta = phi_ele1 .- phL.phi_ele[:,1]
                    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))

                    if occursin("Butler",electrolysis_reaction)

                        # For small cells
                        if bulk_conductivity == 0
                            BC_phi_ele.left.val .= -i_butler./vecb_L(elec_condD, grid)
                        elseif bulk_conductivity == 1
                            # Recommended as long as cell merging not implemented:
                            # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
                            # which produces higher concentration in front of the contact line
                            BC_phi_ele.left.val .= -i_butler./elec_cond[:,1]
                        elseif bulk_conductivity == 2
                            BC_phi_ele.left.val .= -i_butler./vecb_L(elec_condD, grid)

                            iLS = 1 #TODO end ? if several LS ?
                            for j in 1:grid.ny
                                if vecb_L(grid.LS[iLS].geoL.cap[:,5],grid)[j] < ϵ
                                    BC_phi_ele.left.val[j] = -i_butler[j]/elec_cond[j,1] 
                                end
                            end

                        end

                        # if heat
                        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,phL.T)./elec_cond[:,1]
                        # else
                        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,temperature0)./elec_cond[:,1]
                            
                        #     # for iscal=1:nb_transported_scalars
                        #     #     BC_trans_scal[iscal].left.val = butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,temperature0)./(2*Faraday*diffusion_coeff[iscal])
                        #     #     if iscal==1 || iscal==2
                        #     #         BC_trans_scal[iscal].left.val .*=-1 #H2O
                        #     #     end
                        #     # end
                        # end    

                    # elseif electrolysis_reaction == ""
                    #     # BC_phi_ele.left.val = -butler_volmer_concentration.(alpha_a,alpha_c,Faraday,i0,phL.phi_ele[:,1],phi_ele1,Ru,temperature0)./elec_cond
                    
                        
                        # print("\n before ",BC_phi_ele.left.val)

                        # print("\n vecb_L(elec_condD, grid)",vecb_L(elec_condD, grid) )

                        # print("\n vecb_L(elec_condD, grid) after res \n ", vecb_L(phL.trans_scalD[:,2], grid) )


                        # for iLS in 1:nLS
                        #     # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
                        #     for i = 1:grid.ny
                        #         # print("vecb cap",vecb_L(grid.LS[iLS].geoL.cap[:,5],grid))
                        #         II = CartesianIndex(i,1)
                        #         if grid.LS[iLS].geoL.cap[II,5] < 1e-12
                        #             BC_phi_ele.left.val[i] = 1.0
                        #         end
                        #     end
                        # end


                        # grid.LS[iLS].geoL.cap[II,1:4]

                        for iLS in 1:nLS
                            # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
                            for i = 1:grid.ny
                                # print("vecb cap",vecb_L(grid.LS[iLS].geoL.cap[:,5],grid))
                                II = CartesianIndex(i,1)
                                if grid.LS[iLS].geoL.cap[II,1] < 1e-12
                                    BC_phi_ele.left.val[i] = 1.0
                                end
                            end
                        end

                        # print("\n after ",BC_phi_ele.left.val)


                    end

                    #debugphi

                    # print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                    # "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                    # "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
        

                    #TODO kill_dead_cells! ?
                    kill_dead_cells!(phL.phi_ele, grid, LS[1].geoL)
                    veci(phL.phi_eleD,grid,1) .= vec(phL.phi_ele)
                    #################################################################"

                    # print("\n elec_condD: ",any(isnan, elec_condD),"\n BC_phi_ele.left.val: ",any(isnan, BC_phi_ele.left.val),"\n")

                    # printstyled(color=:green, @sprintf "\n BC phi : %.2e \n" BC_phi_ele.int.val)

                    #TODO BC several LS
                    #Poisson with variable coefficient
                    rhs_phi_ele = set_poisson_variable_coeff(
                        [BC_phi_ele.int], num, grid, grid_u, grid_v, a0_p, opC_pL, opC_uL, opC_vL,
                        Aphi_eleL, 
                        # elec_Lpm1_L, elec_bc_Lpm1_L, elec_bc_Lpm1_b_L, 
                        BC_phi_ele,
                        true,elec_condD
                    )

                    # print("\n Aphi_eleL: ",any(isnan, Aphi_eleL),"\n rhs_phi_ele: ",any(isnan, rhs_phi_ele),"\n")

                    # print("\n veci rhs_phi_ele 2: ",any(isnan, veci(rhs_phi_ele,grid,1)),"\n veci rhs_phi_ele 1 : ",any(isnan, veci(rhs_phi_ele,grid,1)),"\n")
            

                    # print("\n \n elec_condD",vecb_L(elec_condD, grid))
                    # print("\n \n BC_phi_ele.left.val",BC_phi_ele.left.val)

                    

                    # print("\n \n vecb_L",vecb_L(rhs_phi_ele, grid))

                    # print("\n \n vecb_R",vecb_R(rhs_phi_ele, grid))

                    # print("\n \n vecb_T",vecb_T(rhs_phi_ele, grid))

                    # print("\n \n vecb_B",vecb_B(rhs_phi_ele, grid))

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

                    @time res_phi_ele .= Aphi_eleL \ rhs_phi_ele

                    #TODO or use mul!(rhs, BTL, phL.TD, 1.0, 1.0) like in :

                    # kill_dead_cells!(phL.T, grid, LS[1].geoL)
                    # veci(phL.TD,grid,1) .= vec(phL.T)
                    # rhs = set_heat!(
                    #     BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                    #     ATL, BTL,
                    #     opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    #     periodic_x, periodic_y, heat_convection, advection, BC_int
                    # )
                    # mul!(rhs, BTL, phL.TD, 1.0, 1.0)

                    # phL.TD .= ATL \ rhs
                    # phL.T .= reshape(veci(phL.TD,grid,1), grid)


                    # phL.phi_eleD .= Aphi_eleL \ rhs_phi_ele
                    phL.phi_eleD .= res_phi_ele

                    phL.phi_ele .= reshape(veci(phL.phi_eleD,grid,1), grid)




                    if any(isnan, phL.phi_eleD)
                        print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                        "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                        "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
            


                        print("\n phL.phi_eleD: ",any(isnan, phL.phi_eleD),"\n phL.phi_ele: ",any(isnan, phL.phi_ele),"\n")

                        print("\n Aphi_eleL: ",any(isnan, Aphi_eleL),"\n rhs_phi_ele: ",any(isnan, rhs_phi_ele),"\n")


                
                        
                        print("\n \n vecb_L",vecb_L(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_R",vecb_R(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_T",vecb_T(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_B",vecb_B(phL.phi_eleD[:,1], grid))

                        plot_python_bc(num,(grid.x[:,1] - 0.5 * grid.dx[:,1])/num.plot_xscale, grid.y[:,1]/num.plot_xscale,vecb_L(phL.phi_eleD, grid),"vecb_L_phi_after_res",num.plot_prefix,grid)
                


                        # Print values on line
                        # for iplot in 1:ny
                        #     II = CartesianIndex(iplot, 1) #(id_y, id_x)
                        #     pII = lexicographic(II, grid.ny)
        
                        #     # printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e \n" iplot grid.y[iplot]/num.plot_xscale mass_flux_vec1[pII] mass_flux_vecb[pII] mass_flux_veci[pII])
                        #     printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e \n" iplot grid.y[iplot]/num.plot_xscale mass_flux_vec1[pII] mass_flux_vecb[pII] mass_flux_veci[pII])
                        # end

                        vecbphi = vecb_L(phL.phi_eleD[:,1], grid)
                        # Print values on line
                        id_x = 1

                        for iplot in 1:ny
                            II = CartesianIndex(iplot, id_x) #(id_y, id_x)
                            # pII = lexicographic(II, grid.ny)
        
                            # printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e \n" iplot grid.y[iplot]/num.plot_xscale mass_flux_vec1[pII] mass_flux_vecb[pII] mass_flux_veci[pII])
                            printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e %.2e\n" iplot grid.y[iplot]/num.plot_xscale vecbphi[iplot] grid.LS[1].u[iplot,id_x] grid.LS[1].geoL.cap[II,5] phL.trans_scal[iplot,id_x,1])
                    
                        end

                        #############################################################################################################
                        #Print values on line
                        id_x = 1
                        for iplot in 1:ny
                            printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n" iplot grid.y[iplot]/num.plot_xscale phL.saved_scal[iplot,id_x,2] phL.saved_scal[iplot,id_x,3] phL.saved_scal[iplot,id_x,4] phL.trans_scal[iplot,id_x,1] grid.LS[1].u[iplot,id_x] phL.trans_scalD[iplot,id_x,1])
                        end
                        #############################################################################################################


                    end

                
                    # TODO compute magnitude of exchange current
                    # gradient!(::Neumann, Ox, Oy, Bx, By, HNx, HNy, Divx, Divy, dcap, n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
                    # TODO add post-treatment variables

                    #TODO update BC concentration

                    # compute_grad_phi_ele!(grid, phL, V, periodic_x, periodic_y) #TODO current
                    compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, opC_pL, opC_pS) #TODO current

                    # scal_magnitude

                    if electrolysis && nb_transported_scalars>1
                        if heat 
                            elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru.*phL.T) #phL.T
                        else
                            elec_cond = 2*Faraday^2 .*phL.trans_scal[:,:,2].*diffusion_coeff[2]./(Ru*temperature0) 
                        end
                    else 
                        elec_cond = ones(grid)
                        printstyled(color=:green, @sprintf "\n conductivity one")

                    end 

                    phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude



                    ####################################################################################################

                    if nb_transported_scalars == 0
                        phL.saved_scal[:,:,1],phL.saved_scal[:,:,2],phL.saved_scal[:,:,3],phL.saved_scal[:,:,4]=compute_mass_flux!(num,grid, grid_u, 
                        grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,0,LS[1].geoL)
                        print()
                        printstyled(color=:green, @sprintf "\n Flux: %.2e \n" -diffusion_coeff[1] * sum(phL.saved_scal[:,:,1]))

                    end

                end #electric_potential

            end
        end
        ####################################################################################################
        
        #    grid.LS[i].α  which is the angle of the outward point normal with respect to the horizontal axis

        for iLS in 1:nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
            elseif is_fs(BC_int[iLS])
                V .= 0.0
                printstyled(color=:green, @sprintf "\n V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))

                if electrolysis_phase_change_case!="none"                   
                    if electrolysis_phase_change_case == "levelset"

                        # plot_electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)

                        update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y, Vmean, phL.trans_scalD[:,1],diffusion_coeff[1],concentration0[1])
                        # ,opC_pL)
                    end
                    # update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
                    printstyled(color=:green, @sprintf "\n V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))

                else
                    update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
                end

            
            elseif (electrolysis && occursin("Khalighi",electrolysis_phase_change_case))

                if ((num.current_i-1)%show_every == 0) 
                    print_electrolysis_statistics(nb_transported_scalars,grid,phL)
                end

                previous_radius = current_radius


                # print("\n test left H2", vecb_L(phL.trans_scalD[:,1], grid))
                # print("\n test left potential", vecb_L(phL.phi_eleD, grid))

                # plot_python_bc(num,(grid.x[:,1] - 0.5 * grid.dx[:,1])/num.plot_xscale, grid.y[:,1]/num.plot_xscale,vecb_L(phL.trans_scalD[:,1], grid),"vecb_L_H2_2",num.plot_prefix,grid)
                # plot_python_bc(num,(grid.x[:,1] - 0.5 * grid.dx[:,1])/num.plot_xscale, grid.y[:,1]/num.plot_xscale,vecb_L(phL.phi_eleD, grid),"vecb_L_phi_2",num.plot_prefix,grid)
        


                # print("\n test left H2", vecb_L(phL.trans_scalD[:,1], grid))

                if nb_saved_scalars==1
                    phL.saved_scal[:,:,1],_,_,_=compute_mass_flux!(num,grid, grid_u, 
                    grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,1,LS[1].geoL)
                else
                    phL.saved_scal[:,:,1],phL.saved_scal[:,:,2],phL.saved_scal[:,:,3],phL.saved_scal[:,:,4]=compute_mass_flux!(num,grid, grid_u, 
                grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,1,LS[1].geoL)
                end
                # print("\n testvalflux1")
                # for jplot in 1:ny
                #     for iplot in 1:nx
                #         if phL.saved_scal[jplot,iplot,1] !=0.0
                #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e\n" iplot jplot phL.saved_scal[jplot,iplot,1])
                #         end
                #     end
                # end



                # for iscal=1:nb_saved_scalars
                #     @views kill_dead_cells!(phL.saved_scal[:,:,iscal], grid, LS[1].geoL)
                # end




                # printstyled(color=:cyan, @sprintf "\n div(0,grad): %.2e %.2e %.2e %.2e \n" sum(phL.saved_scal[:,:,1]) sum(phL.saved_scal[:,:,2]) sum(phL.saved_scal[:,:,3]) sum(phL.saved_scal[:,:,4]))

                varfluxH2 = phL.saved_scal[:,:,1]
                
                # # #############################################################################################################
                # # #Print values on line
                # id_x = 1
                # for iplot in 1:ny
                #     printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n" iplot grid.y[iplot]/num.plot_xscale phL.saved_scal[iplot,id_x,2] phL.saved_scal[iplot,id_x,3] phL.saved_scal[iplot,id_x,4] phL.trans_scal[iplot,id_x,1] grid.LS[1].u[iplot,id_x] veci(phL.trans_scalD,grid,2)[iplot,id_x,1])
                # end
                # # #############################################################################################################
                # print("\n testvalflux2")

                # for jplot in 1:ny
                #     for iplot in 1:nx
                #         if phL.saved_scal[jplot,iplot,1] !=0.0
                #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e\n" iplot jplot phL.saved_scal[jplot,iplot,1])
                #         end
                #     end
                # end


                # varfluxH2O=compute_mass_flux!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,3)

                # Print values on line
                # for iplot in 1:ny
                #     II = CartesianIndex(iplot, 1) #(id_y, id_x)
                #     pII = lexicographic(II, grid.ny)
                #     printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e \n" iplot grid.y[iplot]/num.plot_xscale mass_flux_vec1[pII] mass_flux_vecb[pII] mass_flux_veci[pII])
                # end
       
                #TODO mass_fluxL
                # mass_fluxL=-... * diffusion_coeff[1] 



                # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 
                varnH2 = -sum(varfluxH2) * diffusion_coeff[1] 

                #TODO mode_2d==0 flux corresponds to cylinder of length 1
                #2D cylinder reference length
                if mode_2d == 1
                    varnH2 .*= num.ref_thickness_2d
                end

                #Pliquid is the average value of p over the bubble interface plus the ambient operating pressure (P).
                p_liq= num.pres0 + mean(veci(phL.pD,grid,2)) #TODO here one bubble
                # p_g=p_liq + 2 * num.σ / current_radius #3D
                p_g=p_liq + num.σ / current_radius #2D

                printstyled(color=:green, @sprintf "\n Mole: %.2e \n" nH2)

                nH2 += varnH2 * num.τ

                printstyled(color=:green, @sprintf "\n Mole: %.2e \n" nH2)

                if varnH2 < 0.0 
                    print(@sprintf "error nH2 %.2e dnH2 %.2e new nH2 %.2e\n" nH2-varnH2*num.τ varnH2*num.τ nH2 )
                    @error ("error nH2")
                    crashed = true
                    nH2 -= varnH2 * num.τ
                    print("wrong nH2 ")
                    # println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                    return nothing
                end
                
                #TODO using temperature0
                if mode_2d == 0
                    current_radius = cbrt(3.0 * nH2 * num.Ru * temperature0/( 4.0 * pi * p_g) )
                elseif mode_2d == 1
                    current_radius = sqrt(nH2 * num.Ru * temperature0/( pi * p_g * num.ref_thickness_2d) )
                elseif mode_2d == 2
                    current_radius = sqrt(nH2/(concentration0[1] * pi))
                elseif mode_2d == 3
                    current_radius = sqrt(2*nH2/(concentration0[1] * pi))
                end

                printstyled(color=:green, @sprintf "\n radius CFL: %.2e \n" (current_radius-previous_radius)/(num.L0/nx))

                if (current_radius-previous_radius)/(num.L0/nx) > 0.5
                    printstyled(color=:red, @sprintf "\n radius CFL: %.2e \n" (current_radius-previous_radius)/(num.L0/nx))
                    @error ("CFL radius")
                    crashed = true
                    return num.current_i
                end

                if nb_saved_scalars>3
                    printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n" nx num.τ num.L0/nx (current_radius-previous_radius)/(num.L0/nx) sum(phL.saved_scal[:,:,1]) sum(phL.saved_scal[:,:,2]) sum(phL.saved_scal[:,:,3]) sum(phL.saved_scal[:,:,4]))
                else
                    printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e %.2e\n" nx num.τ num.L0/nx (current_radius-previous_radius)/(num.L0/nx) sum(phL.saved_scal[:,:,1]))
                end
                printstyled(color=:green, @sprintf "\n n(H2): %.2e added %.2e old R %.2e new R %.2e \n" nH2 varnH2*num.τ previous_radius current_radius)
                printstyled(color=:green, @sprintf "\n p0: %.2e p_liq %.2e p_lapl %.2e \n" num.pres0 p_liq p_g)

                if mode_2d == 3
                    grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (current_radius) * ones(ny, nx)                  
                else
                    grid.LS[1].u .= sqrt.((grid.x .- num.xcoord .- current_radius .+ num.R ).^ 2 + (grid.y .- num.ycoord) .^ 2) - (current_radius) * ones(ny, nx)
                end
                # init_franck!(grid, TL, R, T_inf, 0)
                # u

            elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius")

                #CFL 0.5
                current_radius = current_radius + grid.dx[1,1]/2

                grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (current_radius) * ones(ny, nx)                  


            end #phase change

        end #iLS

        if verbose && adaptative_t
            println("num.τ = $num.τ")
        end

        if advection
            for (iLS, bc) in enumerate(BC_int)
                if is_stefan(bc)
                    IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u)), grid)
                    # u .= sqrt.((x .- num.current_i*Δ/1).^ 2 + y .^ 2) - (0.5) * ones(nx, ny);
                elseif is_fs(bc)
                    rhs_LS .= 0.0
                    LS[iLS].A.nzval .= 0.0
                    LS[iLS].B.nzval .= 0.0
                    IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    utmp .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    rhs_LS .= 0.0
                    S2IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, utmp, LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # Project velocities to the normal and use advecion scheme for advection just
                    # in the normal direction
                    # tmpVx = zeros(grid)
                    # tmpVy = zeros(grid)
                    # V .= 0.0
                    # @inbounds @threads for II in grid.LS[iLS].MIXED
                    #     cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                    #     cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                    #     tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))

                    #     cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                    #     cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                    #     tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                    #     tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                    #     β = atan(tmpVy[II], tmpVx[II])
                    #     if grid.LS[iLS].α[II] > 0.0 && β < 0.0
                    #         β += 2π
                    #     end
                    #     if grid.LS[iLS].α[II] < 0.0 && β > 0.0
                    #         β -= 2π
                    #     end

                    #     V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
                    # end

                    # i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
                    # field_extension!(grid, grid.LS[iLS].u, V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                    # rhs_LS .= 0.0
                    # IIOE_normal!(grid, LS[iLS].A, LS[iLS].B, LS[iLS].u, V, CFL_sc, periodic_x, periodic_y)
                    # BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    # Impose contact angle if a wall is present
                    # rhs_LS .= 0.0
                    # LS[iLS].A.nzval .= 0.0
                    # LS[iLS].B.nzval .= 0.0
                    # for II in grid.ind.all_indices
                    #     pII = lexicographic(II, grid.ny)
                    #     LS[iLS].A[pII,pII] = 1.0
                    #     LS[iLS].B[pII,pII] = 1.0
                    # end
                    # BC_LS_interior!(num, grid, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    # LS[iLS].u .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)
                end
            end
            if analytical
                u[ind.b_top[1]] .= sqrt.(x[ind.b_top[1]] .^ 2 + y[ind.b_top[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[ind.b_bottom[1]] .= sqrt.(x[ind.b_bottom[1]] .^ 2 + y[ind.b_bottom[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[ind.b_left[1]] .= sqrt.(x[ind.b_left[1]] .^ 2 + y[ind.b_left[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[ind.b_right[1]] .= sqrt.(x[ind.b_right[1]] .^ 2 + y[ind.b_right[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
            elseif nb_reinit > 0
                if auto_reinit && (num.current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg, rl_rg_v = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)
                            println("$(ls_rg)")
                            if ls_rg >= δreinit || num.current_i == 1
                                println("yes")
                                RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                                
                                ls_rg, rl_rg_v = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)
                                println("$(ls_rg) ")
                            end
                        end
                    end
                elseif (num.current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                        end
                    end
                # elseif nLS > 1
                #     for iLS in 1:nLS
                #         if !is_wall(BC_int[iLS])
                #             RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, 2nb_reinit, periodic_x, periodic_y, BC_u, BC_int, true)
                #         end
                #     end
                        end
                    end

            # Numerical breakup
            if free_surface && breakup
                count, id_break = breakup_n(LS[1].u, nx, ny, dx, dy, periodic_x, periodic_y, NB_indices, 5e-2)
                println(count)
                if count > count_limit_breakup
                    println("BREAK UP!!") 
                    breakup_f(grid, LS[1].u, id_break)
                    RK2_reinit!(ls_scheme, grid, ind, 1, LS[1].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                end
            end
        end

        if verbose
            if (num.current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) | t = %.2e \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations current_t)
                printstyled(color=:green, @sprintf "\n CFL : %.2e CFL : %.2e num.τ : %.2e\n" CFL max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)*num.τ/Δ num.τ)
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
                        normpS = norm(phS.p.*num.τ)
                        print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                    end
                    if ns_liquid_phase
                        # normuL = norm(phL.u)
                        # normvL = norm(phL.v)
                        # normpL = norm(phL.p.*num.τ)
                        # print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                        if electrolysis
                            print_electrolysis_statistics(nb_transported_scalars,grid,phL) 
                        end 
                    end
                end
            end
        end


        # if levelset && (advection || num.current_i<2 || electrolysis_advection)
        if levelset && (advection || num.current_i<2)
            try
                NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y) 
            catch errorLS
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                printstyled(color=:red, @sprintf "\n LS not updated \n")
                print(errorLS)
                return current_i
            end
            # printstyled(color=:red, @sprintf "\n levelset 4:\n")
            # println(grid.LS[1].geoL.dcap[1,1,:])

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

            if iszero(num.current_i%save_every) || num.current_i==max_iterations
                snap = num.current_i÷save_every+1
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
                # if save_length
                #     fwd.length[snap] = arc_length2(LS[1].geoS.projection, LS[1].MIXED)
                # end
            end
        end

        if navier_stokes
            # if !advection
            #     @time no_slip_condition!(num, grid, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
            #     # grid_u.V .= Δ / (1 * num.τ)
            #     # grid_v.V .= 0.0
            # end

            if ns_solid_phase
                geoS = [LS[iLS].geoS for iLS in 1:_nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:_nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:_nLS]
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                    BC_uS, BC_vS, BC_pS,
                    opC_pS, opC_uS, opC_vS, opS,
                    AuS, BuS, AvS, BvS, AϕS, AuvS, BuvS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, num.current_i, Ra, navier,pres_free_surfaceS,jump_mass_fluxS,mass_fluxS
                )
            end
            if ns_liquid_phase
                geoL = [LS[iLS].geoL for iLS in 1:_nLS]
                geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:_nLS]
                geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:_nLS]
                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    opC_pL, opC_uL, opC_vL, opL,
                    AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, advection, num.current_i, Ra, navier,pres_free_surfaceL,jump_mass_fluxL,mass_fluxL
                )
                # if num.current_i == 1
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

        # cD, cL, D, L = force_coefficients!(num, grid, grid_u, grid_v, opL, fwd, phL; step = num.current_i+1, saveCoeffs = false)

        #update time
        current_t += num.τ

        # if iszero(num.current_i%save_every) || num.current_i==max_iterations
        #     snap = num.current_i÷save_every+1
        #     if num.current_i==max_iterations
        #         snap = size(fwd.T,1)
        #     end
        #     fwd.t[snap] = current_t
        #     @views fwd.V[snap,:,:] .= V
        #     for iLS in 1:_nLS
        #         @views fwd.u[iLS,snap,:,:] .= LS[iLS].u
        #         @views fwd.ux[iLS,snap,:,:] .= grid_u.LS[iLS].u
        #         @views fwd.uy[iLS,snap,:,:] .= grid_v.LS[iLS].u
        #         @views fwd.κ[iLS,snap,:,:] .= LS[iLS].κ
        #     end

        #     if heat_solid_phase && heat_liquid_phase
        #         @views fwd.T[snap,:,:] .= phL.T.*LS[end].geoL.cap[:,:,5] .+ phS.T.*LS[end].geoS.cap[:,:,5]
        #     end
        #     if heat_solid_phase
        #         @views fwd.T[snap,:,:] .= phS.T
        #         @views fwdS.T[snap,:,:] .= phS.T
        #         @views fwdS.TD[snap,:] .= phS.TD
        #     end
        #     if heat_liquid_phase
        #         @views fwd.T[snap,:,:] .= phL.T
        #         @views fwdL.T[snap,:,:] .= phL.T
        #         @views fwdL.TD[snap,:] .= phL.TD
        #     end

        #     #TODO
        #     # if electrolysis_solid_phase && electrolysis_liquid_phase
        #     #     @views fwd.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal].*LS[end].geoL.cap[:,:,5] .+ phS.trans_scal[:,:,iscal].*LS[end].geoS.cap[:,:,5]
        #     # end

        #     if electrolysis_solid_phase #TODO

        #         for iscal=1:nb_transported_scalars
        #             @views fwd.trans_scal[snap,:,:,iscal] .= phS.trans_scal[:,:,iscal]
        #             @views fwdS.trans_scal[snap,:,:,iscal] .= phS.trans_scal[:,:,iscal]
        #             @views fwdS.trans_scalD[snap,:,iscal] .= phS.trans_scalD[:,iscal]
        #         end
        #     end

        #     if electrolysis_liquid_phase #TODO

        #         for iscal=1:nb_transported_scalars
        #             @views fwd.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal].*LS[end].geoL.cap[:,:,5] .+ phS.trans_scal[:,:,iscal].*LS[end].geoS.cap[:,:,5]

        #             # @views fwd.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal]
        #             @views fwdL.trans_scal[snap,:,:,iscal] .= phL.trans_scal[:,:,iscal]
        #             @views fwdL.trans_scalD[snap,:,iscal] .= phL.trans_scalD[:,iscal]
                    
        #             # printstyled(color=:cyan, @sprintf "\n write scal \n")
        #             # print(minimum(fwdL.trans_scal[snap,:,:,iscal]), "phL ", minimum(phL.trans_scal[:,:,iscal]))
        #         end
                
        #         # @views fwd.mass_flux[snap,:,:] .= phL.mass_flux #reshape(varfluxH2, grid)

        #         for iscal=1:nb_saved_scalars
        #             @views fwd.saved_scal[snap,:,:,iscal] .= phL.saved_scal[:,:,iscal] #varflux[iscal] #reshape(varfluxH2, grid)
        #         end

        #         @views fwdL.phi_ele[snap,:,:] .= phL.phi_ele

        #         @views fwdL.phi_eleD[snap,:] .= phL.phi_eleD

        #         @views fwdL.i_current_mag[snap,:,:] .= phL.i_current_mag
        #         @views fwdS.Eu[snap,:,:] .= phS.Eu
        #         @views fwdS.Ev[snap,:,:] .= phS.Ev
        #         @views fwdL.Eu[snap,:,:] .= phL.Eu
        #         @views fwdL.Ev[snap,:,:] .= phL.Ev


        #     end
        #     if ns_solid_phase
        #         @views fwdS.p[snap,:,:] .= phS.p
        #         @views fwdS.pD[snap,:] .= phS.pD
        #         @views fwdS.ϕ[snap,:,:] .= phS.ϕ
        #         @views fwdS.u[snap,:,:] .= phS.u
        #         @views fwdS.v[snap,:,:] .= phS.v
        #         @views fwdS.ucorrD[snap,:,:] .= phS.ucorrD
        #         @views fwdS.vcorrD[snap,:,:] .= phS.vcorrD
        #     end
        #     if ns_liquid_phase
        #         @views fwdL.p[snap,:,:] .= phL.p
        #         @views fwdL.pD[snap,:]  .= phL.pD
        #         @views fwdL.ϕ[snap,:,:] .= phL.ϕ
        #         @views fwdL.u[snap,:,:] .= phL.u
        #         @views fwdL.v[snap,:,:] .= phL.v
        #         @views fwdL.vD[snap,:]  .= phL.vD

        #         @views fwdL.ucorrD[snap,:,:] .= phL.ucorrD
        #         @views fwdL.vcorrD[snap,:,:] .= phL.vcorrD
        #         # @views fwd.Cd[snap] = cD
        #         # @views fwd.Cl[snap] = cL
        #         @views fwd.radius[snap] = current_radius
        #     end
        #     if advection
        #         fwdS.Vratio[snap] = volume(LS[end].geoS) / V0S
        #         fwdL.Vratio[snap] = volume(LS[end].geoL) / V0L
        #     end
        # end
        # @views fwd.Cd[num.current_i+1] = cD
        # @views fwd.Cl[num.current_i+1] = cL
        # # @views fwd.radius[num.current_i+1] = current_radius




        if electrolysis

            ####################################################################################################
            #PDI (IO)
            ####################################################################################################

            if num.io_pdi>0

                try
                    printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                    time = current_t #Cdouble
                    nstep = num.current_i
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
            
                    Eus,Evs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.Eu, phL.Ev)
            
                    us,vs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.u,phL.v)
            
                    # print("\n before write \n ")
            
                    iLSpdi = 1 # all LS iLS = 1 # or all LS ?


                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name

                    @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                    "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                    C_NULL::Ptr{Cvoid})::Cvoid
            
                    # print("\n after write \n ")
            
                    # @ccall "libpdi".PDI_finalize()::Cvoid
            
                    # printstyled(color=:red, @sprintf "\n PDI test end\n" )
            
                catch error
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                    print(error)
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                end

                # try
                #     printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                #     time = current_t #Cdouble
                #     nstep = num.current_i
                #     print("\n nstep ",typeof(nstep))
                #     # pdi_array =zeros(nx,ny)
            
                #     # for j in 1:grid.ny
                #     #     for i in 1:grid.nx
                #     #         pdi_array[j,i]=1000*i+j
                #     #     end
                #     # end
            
                #     print("\n before write \n ")
            
                #     # @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                #     #             "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Clonglong,
                #     #             "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Clonglong,
                #     #             "main_field"::Cstring, pdi_array::Ptr{Cdouble}, 
                #     #             PDI_OUT::Clonglong,
                #     #             C_NULL::Ptr{Cvoid})::Cvoid

                #     @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                #     "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Clonglong,
                #     "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Clonglong,
                #     "u"::Cstring, phL.u::Ptr{Cdouble}, 
                #     "v"::Cstring, phL.v::Ptr{Cdouble}, 
                #     "trans_scal"::Cstring, phL.trans_scal::Ptr{Cdouble}, 
                #     "phi_ele"::Cstring, phL.v::Ptr{Cdouble}, 
                #     PDI_OUT::Clonglong,
                #     C_NULL::Ptr{Cvoid})::Cvoid
            
                #     print("\n after write \n ")
            
                #     @ccall "libpdi".PDI_finalize()::Cvoid
            
                #     printstyled(color=:red, @sprintf "\n PDI test end\n" )
            
                # catch error
                #     printstyled(color=:red, @sprintf "\n PDI error \n")
                #     print(error)
                #     printstyled(color=:red, @sprintf "\n PDI error \n")
                # end

            end #if io_pdi

            ####################################################################################################


            if crashed #due to nH2<0...
                return num.current_i
            end

        
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8 || norm(phL.trans_scal) > 1e8 || norm(phL.phi_ele) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                
                print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
    

                print_electrolysis_statistics(nb_transported_scalars,grid,phL)

                crashed=true
                # return nothing
                return num.current_i

            end
        else
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
               
                crashed=true
                return nothing
                
            end

        end

        num.current_i += 1

        # if adaptative_t
           
        # elseif adaptative_t
        #     num.τ = min(CFL*Δ^2*Re, CFL*Δ/max(
        #         abs.(V)..., abs.(grid_u.V)..., abs.(grid_v.V)..., 
        #         abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)
        #     )
        # end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) | t = %.2e \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations current_t)
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
                    normpS = norm(phS.p.*num.τ)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*num.τ)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    if electrolysis
                       print_electrolysis_statistics(nb_transported_scalars,grid,phL)
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
    # elseif flapping
    #     return xc, yc
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

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, num.τ, L0, NB, max_iterations, num.current_i, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
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

    num.current_i = max_iterations + 1

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

    laplacian!(dir, num, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
    laplacian!(dir, num, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])

    while num.current_i > 1

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

                    laplacian!(dir, num, opS.LT, opS.CUTT, bcSx, bcSy, geoS.dcap, ny, BC_TS, inside, LIQUID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoS, opS)
                    TS .= reshape(gmres(opS.A,(opS.B*vec(TS) + 2.0*num.τ*opS.CUTT)), (ny,nx))
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

                    laplacian!(dir, num, opL.LT, opL.CUTT, bcLx, bcLy, geoL.dcap, ny, BC_TL, inside, SOLID,
                                MIXED, b_left[1], b_bottom[1], b_right[1], b_top[1])
                    crank_nicolson!(num, grid, geoL, opL)
                    TL .= reshape(gmres(opL.A,(opL.B*vec(TL) + 2.0*num.τ*opL.CUTT)), (ny,nx))
                end
            catch
                @error ("Unphysical temperature field, iteration $num.current_i")
                break
            end
        end

        if verbose
            if num.current_i%show_every == 0
                try
                    printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations)
                    printstyled(color=:green, @sprintf "\n CFL : %.2e CFL : %.2e num.τ : %.2e\n" CFL max(abs.(V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)*num.τ/Δ num.τ)

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

        num.current_i -= 1
        κ .= κsave[num.current_i,:,:]
        u .= usave[num.current_i,:,:]
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) \n" (num.current_i-1) 100*(num.current_i-1)/max_iterations)
            print(@sprintf "V_mean = %.2f  V_max = %.2f  V_min = %.2f  V_stdev = %.5f\n" mean(V[MIXED]) findmax(V[MIXED])[1] findmin(V[MIXED])[1] std(V[MIXED]))
            print(@sprintf "κ_mean = %.2f  κ_max = %.2f  κ_min = %.2f  κ_stdev = %.5f\n" mean(κ[MIXED]) findmax(κ[MIXED])[1] findmin(κ[MIXED])[1] std(κ[MIXED]))
            print("\n \n")
        catch
            @show (length(MIXED))
        end
    end
end
