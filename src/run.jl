"""

Main function of Flower.jl code to run a simulation

"""
function run_forward!(
    num::Numerical{Float64, Int64},
    grid_p::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    op::DiscreteOperators{Float64, Int64}, 
    phS::Union{Phase{Float64},Nothing}, 
    phL::Phase{Float64};
    periodic_x::Bool = false,
    periodic_y::Bool = false,
    BC_TS = Boundaries(),
    BC_TL = Boundaries(),
    BC_pS = Boundaries(),
    BC_pL = Boundaries(),
    BC_uS = BoundariesInt(),
    BC_uL = BoundariesInt(),
    BC_vS = BoundariesInt(),
    BC_vL = BoundariesInt(),
    BC_u = Boundaries(),
    BC_trans_scal = Vector{BoundariesInt}(),
    BC_phi_ele = BoundariesInt(),
    BC_int::Vector{<:BoundaryCondition} = [WallNoSlip()],
    time_scheme::TemporalIntegration = CN,
    ls_scheme::LevelsetDiscretization = weno5,
    auto_reinit::Int64 = 0,
    heat::Bool = false,
    heat_convection::Bool = false,
    heat_liquid_phase::Bool = false,
    heat_solid_phase::Bool = false,
    navier_stokes ::Bool= false,
    ns_advection::Bool = false,
    ns_liquid_phase::Bool = false,
    ns_solid_phase::Bool = false,
    hill::Bool = false,
    Vmean::Bool = false,
    levelset::Bool = true,
    analytical::Bool = false,
    verbose::Bool = false,
    show_every::Int64 = 100,
    save_radius::Bool = false,
    adaptative_t::Bool = false,
    breakup::Int64 = 0,
    Ra::Float64 = 0.0,
    electrolysis::Bool = false,
    electrolysis_convection::Bool = false,  
    electrolysis_liquid_phase::Bool = false,
    electrolysis_solid_phase::Bool = false,
    electrolysis_phase_change_case::String = "Khalighi",
    imposed_velocity::String = "none",
    adapt_timestep_mode::Int64 = 0,
    non_dimensionalize::Int64=1,
    mode_2d::Int64=0,
    test_laplacian::Bool = false,
    )

    #region Initialize simulation parameters

    λ = 1 #for Stefan velocity
    speed = 0.0

    radius_pdi = [0.0] #radius for pdi, list so that value is mutable

    if num.time > num.nucleation_time && num.phase_change_method >0 #TODO more precisely no mass transfer but velocity     
        num.phase_change_currently_activated = 1
    end

    if num.mu1 == 0.0 && num.mu2 == 0.0 && num.mu_one_fluid_average !=0
        @error ("Resetting num.mu_one_fluid_average to 0")
        num.mu_one_fluid_average = 0
    end

    # if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
    #     if CL BC 
    #     # grid_u.LS[1].cl
    #     @error("BC CL with one fluid")
    #     end
    # end
    
    #TODO Re
    #TODO homogenize set_poisson... functions which use -b - in the original implementation and not +b + for the Robin BC

    # Defining epsilon length and volume to handle special cells
    if num.epsilon_mode == 1 || num.epsilon_mode ==2
        num.epsilon_dist = eps(0.01) * num.Δ
        num.epsilon_vol = (eps(0.01)*num.Δ)^2
        num.epsilon_dist_mass_transfer_rate =  num.Δ / 10 #TODO improve ex crit vol cf num.epsilon_volume_fraction_phase_change
        #TODO kill dead cells
        #TODO 1e-...
    end

    # if monophasic
    # if length(BC_int) != num.nLS
    #     @error ("You have to specify $(num.nLS) boundary conditions.")
    #     return nothing
    # end

    num.status =  0 # used to stop the simulation
    num.current_iter = 0 #1
    num.time = 0.
    free_surface = false
    stefan = false
    navier = false
    ls_advection = true
    if any(is_fs, BC_int)
        free_surface = true
    end
    if any(is_stefan, BC_int)
        stefan = true
    end
    if any(is_navier_cl, BC_int) || any(is_navier, BC_int)
        navier = true
    end

    if num.nNavier > 1
        @warn ("When using more than 1 Navier BC, the interfaces shouldn't cross")
    end

   

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan || (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) && (num.activate_interface != 0) #|| electrolysis_phase_change_case !="none"
        # num.activate_interface != 0 : no interface, do not advect (monophasic)
        advection = true
    else
        advection = false
    end

    if num.solve_Navier_Stokes_liquid_phase == 0
        advection = false
        electrolysis_advection = false

    else
        if electrolysis
            electrolysis_advection = true
        end
    end

    if num.advection_LS_mode == 16 || num.advection_LS_mode == 16 #test Cipriano 2024 's method
        extend_liquid_velocity = true
    else
        extend_liquid_velocity = false
    end

    # The count threshold shouldn't be smaller than 2
    count_limit_breakup = 6
    
    if num.verbosity > 0
        printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e\n" num.CFL num.timestep_n)
    end

    #compute initial curvature (TODO other interfaces than circle)
    num.mean_curvature = 1.0/num.R 
    num.current_radius = num.R

    # if num.surface_tension == 0
    #     compute_surface_tension_VOF!(num,grid_p, grid_u, grid_v, opC_p, opC_u, opC_v, 
    #     volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)
    # elseif num.surface_tension == 1
    #     compute_surface_tension_LS!(num,grid_p, grid_u, grid_v, opC_p, opC_u, opC_v, 
    #     volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0,
    #     levelset_1D, levelset_heavyside_2D, normal_and_dirac_u, normal_and_dirac_v,
    #     normal_u, normal_v, curvature_u, curvature_v)
    # end

    
    num.timestep_n = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)

    pres_free_surfaceS = 0.0
    # pres_free_surfaceL = 0.0
    pres_free_surfaceL = num.pres0

    if occursin("levelset",electrolysis_phase_change_case)
        jump_mass_transfer_rateS = false 
        jump_mass_transfer_rateL = true
    else
        jump_mass_transfer_rateS = false 
        jump_mass_transfer_rateL = false
    end

    mass_transfer_rateS = 0.0
    
    iRe = 1.0 / num.Re

    # CFL_sc is not a CFL, it is supposed to be the timestep divided by the cell volume
    CFL_sc = num.timestep_n / num.Δ^2

    irho = 1.0

    num.mu_cin1 = num.mu1 / num.rho1
    num.mu_cin2 = num.mu2 / num.rho2

    if non_dimensionalize==0
        #force L=1 u=1
        Re=num.rho1/num.mu1
        iRe = 1.0/Re
        irho=1.0/num.rho1
        num.visc_coeff=iRe
    else 
        num.visc_coeff=iRe
        Re = num.Re
    end

    if num.verbosity > 0
        printstyled(color=:green, @sprintf "\n Re : %.2e %.2e\n" Re num.visc_coeff)
        printstyled(color=:magenta, @sprintf "\n CFL_sc : %.2e\n" CFL_sc)
    end

    # Electrolysis
    num.current_radius = 0.0

    # TODO kill_dead_cells! for [:,:,iscal]
    #region init number of moles
    if electrolysis
        if electrolysis_phase_change_case != "none"
            num.current_radius = num.R

            printstyled(color=:green, @sprintf "\n radius: %.2e \n" num.current_radius)

            p_liq= num.pres0 + mean(veci(phL.pD,grid_p,2)) #TODO here one bubble
            p_g=p_liq + 2 * num.σ / num.current_radius

            #TODO using num.temperature0
            if mode_2d==0
                num.nH2 = p_g * 4.0 / 3.0 * pi * num.current_radius ^ 3 / (num.temperature0 * num.Ru) 
            elseif mode_2d == 1 #reference thickness for a cylinder
                num.nH2 = p_g * pi * num.current_radius ^ 2 * num.ref_thickness_2d / (num.temperature0 * num.Ru) 
            elseif mode_2d==2 #mol/meter
                num.nH2=num.concentration0[num.index_phase_change]* pi * num.current_radius ^ 2
            elseif mode_2d==3 #mol/meter half circle
                num.nH2=1.0/2.0*num.concentration0[num.index_phase_change]* pi * num.current_radius ^ 2
            end
            # num.nH2 = 4.0/3.0 * pi * num.current_radius^3 * num.rho2 / num.MWH2

            printstyled(color=:green, @sprintf "\n Mole: %.2e \n" num.nH2)

            printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" num.concentration0[num.index_phase_change]*4.0/3.0*pi*num.current_radius^3 p_g*4.0/3.0*pi*num.current_radius^3/(num.temperature0*num.Ru))

        end
    end # if electrolysis    
    #endregion init number of moles

    #endregion Initialize simulation parameters

    local NB_indices;

    #region Allocations

    if num.solve_solid == 1
        local Cum1S = fzeros(grid_u)
        local Cvm1S = fzeros(grid_v)
        local Mm1_S
        local Mum1_S
        local Mvm1_S
    end


    local Cum1L = fzeros(grid_u)
    local Cvm1L = fzeros(grid_v)
    local Mm1_L
    local Mum1_L
    local Mvm1_L
   
    if extend_liquid_velocity #num.advection_LS_mode == 16
        # *_ext_vel : variables for the second system of Navier-Stokes equations
        local Cum1L_ext_vel = fzeros(grid_u)
        local Cvm1L_ext_vel = fzeros(grid_v)
        # local Mm1_L_ext_vel
        # local Mum1_L_ext_vel
        # local Mvm1_L_ext_vel
        p_ext_vel = zeros(grid_p)
        pD_ext_vel = fzeros(grid_p)
        phi_ext_vel = zeros(grid_p)
        u_ext_vel = zeros(grid_u)
        v_ext_vel= zeros(grid_v)
        u_predictionD_ext_vel= fzeros(grid_u)
        v_predictionD_ext_vel= fzeros(grid_v)
        uD_ext_vel= fzeros(grid_u)
        vD_ext_vel= fzeros(grid_v)
        u_prediction_ext_vel= zeros(grid_u)
        v_prediction_ext_vel= zeros(grid_v)
        uT_ext_vel = zeros(grid_p)
        pres_grad_x = fzeros(grid_u)
        pres_grad_y = fzeros(grid_v)
    else 
        u_ext_vel = nothing
        v_ext_vel= nothing

    end

    θ_out = zeros(grid_p, 4)
    utmp = copy(grid_p.LS[1].u)
    rhs_LS = fzeros(grid_p)

    #vectors used reset at start of functions to limit allocations
    tmp_vec_u = zeros(grid_u) 
    tmp_vec_v = zeros(grid_v) 
    tmp_vec_u0 = zeros(grid_u) 
    tmp_vec_v0 = zeros(grid_v)
   
    tmp_vec_u1 = zeros(grid_u) 
    tmp_vec_v1 = zeros(grid_v)

    tmp_vec_p = zeros(grid_p) 
    tmp_vec_p0 = zeros(grid_p) 
    tmp_vec_p1 = zeros(grid_p) 

    tmp_vec_1D = fnzeros(grid_p,num)
    tmp_vec_1D_2 = fnzeros(grid_p,num)

    # tmp_vec_1D_u = fnzeros(grid_p,num)
    tmp_vec_1D_v = fnzeros(grid_v,num)
    tmp_vec_1D_v0 = fnzeros(grid_v,num)

    tmp_vec_1D_p = fnzeros(grid_p,num)
    tmp_vec_1D_p0 = fnzeros(grid_p,num)


    # nb_gaz_acceptors = zeros(grid_p)

    #region electrolysis
    #TODO harmonic conductivity
    if electrolysis
        if num.nb_transported_scalars>1
            elec_cond = zeros(grid_p)
            elec_condD = fnzeros(grid_p,num)

            if heat 
                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
                elec_cond .= reshape(vec1(elec_condD,grid_p),grid_p)
                # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
                # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
            else
                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
                elec_cond .= reshape(vec1(elec_condD,grid_p),grid_p)
                # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
                # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
            end
        else 
            elec_cond = ones(grid_p)
            elec_condD = fnones(grid_p,num)
            printstyled(color=:green, @sprintf "\n conductivity one")
        end 

        if num.electrolysis_reaction == "Butler_no_concentration"
            i_butler = zeros(grid_p.ny) #left wall
        elseif num.electrolysis_reaction == "fixed_current"
            i_butler = zeros(grid_p.nx) #bottom wall

        end

    end #electrolysis
    #endregion electrolysis

    if levelset

        # At every iteration, update_all_ls_data is called twice, 
        # once inside run.jl and another one (if there's advection of the levelset) inside set_heat!. 
        # The difference between both is a flag as last argument, inside run.jl is 
        # implicitly defined as true and inside set_heat is false. 
        # If you're calling your version of set_heat several times, then you're calling the version with the flag set to false, but for the convective term it has to be set to true, so that's why
        # The flag=true, the capacities are set for the convection, the flag=false they are set for the other operators

        NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y)

        # printstyled(color=:red, @sprintf "\n levelset:\n")
        # println(grid_p.LS[1].geoL.dcap[1,1,:])

        if save_radius
            n_snaps = iszero(num.max_iterations%num.save_every) ? num.max_iterations÷num.save_every+1 : num.max_iterations÷num.save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid_p, grid_p.LS[1])
        end
        if hill
            local radius = zeros(num.max_iterations+1)
            a = zeros(length(grid_p.LS[1].MIXED))
            for i in eachindex(grid_p.LS[1].MIXED)
                a[i] = grid_p.LS[1].geoL.projection[grid_p.LS[1].MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        grid_p.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_u.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_v.LS[1].MIXED = [CartesianIndex(-1,-1)]
    end

    # if save_length
    #     fwd.length[1] = arc_length2(grid_p.LS[1].geoS.projection, grid_p.LS[1].MIXED)
    # end

    #endregion


    #region Initialisation of bulk and interfacial values
    # Initialisation of bulk and interfacial values
    
    #TODO which grid_p.LS
 
    #TODO perio, intfc, ... check init_fields_2!

    #No scalar in "solid" phase
    # @views init_fields_multiple_levelsets!(num,phS.trans_scalD[:,iscal],phS.trans_scal[:,:,iscal],HS,BC_trans_scal[iscal],grid_p,num.concentration0[iscal])
    # @views phS.trans_scal[:,:,iscal] .= num.concentration0[iscal]

    # TODO reset zero
    tmp_vec_u .= 0.0

    if num.solve_solid == 1
        #from LS 1 to centroid grid_u.LS[end].geoS
        get_height!(grid_u.LS[1],grid_u.ind,grid_u.dx,grid_u.dy,grid_u.LS[end].geoS,tmp_vec_u) #here tmp_vec_u solid

        init_fields_multiple_levelsets!(num,phS.uD,phS.u,tmp_vec_u,BC_uS,grid_u,num.uD,"uS")
        # init_fields_multiple_levelsets!(num,phS.u_predictionD,phS.u,HSu,BC_uS,grid_u,num.uD)
    end

    get_height!(grid_u.LS[1],grid_u.ind,grid_u.dx,grid_u.dy,grid_u.LS[end].geoL,tmp_vec_u)  #here tmp_vec_u liquid

    init_fields_multiple_levelsets!(num,phL.uD,phL.u,tmp_vec_u,BC_uL,grid_u,num.uD,"uL")
    # init_fields_multiple_levelsets!(num,phL.u_predictionD,phL.u,HLu,BC_uL,grid_u,num.uD)


    # TODO reset zero
    tmp_vec_v .= 0.0
    if num.solve_solid == 1
        get_height!(grid_v.LS[1],grid_v.ind,grid_v.dx,grid_v.dy,grid_v.LS[end].geoS,tmp_vec_v) 

        init_fields_multiple_levelsets!(num,phS.vD,phS.v,tmp_vec_v,BC_vS,grid_v,num.vD,"vS")
        # init_fields_multiple_levelsets!(num,phS.v_predictionD,phS.v,HSv,BC_vS,grid_v,num.vD)
    end

    get_height!(grid_v.LS[1],grid_v.ind,grid_v.dx,grid_v.dy,grid_v.LS[end].geoL,tmp_vec_v)

    init_fields_multiple_levelsets!(num,phL.vD,phL.v,tmp_vec_v,BC_vL,grid_v,num.vD,"vL")
    # init_fields_multiple_levelsets!(num,phL.v_predictionD,phL.v,HLv,BC_vL,grid_v,num.vD)


    # TODO reset zero
    tmp_vec_p .= 0.0

    if num.solve_solid == 1

        get_height!(grid_p.LS[1],grid_p.ind,grid_p.dx,grid_p.dy,grid_p.LS[end].geoS,tmp_vec_p) #here tmp_vec_p solid

        init_fields_multiple_levelsets!(num,phS.pD,phS.p,tmp_vec_p,BC_pS,grid_p,num.pres_intfc,"pS")

        if heat
            init_fields_multiple_levelsets!(num,phS.TD,phS.T,tmp_vec_p,BC_TS,grid_p,num.θd,"TS")
        end

        #Electrolysis
        if electrolysis && num.electrical_potential > 0
            init_fields_multiple_levelsets!(num,phS.phi_eleD,phS.phi_ele,tmp_vec_p,BC_phi_ele,grid_p,num.phi_ele0,"phiS")
        end
    end

    get_height!(grid_p.LS[1],grid_p.ind,grid_p.dx,grid_p.dy,grid_p.LS[end].geoL,tmp_vec_p) #here tmp_vec_p liquid

    init_fields_multiple_levelsets!(num,phL.pD,phL.p,tmp_vec_p,BC_pL,grid_p,num.pres_intfc,"pL")

    if heat
        init_fields_multiple_levelsets!(num,phL.TD,phL.T,tmp_vec_p,BC_TL,grid_p,num.θd,"TL")
    end

    # Electrolysis
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check %s %s %s %s %.2e %.2e %2i\n" heat heat_convection electrolysis electrolysis_convection num.timestep_n num.θd num.nb_transported_scalars)

        for iscal=1:num.nb_transported_scalars
            @views phL.trans_scal[:,:,iscal] .= num.concentration0[iscal]
            @views init_fields_multiple_levelsets!(num,phL.trans_scalD[:,iscal],phL.trans_scal[:,:,iscal],
            tmp_vec_p,BC_trans_scal[iscal],grid_p,num.concentration0[iscal],"scalL")
            # tmp_vec_p,BC_trans_scal[iscal],grid_p,num.concentration0[iscal],"testscalL")

        end

        if num.electrical_potential > 0 #TODO init phi =0 or with Neumann for BC of concentration? for now not done since "phiL" given in arg
            init_fields_multiple_levelsets!(num,phL.phi_eleD,phL.phi_ele,tmp_vec_p,BC_phi_ele,grid_p,num.phi_ele0,"phiL")
            if num.electrical_potential == 3
                vecb(phL.phi_eleD,grid_p) .= 0.0
            end
        end
    end  
    
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check init_fields_2!\n")
        # print_electrolysis_statistics(num,grid_p,phL)

        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
        # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        if num.solve_solid == 1

            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                norm(phL.u) > 1e8 || norm(phL.T) > 1e8 || 
                norm(phS.u) > 1e8 || norm(phS.T) > 1e8 
                # norm(phL.trans_scal) > 1e8 
                || norm(phL.phi_ele) > 1e8 
                # ||
                # any(phL.trans_scal .<0)
                )
                println(@sprintf "\n CRASHED start \n")
                
                print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) ,
                "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , 
                "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 ,
                "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8,"\n any(phL.trans_scal .<0): ", any(phL.trans_scal .<0))

                num.status =  1
                return num.current_iter

            end
        else 

            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                norm(phL.u) > 1e8 || norm(phL.T) > 1e8 || 
                # norm(phL.trans_scal) > 1e8 || 
                norm(phL.phi_ele) > 1e8 
                # || any(phL.trans_scal .<0)
                )
                println(@sprintf "\n CRASHED start \n")

                # println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                
                print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) ,
                "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                "\n phL.u: ",norm(phL.u) > 1e8, "\n phL.T: ",norm(phL.T) > 1e8 , 
                "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 ,
                "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8,"\n any(phL.trans_scal .<0): ", any(phL.trans_scal .<0))

                num.status =  1
                return num.current_iter

            end

        end

        interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)

        PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
            "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
            "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
            "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
            "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,
            "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
            "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
            # "normal_x"::Cstring, normal_x::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "normal_y"::Cstring, normal_y::Ptr{Cdouble}, PDI_OUT::Cint,  
            # grid_u.LS[iLS].α
            # "normal_angle"::Cstring, grid_p.LS[num.iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
            "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
            "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
            "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
            # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
            # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
            # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint




    end

    if num.solve_solid == 1
        kill_dead_cells!(phS.T, grid_p, grid_p.LS[end].geoS)
    end
    kill_dead_cells!(phL.T, grid_p, grid_p.LS[end].geoL)


    #TODO check timestep coefficients num.n-1 

    #Electrolysis
    # TODO kill_dead_cells! for [:,:,iscal]
    if electrolysis
        for iscal=1:num.nb_transported_scalars
            # @views kill_dead_cells!(phS.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoS) #TODO
            # @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL) 
            # @views kill_dead_cells_val!(phS.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoS) #TODO
            # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,num.concentration0[iscal]) 
            # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,0.0) 

            if num.kill_dead_cells == 0
                @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,0.0)
            else
                @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,num.concentration0[iscal])
            end

            @views veci(phL.trans_scalD[:,iscal],grid_p,1) .= vec(phL.trans_scal[:,:,iscal])


        end
    end #if electrolysis    

    #endregion
    





    # print("\n vecb_L(elec_condD, grid_p) after kill \n ", vecb_L(phL.trans_scalD[:,2], grid_p) )


    # if num.nb_transported_scalars>0
    #     printstyled(color=:green, @sprintf "\n after kill \n")
    #     print_electrolysis_statistics(num,phL)
    #     printstyled(color=:green, @sprintf "\n average T %s\n" average!(phL.T, grid_p, grid_p.LS[1].geoL, num))
    # end

    #region Allocations

    # Set matrices/operators
    if is_Forward_Euler(time_scheme) || is_Crank_Nicolson(time_scheme)
        NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        # printstyled(color=:red, @sprintf "\n levelset 2:\n")
        # println(grid_p.LS[1].geoL.dcap[1,1,:])

        if navier_stokes || heat || electrolysis
            if num.solve_solid == 1
                geoS = [grid_p.LS[iLS].geoS for iLS in 1:num._nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:num._nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:num._nLS]
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                    num, grid_p, geoS, grid_u, geo_uS, grid_v, geo_vS,
                    op.opC_pS, op.opC_uS, op.opC_vS, periodic_x, periodic_y
                )
            end

            geoL = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                num, grid_p, geoL, grid_u, geo_uL, grid_v, geo_vL,
                op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y
            )
        end

        # Contains volumes
        if num.solve_solid == 1
            Mm1_S = copy(op.opC_pS.M)
            Mum1_S = copy(op.opC_uS.M)
            Mvm1_S = copy(op.opC_vS.M)
        end

        # M defined in set_cutcell_matrices!
        Mm1_L = copy(op.opC_pL.M)
        Mum1_L = copy(op.opC_uL.M)
        Mvm1_L = copy(op.opC_vL.M)

        if num.one_fluid_model == 1
            @error("\n Mum1L needs to be redefined for one-fluid ")
       
        end


        if navier_stokes || heat || electrolysis

            #Allocations
            #TODO pre-allocate at start to save up allocations
            #TODO optimize allocations

            # #Preallocate for mass flux computations
            # if electrolysis_phase_change_case != "None"
            #     mass_transfer_rate_vec1 = fzeros(grid_p)
            #     mass_transfer_rate_vecb = fzeros(grid_p)
            #     mass_transfer_rate_veci = fzeros(grid_p)
            #     mass_transfer_rate = zeros(grid_p)
            # else
            #     mass_transfer_rate = 0.0
            # end
            mass_transfer_rate_vec1 = fzeros(grid_p)
            mass_transfer_rate_vecb = fzeros(grid_p)
            mass_transfer_rate_veci = fzeros(grid_p)
            mass_transfer_rate = zeros(grid_p)

            mass_transfer_rate_redistributed = zeros(grid_p)
            nb_gaz_acceptors = zeros(Int64,grid_p.ny,grid_p.nx)
            # nb_gaz_acceptors = zeros(grid_p)

            #Allocations for scalar grid_p
            ni = grid_p.nx * grid_p.ny
            nb = 2 * grid_p.nx + 2 * grid_p.ny
            nt = (num.nLS + 1) * ni + nb

            a1_p = spdiagm(ni,ni,zeros(ni))

            Ascal = spzeros(nt, nt)
            Bscal = spzeros(nt, nt)
            rhs_scal = fnzeros(grid_p, num)
            F_residual = fnzeros(grid_p, num)


            all_CUTCT = zeros(grid_p.ny * grid_p.nx, num.nb_transported_scalars)

            if num.one_fluid_model == 0
                if num.solve_solid == 1
                    AϕS = spzeros(nt, nt)
                end
                AϕL = spzeros(nt, nt)

            else
                nt_pressure = ni + nb
                AϕL = spzeros(nt_pressure, nt_pressure)
                rhs_phi = zeros(nt_pressure)
            end

            if electrolysis 
                # Aphi_eleL = spzeros(nt, nt)

                # coeffDu = zeros(grid_u)
                # coeffDv = zeros(grid_v)
                coeffDx_interface = zeros(grid_u)
                coeffDy_interface = zeros(grid_v)
            end

            if num.solve_solid == 1
                ATS = spzeros(nt, nt)
                BTS = spzeros(nt, nt)
            end
            ATL = spzeros(nt, nt)
            BTL = spzeros(nt, nt)

            ni = grid_u.nx * grid_u.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny
            nt = (num.nLS + 1) * ni + nb
            if num.solve_solid == 1
            AuS = spzeros(nt, nt)
            BuS = spzeros(nt, nt)
            end
            AuL = spzeros(nt, nt)
            BuL = spzeros(nt, nt)

            ni = grid_v.nx * grid_v.ny
            nb = 2 * grid_v.nx + 2 * grid_v.ny
            nt = (num.nLS + 1) * ni + nb
            if ns_solid_phase
                AvS = spzeros(nt, nt)
                BvS = spzeros(nt, nt)
            end
            AvL = spzeros(nt, nt)
            BvL = spzeros(nt, nt)

            #region Allocate for Navier case
            ni = grid_u.nx * grid_u.ny + grid_v.nx * grid_v.ny
            nb = 2 * grid_u.nx + 2 * grid_u.ny + 2 * grid_v.nx + 2 * grid_v.ny
            nt = (num.nLS - num.nNavier + 1) * ni + num.nNavier * grid_p.nx * grid_p.ny + nb
            #dev was done with nNavier == ?
            
            #when no Navier: nt = (num.nLS + 1) * ni + nb


            # if ((num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) && num.pressure_velocity_coupling != 0)
            #     @error("\nCoupled pressure velocity + one-fluid model error")
            #     return
            # end


            if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
                rho_one_fluid = zeros(grid_p)
                mu_one_fluid  = zeros(grid_p)
                volume_fraction = zeros(grid_p)
                interface_length = zeros(grid_p)
                levelset_one_fluid = zeros(grid_p)

                rho_one_fluid_u = zeros(grid_u)
                # mu_one_fluid_u  = zeros(grid_u)

                rho_one_fluid_v = zeros(grid_v)
                # mu_one_fluid_v  = zeros(grid_v)

                volumic_surface_tension_u = zeros(grid_u)
                volumic_surface_tension_v = zeros(grid_v)

                # Pre-allocate arrays
                levelset_1D = fnzeros(grid_p, num)
                levelset_heavyside_2D = zeros(grid_p)

                convection_u = fzeros(grid_u)
                convection_v = fzeros(grid_v)

                viscosity_coeff_for_du_dx = zeros(grid_u.ny, grid_u.nx+1)
                viscosity_coeff_for_du_dy = zeros(grid_u.ny+1, grid_u.nx)
                viscosity_coeff_for_dv_dx = zeros(grid_v.ny, grid_v.nx+1)
                viscosity_coeff_for_dv_dy = zeros(grid_v.ny+1, grid_v.nx)

                velocity_and_BC_convection_u_x = zeros(grid_u) #Du_x
                velocity_and_BC_convection_u_y = zeros(grid_u)
                velocity_and_BC_convection_v_x = zeros(grid_v)
                velocity_and_BC_convection_v_y = zeros(grid_v)

            end

            if (num.pressure_velocity_coupling == 0 && num.one_fluid_model == 0)
            
                if ns_solid_phase
                    AuvS = spzeros(nt, nt)
                    BuvS = spzeros(nt, nt)
                end

                AuvL = spzeros(nt, nt)
                BuvL = spzeros(nt, nt)

            else
                # We solve for u, v with borders 
                # dev was done with nNavier == ?

                ni_p = grid_p.nx * grid_p.ny
                nb_p = 2 * grid_p.nx + 2 * grid_p.ny
            
                ni_u = grid_u.nx * grid_u.ny
                nb_u = 2 * grid_u.nx + 2 * grid_u.ny
            
                ni_v = grid_v.nx * grid_v.ny
                nb_v = 2 * grid_v.nx + 2 * grid_v.ny

                ni_uv = ni_u + ni_v
                nb_uv = nb_u + nb_v

                if num.pressure_velocity_coupling == 0
                    if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
                        nt = ni_uv + nb_uv
                        ncol_A = ni_uv + nb_uv
                    else
                        nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv 
                        ncol_A = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv
                    end

                    AuvL = spzeros(ncol_A, nt)
                    BuvL = spzeros(ncol_A, nt)
                    rhs_uv = zeros(ncol_A)  
                    # AϕL = spzeros(nt, nt)



                elseif num.pressure_velocity_coupling == 1
                    nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + (num.nLS + 1) * ni_p + nb_p
                    ncol_A = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + ni_p
                    # so 1 * ni + 1 * ni_p +nb + ni_p + nb
                    # u v Navier, pression     
                
                    # AuvL = spzeros(nt, nt)
                    # BuvL = spzeros(nt, nt)

                    AuvL = spzeros(ncol_A, nt)
                    BuvL = spzeros(ncol_A, nt)
                    rhs_uv = zeros(ncol_A)  

                elseif num.pressure_velocity_coupling == 2
                    nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + (num.nLS + 1) * ni_p + nb_p
                    ncol_A = nt
                    
                    # so 1 * ni + 1 * ni_p +nb + ni_p + nb
                    # u v Navier, pression     
                
                    # AuvL = spzeros(nt, nt)
                    # BuvL = spzeros(nt, nt)

                    AuvL = spzeros(ncol_A, nt)
                    BuvL = spzeros(ncol_A, nt)
                    rhs_uv = zeros(ncol_A)  

                elseif num.pressure_velocity_coupling == 3 #no BC for pressure
                   


                    if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
                        nt = ni_uv + nb_uv + ni_p
                        ncol_A = ni_uv + nb_uv + ni_p
                    else
                        nt = (num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + nb_uv + ni_p
                        ncol_A = nt               
                    end

                    AuvL = spzeros(ncol_A, nt)
                    BuvL = spzeros(ncol_A, nt)
                    rhs_uv = zeros(ncol_A)  
                    # AϕL = spzeros(nt, nt)


                elseif num.pressure_velocity_coupling == 4 # BC for pressure on interfaces (bubble)
                    
                    n_phase = 2
                    
                    nt = n_phase * ((num.nLS - num.nNavier + 1) * ni_uv + num.nNavier * ni_p + (num.nLS + 1) * ni_p ) + nb_uv

                    ncol_A = nt
                    
                    # so 1 * ni + 1 * ni_p +nb + ni_p + nb
                    # u v Navier, pression     
                
                    # AuvL = spzeros(nt, nt)
                    # BuvL = spzeros(nt, nt)

                    AuvL = spzeros(ncol_A, nt)
                    BuvL = spzeros(ncol_A, nt)
                    rhs_uv = zeros(ncol_A)  

                # elseif num.pressure

                end

            end

            #endregion Allocate for Navier case

            #endregion Allocations

            #TODO remove allocations


            # Adapt cell volume W for gradients 
            # cf 4/3 factor for Laplacian
            if num.laplacian == 1
                AvLcopy = copy(AvL)
                Lvm1_L,bc_Lvm1_L,bc_Lvm1_b_L=compute_divergence!(num, 
                # grid_p, 
                # grid_u, 
                grid_v, 
                op.opC_vL,
                AvLcopy, 
                # rhs_scal,
                # tmp_vec_p, #a0
                tmp_vec_1D_v0,
                tmp_vec_1D_v,
                Lvm1_L, 
                bc_Lvm1_L, 
                bc_Lvm1_b_L
                # tmp_vec_u0,
                # tmp_vec_v0,
                # tmp_vec_1D,
                # ls_advection
                )

                ApLcopy = copy(Ascal)
                Lpm1_L,bc_Lpm1_L,bc_Lpm1_b_L=compute_divergence!(num, 
                # grid_p, 
                # grid_u, 
                grid_p, 
                op.opC_pL,
                ApLcopy, 
                # rhs_scal,
                # tmp_vec_p, #a0
                tmp_vec_1D_p0,
                tmp_vec_1D_p,
                Lpm1_L, 
                bc_Lpm1_L, 
                bc_Lpm1_b_L
                # tmp_vec_u0,
                # tmp_vec_v0,
                # tmp_vec_1D,
                # ls_advection
                )

            end  
           
            if num.pressure_velocity_coupling == 0
                #TODO why this call without interface initialization ?
                if !navier
                    if (num.solve_solid == 1) && ns_solid_phase
                        _ = FE_set_momentum(
                            num, grid_u, op.opC_uS,
                            AuS, BuS,
                            iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                            true
                        )
                    
                    _ = FE_set_momentum(
                        num, grid_v, op.opC_vS,
                        AvS, BvS,
                        iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                        true
                    )
                    end
                else
                    _ = FE_set_momentum_coupled(
                        BC_int, num, grid_p, grid_u, grid_v,
                        op.opC_pS, op.opC_uS, op.opC_vS,
                        AuvS, BuvS,
                        iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                        iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                        true
                    )
                end
            elseif num.pressure_velocity_coupling > 1
                if num.pressure_velocity_coupling == 4
                    # Coupled resolution of u and v, going to two-phases for pressure
                    _ = FE_set_momentum_coupled_two_phases(
                    BC_int, num, grid_p, grid_u, grid_v,
                    op,
                    AuvL, BuvL,rhs_uv,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true,BC_pL,phL,phS
                    )
                else
                    if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
                        # _ = FE_set_momentum_coupled2_one_fluid(
                        #     BC_int, num, grid_p, grid_u, grid_v,
                        #     op.opC_pL, op.opC_uL, op.opC_vL,
                        #     AuvL, BuvL,rhs_uv,
                        #     diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_uL,
                        #     diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_vL,
                        #     cross_term_diffusion_bulk_d_dv_dx_dy,cross_term_diffusion_bulk_d_du_dy_dx,
                        #     cross_term_diffusion_bulk_d_dv_dx_dy_border,cross_term_diffusion_bulk_d_du_dy_dx_border,
                        #     rho_one_fluid_u,rho_one_fluid_v,
                        #     true,
                        #     BC_pL,phL
                        #     )
                    else

                        # Coupled resolution of u and v
                        _ = FE_set_momentum_coupled2(
                        BC_int, num, grid_p, grid_u, grid_v,
                        op.opC_pL, op.opC_uL, op.opC_vL,
                        AuvL, BuvL,rhs_uv,
                        iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                        iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                        true,BC_pL,phL
                        )
                    end
                end


            end


            #TODO remove alloc a0_p
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid_p))
            end
            
            
            if !advection && (num.solve_solid == 1)
            #call to set_heat! is there to set up the matrices for the heat equation. 
            #If the level-set is not advected, then after this call there is no need to update these matrices anymore

                _ = set_poisson(
                    BC_int, num, grid_p, a0_p, op.opC_pS, op.opC_uS, op.opC_vS,
                    AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                    true
                )

                set_heat!(
                    BC_int[1], num, grid_p, op.opC_TS, grid_p.LS[1].geoS, phS, num.θd, BC_TS, grid_p.LS[1].MIXED, grid_p.LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    op.opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, true, BC_int
                )
            end

            if test_laplacian
                num.exact_laplacian = test_laplacian_pressure(num,grid_v,phL.vD,op.opC_pL, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L)
                return
            end

            # Segregated resolution for u and v since the viscosity is assumed to be constant in a phase 
            if num.pressure_velocity_coupling == 0
                if !navier
                    _ = FE_set_momentum(
                        num, grid_u, op.opC_uL,
                        AuL, BuL,
                        iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                        true
                    )
                    _ = FE_set_momentum(
                        num, grid_v, op.opC_vL,
                        AvL, BvL,
                        iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                        true
                    )
                else
                    # Coupled resolution of u and v if navier BC is activated
                    _ = FE_set_momentum_coupled(
                        BC_int, num, grid_p, grid_u, grid_v,
                        op.opC_pL, op.opC_uL, op.opC_vL,
                        AuvL, BuvL,
                        iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                        iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                        true
                    )
                end

            elseif num.pressure_velocity_coupling > 1


                if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)
                        # _ = FE_set_momentum_coupled2_one_fluid(
                        #     BC_int, num, grid_p, grid_u, grid_v,
                        #     op.opC_pL, op.opC_uL, op.opC_vL,
                        #     AuvL, BuvL,rhs_uv,
                        #     diffusion_bulk_u, diffusion_LS_u, diffusion_border_u, Mum1, BC_uL,
                        #     diffusion_bulk_v, diffusion_LS_v, diffusion_border_v, Mvm1, BC_vL,
                        #     cross_term_diffusion_bulk_d_dv_dx_dy,cross_term_diffusion_bulk_d_du_dy_dx,
                        #     cross_term_diffusion_bulk_d_dv_dx_dy_border,cross_term_diffusion_bulk_d_du_dy_dx_border,
                        #     rho_one_fluid_u,rho_one_fluid_v,
                        #     true,
                        #     BC_pL,phL
                        #     )
                else
                    # Coupled resolution of u and v
                    _ = FE_set_momentum_coupled2(
                    BC_int, num, grid_p, grid_u, grid_v,
                    op.opC_pL, op.opC_uL, op.opC_vL,
                    AuvL, BuvL,rhs_uv,
                    iRe.*Lum1_L, iRe.*bc_Lum1_L, iRe.*bc_Lum1_b_L, Mum1_L, BC_uL,
                    iRe.*Lvm1_L, iRe.*bc_Lvm1_L, iRe.*bc_Lvm1_b_L, Mvm1_L, BC_vL,
                    true,BC_pL,phL
                    )
                end #one-fluid
                
            end

            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid_p))
            end

            if num.one_fluid_model == 0
                _ = set_poisson(
                    BC_int, num, grid_p, a0_p, op.opC_pL, op.opC_uL, op.opC_vL,
                    AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                    true
                )
            else
                _ = set_poisson_one_fluid(
                BC_int, num, grid_p, a0_p, op.opC_pL, op.opC_uL, op.opC_vL,
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )
            end

            if num.nLS>1 #monophasic
                #call to set_heat! is there to set up the matrices for the heat equation. 
                #If the level-set is not advected, then after this call there is no need to update these matrices anymore
                set_heat!(
                    BC_int[1], num, grid_p, op.opC_TL, grid_p.LS[1].geoL, phL, num.θd, BC_TL, grid_p.LS[1].MIXED, grid_p.LS[1].geoL.projection,
                    ATL, BTL,rhs_scal,
                    op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, true, BC_int
                )
            end
         
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    if heat_convection || electrolysis_convection
        NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y)
        # println(grid_p.LS[1].geoL.dcap[1,1,:])
    end

    # V0S = volume(grid_p.LS[end].geoS)
    # V0L = volume(grid_p.LS[end].geoL)

    if num.debug == "allocations_start"
        print("\n STOP allocations")
        return
    end

    #region restart    
    #  - file: decl_hdf5_test_02_r${rank}.h5
    #   when: $input=1
    #   read: [ reals, values ]

    PDI_status = @ccall "libpdi".PDI_multi_expose("restart"::Cstring,
    "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
    "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
    "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
    "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
    "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,
    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
    "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
    "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,        
    "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
    "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
    "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint,        
    C_NULL::Ptr{Cvoid})::Cint

    #endregion restart 

    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
    
    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
        "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
        "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
        "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,
        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,        
        "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
        "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint,        
        C_NULL::Ptr{Cvoid})::Cint

        #compute numerical radius 
        PDI_status = @ccall "libpdi".PDI_multi_expose("compute_radius"::Cstring,
        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "mesh_p_x"::Cstring, grid_p.x::Ptr{Cdouble}, PDI_OUT::Cint,
        "mesh_p_y"::Cstring, grid_p.y::Ptr{Cdouble}, PDI_OUT::Cint,
        "radius_vec"::Cstring, radius_pdi::Ptr{Cdouble}, PDI_INOUT::Cint,                             
        C_NULL::Ptr{Cvoid})::Cint
        num.current_radius = radius_pdi[1]

        
        print("\n radius pdi ", radius_pdi[1]) 

        # slice = levelset_p[nx//2,:]
        # x_1D = mesh_p_y[nx//2,:]

        # try
    #     radius_vertical = compute_radius_from_levelset_slice(grid_p.LS[num.iLSpdi].u[:,div(grid_p.nx,2)],
    #     grid_p.y[:,div(grid_p.nx,2)])
    #     # catch error
    #     # volume_cell = grid_p.LS[num.iLSpdi].geoS.cap[:,:,5] #bubble
    #     volume_cell = grid_p.LS[num.iLSpdi].geoL.cap[:,:,5] #drop

    #     center_of_mass_x, center_of_mass_y = calculate_centroid(grid_p.x, grid_p.y, volume_cell)

        
    #     indices_bubble_mass_center = find_slice_coord_bubble_mass_center(center_of_mass_x,center_of_mass_y,
    #     num,grid_p)
    #     radius_horizontal = compute_radius_from_levelset_slice(grid_p.LS[num.iLSpdi].u[indices_bubble_mass_center[1],:],
    #     grid_p.y[indices_bubble_mass_center[1],:]) 

    #     if isnothing(radius_vertical)
    #         print("\n error radius vertical")
    #         if isnothing(radius_horizontal)
    #             print("\n error radius horizontal and vertical")
    #         else
    #             num.current_radius = radius_horizontal
    #         end

    #     else
    #         if isnothing(radius_horizontal)
    #             print("\n error radius horizontal")
    #         else
    #             num.current_radius = max(radius_horizontal, radius_vertical)
    #         end
    #     end
    #     # end

    #     num.current_radius = maximum(compute_radii_from_slices(
    # grid_p.LS[num.iLSpdi].u,
    # x_grid::AbstractMatrix,
    # y_grid::AbstractMatrix,
    # slices::Vector{Tuple{Union{Int, Colon}, Union{Int, Colon}}};
    # center_x,center_y)

    
    compute_bubble_drop_radius(num, grid_p)


    if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) 
        # rise_velocity_y =0.0
        # PDI_status = @ccall "libpdi".PDI_multi_expose("post_processing_rising_bubble_first_share"::Cstring,
        # "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        # # "rho_one_fluid"::Cstring, rho_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "rho_one_fluid_u"::Cstring, rho_one_fluid_u::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "rho_one_fluid_v"::Cstring, rho_one_fluid_v::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "mu_one_fluid"::Cstring, mu_one_fluid::Ptr{Cdouble}, PDI_OUT::Cint,
        # "rise_velocity_y"::Cstring, rise_velocity_y::Ref{Cdouble}, PDI_OUT::Cint,  
        # "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
        # # "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "volume_cell"::Cstring, grid_p.LS[end].geoS.dcap[:,:,5]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # # "mesh_p_x"::Cstring, grid_p.x::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "mesh_p_y"::Cstring, grid_p.y::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "dcap_1"::Cstring, grid_p.LS[num.iLSpdi].geoS.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # # "dcap_2"::Cstring, grid_p.LS[num.iLSpdi].geoS.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # # "dcap_3"::Cstring, grid_p.LS[num.iLSpdi].geoS.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # # "dcap_4"::Cstring, grid_p.LS[num.iLSpdi].geoS.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # C_NULL::Ptr{Cvoid})::Cint

        update_one_fluid_density_viscosity(num,grid_p,grid_u,grid_v,volume_fraction,levelset_one_fluid,rho_one_fluid,
                                                    rho_one_fluid_u,rho_one_fluid_v,mu_one_fluid,tmp_vec_p0)

        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_timestep"::Cstring,
        # "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        # "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
        # "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # #region initial values
        # PDI_status = @ccall "libpdi".PDI_multi_expose("post_processing_rising_bubble"::Cstring,
        # #endregion initial values          
                                           
    end

    
    conservation = compute_conservation_mass(num,phL, grid_p ,grid_u, grid_v, rho_one_fluid)

    # Compute divergence of velocity
    Duv = op.opC_pL.AxT * vec1(phL.uD,grid_u) .+ op.opC_pL.Gx_b * vecb(phL.uD,grid_u) .+
    op.opC_pL.AyT * vec1(phL.vD,grid_v) .+ op.opC_pL.Gy_b * vecb(phL.vD,grid_v)
    for iLS in 1:num.nLS
        if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS]) #otherwise normal velocity null if no blowing
            Duv .+= op.opC_pL.Gx[iLS] * veci(phL.uD,grid_u,iLS+1) .+ 
                    op.opC_pL.Gy[iLS] * veci(phL.vD,grid_v,iLS+1)
        end
    end

    PDI_status = @ccall "libpdi".PDI_multi_expose("print_conservation"::Cstring,
    "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
    "conservation"::Cstring, conservation::Ref{Cdouble}, PDI_OUT::Cint,
    "velocity_divergence"::Cstring, Duv::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

    # PDI_status = @ccall "libpdi".PDI_multi_expose("check_conservation"::Cstring,
    # "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
    # "conservation"::Cstring, conservation::Ref{Cdouble}, PDI_OUT::Cint,
    # "velocity_divergence"::Cstring, Duv::Ptr{Cdouble}, PDI_OUT::Cint,
    # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
    # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
    # C_NULL::Ptr{Cvoid})::Cint



    simulation_finished = false

    #TODO variable time steps 

    #region time loop
    while (num.current_iter < num.max_iterations + 1) && (num.time < num.end_time) && (num.stop_simulation == 0) 

        #update time
        num.time += num.timestep_n
        num.current_iter += 1


        #region start iter

        #region Adapt timestep

        
        num.timestep_n = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
        printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e num.timestep_n : %.2e\n" num.CFL num.timestep_n num.timestep_n)
        print("\n num.stop_simulation start loop ",num.stop_simulation)
        #endregion adapt time
        
        if grid_p.LS[1].geoL.dcap[1,1,:] == 0.0
            printstyled(color=:red, @sprintf "\n Error operator null \n")
            return
        end

        # Print information at the start of temporal iteration and check definition of operators 
        # cf update_all_ls_data (true VS false)
        PDI_status = @ccall "libpdi".PDI_multi_expose("print_start_temporal_iteration"::Cstring,
        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
        # "levelset_p"::Cstring, grid_p.LS[num.index_levelset_pdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "dcap"::Cstring, grid_p.LS[num.index_levelset_pdi].geoL.dcap[:,:,:]::Ptr{Cdouble}, PDI_OUT::Cint,
        # grid_p.LS[1].geoL.dcap[1,1,:]
        C_NULL::Ptr{Cvoid})::Cint
 



        # if num.io_pdi>0

        #     #or permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1)) (3, 1, 2)
        #     try                            
        #         # num.iLSpdi = 1 # TODO all grid_p.LS                
        #         PDI_status = @ccall "libpdi".PDI_multi_expose("write_capacities"::Cstring,                    
        #         # "dcap"::Cstring, permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1))::Ptr{Cdouble}, PDI_OUT::Cint,
        #         "dcap_1"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
        #         "dcap_2"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
        #         "dcap_3"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
        #         "dcap_4"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint,

        #         C_NULL::Ptr{Cvoid})::Cint                           
        #     catch error
        #         printstyled(color=:red, @sprintf "\n PDI error \n")
        #         print(error)
        #         printstyled(color=:red, @sprintf "\n PDI error \n")
        #     end
        # end #if io_pdi

        #endregion start iter

        # PDI (IO)
        if electrolysis

            #TODO not necessary to expose everything now for ex only grid_p.LS ? and the rest later
            
            intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num, intfc_seg_num = convert_interfacial_D_to_segments(num,grid_p,phL.TD,1,2)
          
            #TODO when to write elec dat, ...
    
            if num.io_pdi>0
    
                try
                    # printstyled(color=:red, @sprintf "\n PDI test \n" )
            
            
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    
                    if num.electrical_potential > 0

                        # print("\n type of elec_cond ", typeof(elec_cond)," \n")
                        try
                            compute_grad_phi_ele!(num, grid_p, grid_u, grid_v,grid_u.LS[end], grid_v.LS[end], 
                            phL, 
                            # phS, 
                            op.opC_pL, 
                            # op.opC_pS, 
                            elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                        catch
                            @error("compute_grad_phi_ele!")
                            break
                        end
                    end
                    # #store in us, vs instead of Eus, Evs
                    # interpolate_grid_liquid!(grid_p,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)
    
                    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # C_NULL::Ptr{Cvoid})::Cint
                    
                    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
            
                    # print("\n before write \n ")
            
                    # num.iLSpdi = 1 # all grid_p.LS iLS = 1 # or all grid_p.LS ?
    
                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name
                    
                    # printstyled(color=:magenta, @sprintf "\n PDI write_data_start_loop %.5i \n" num.current_iter)
    
                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data_start_loop"::Cstring,
                    "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,                  
                    "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
                    "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
                    "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
                    "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
                    "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
                    "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
                    "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint
            
                catch error
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                    print(error)
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                end
    
            end #if io_pdi
        end #if electrolysis
        

        if !stefan
            grid_p.V .= speed #*ones(grid_p.ny, grid_p.nx)
        end

        #region Heat equation
        # Solve heat equation
        if heat
            if heat_solid_phase 
                kill_dead_cells!(phS.T, grid_p, grid_p.LS[1].geoS)
                veci(phS.TD,grid_p,1) .= vec(phS.T)
                set_heat!(
                    BC_int[1], num, grid_p, op.opC_TS, grid_p.LS[1].geoS, phS, num.θd, BC_TS, grid_p.LS[1].MIXED, grid_p.LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    op.opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs_scal
                phS.T .= reshape(veci(phS.TD,grid_p,1), grid_p)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid_p, grid_p.LS[1].geoL)
                veci(phL.TD,grid_p,1) .= vec(phL.T)
                set_heat!(
                    BC_int[1], num, grid_p, op.opC_TL, grid_p.LS[1].geoL, phL, num.θd, BC_TL, grid_p.LS[1].MIXED, grid_p.LS[1].geoL.projection,
                    ATL, BTL,rhs_scal,
                    op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs_scal
                phL.T .= reshape(veci(phL.TD,grid_p,1), grid_p)

            end
        end # #if heat
        #endregion heat equation

        
        #region Electrolysis 
        if electrolysis && (num.solve_potential == 1 ) && (num.solve_species == 1)
            if electrolysis_liquid_phase

                #region Electrolysis: Poisson  
                if num.electrical_potential > 0
                    solve_poisson_loop!(num, grid_p, grid_u, grid_v, op, Ascal, rhs_scal,F_residual,
                        tmp_vec_p,tmp_vec_p0,tmp_vec_p1, a1_p, BC_phi_ele, phL, phS,elec_cond,
                        elec_condD, tmp_vec_u, tmp_vec_v, i_butler, ls_advection, heat)
                end
                #endregion Poisson
    

               
                #TODO check BC not overwritten by different scalars
                #TODO check ls advection true when num.n scalars

                # print("\n vecb_L",vecb_L(phL.phi_eleD, grid_p))
                # print("\n phL.phi_ele[:,1]",phL.phi_ele[:,1])

                #TODO phL.phi_ele[:,1] or vecb_L(phL.phi_eleD, grid_p)
                # we suppose phi(x=0)=... cf Khalighi
                #but here BC
                
                # New start scalar loop

                #region velocity
                # TODO
                interpolate_interface_velocity!(phL,grid_u,grid_v)

                #endregion velocity

                #region Impose velocity
                if imposed_velocity == "zero"
                    phL.u .= 0.0
                    phL.v .= 0.0
                    phL.uD .= 0.0
                    phL.vD .= 0.0
                    phL.p .= 0.0
                    phL.pD .= 0.0


                elseif imposed_velocity == "constant"
                        
                        #Required to modify whole uD vD
                        phL.u .= 0.0
                        phL.v .= BC_vL.bottom.val    
                        phL.uD .= 0.0
                        phL.vD .= BC_vL.bottom.val    

                        phL.pD .= 0.0

                        if ((num.current_iter-1)%show_every == 0) 
                            printstyled(color=:red, @sprintf "\n Imposed velocity v min %.2e max %.2e\n" minimum(phL.vD) maximum(phL.vD))
                            printstyled(color=:red, @sprintf "\n Imposed velocity u min %.2e max %.2e\n" minimum(phL.uD) maximum(phL.uD))
                        end

                elseif imposed_velocity == "Poiseuille_bottom_top"

                    vPoiseuille = Poiseuille_fmax.(grid_v.x,num.v_inlet,num.L0)

                    
                    #Required to modify whole uD vD
                    phL.u .= 0.0
                    phL.v .= vPoiseuille 
 
                    phL.uD .= 0.0 #Dirichlet and Neumann
                    # phL.vD .= BC_vL.bottom.val    
                    # vecb...TODO
                    # ...
                    
                elseif imposed_velocity == "radial"
                    impose_radial_velocity(phS,phL,num)
                end #imposed_velocity
                #endregion Impose velocity

                #region Update current
                if num.electrolysis_reaction == "Butler_no_concentration"
                    
                    update_electrical_current_from_Butler_Volmer!(num,grid_p,heat,phL.phi_eleD,i_butler;phL.T)

                    # #TODO dev multiple levelsets
                    # if heat
                    #     i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid_p),
                    #     num.phi_ele1,num.Ru,phL.T)
                    # else
                    #     if num.nLS == 1
                    #         i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid_p),
                    #         num.phi_ele1,num.Ru,num.temperature0)
                    #     # else
                    #         #imposed by LS 2
                    #         # iLS_elec = 2
                    #         # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid_p,iLS_elec+1),
                    #         # num.phi_ele1,num.Ru,num.temperature0)
                    #     end
                    # end   
                end
                #endregion Update current


                #region Scalar transport: update boundary conditions
                if num.nb_transported_scalars>0

                    # printstyled(color=:magenta, @sprintf "\n num.nb_transported_scalars %.5i " num.nb_transported_scalars)

                    for iscal=1:num.nb_transported_scalars
                        # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[1].geoL,0.0) 
                        # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[1].geoL,num.concentration0[iscal]) 

                        if num.kill_dead_cells == 0
                            @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,0.0)
                        else
                            @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid_p, grid_p.LS[end].geoL,num.concentration0[iscal])
                        end

                        @views veci(phL.trans_scalD[:,iscal],grid_p,1) .= vec(phL.trans_scal[:,:,iscal])

                        if num.electrolysis_reaction == "Butler_no_concentration" && num.nLS == 1

                            #BC for LS 2 in scalar transport : done in scalar loop

                            if iscal==1 || iscal==2
                                inv_stoechiometric_coeff = -1.0/2.0 #H2 and KOH
                            elseif iscal == 3
                                inv_stoechiometric_coeff = 1.0 #H2O consummed
                            end

                            # print("\n BC_trans_scal[iscal].left.val ", BC_trans_scal[iscal].left.val )
                            # print("\n BC_trans_scal[iscal].left.val ", i_butler./(num.Faraday*num.diffusion_coeff[iscal])*inv_stoechiometric_coeff )

                            # BC at left wall
                            # -(-/i_butler) because i=-lambda grad phi and BC at left: -e_x
                            BC_trans_scal[iscal].left.val = i_butler./(num.Faraday*num.diffusion_coeff[iscal])*inv_stoechiometric_coeff

                            # print("\n")
                            # print("\n left BC ", BC_trans_scal[iscal].left.val)

                            # for testn in 1:grid_p.ny
                            #     printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i border %.5e\n" testn grid_p.ny-testn+1 vecb_L(phL.trans_scalD[:,iscal], grid_p)[testn])
                            # end
                        elseif num.electrolysis_reaction == "Butler_no_concentration" && num.nLS == 1
                            if iscal==1 || iscal==2
                                inv_stoechiometric_coeff = -1.0/2.0 #H2 and KOH
                            elseif iscal == 3
                                inv_stoechiometric_coeff = 1.0 #H2O consummed
                            end

                            BC_trans_scal[iscal].bottom.val = i_butler./(num.Faraday*num.diffusion_coeff[iscal])*inv_stoechiometric_coeff


                        end
                    end

                    #TODO convection_Cdivu BC divergence
                    #TODO check num.nb_transported_scalars>1

                    if ((num.current_iter-1)%show_every == 0) 
                        # printstyled(color=:cyan, @sprintf "\n before scalar transport \n")
                        # print_electrolysis_statistics(num,grid_p,phL)
                        
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                           
                        C_NULL::Ptr{Cvoid})::Cint
                    end
                
                    # printstyled(color=:red, @sprintf "\n levelset: before scalar_transport\n")
                    # println(grid_p.LS[1].geoL.dcap[1,1,:])

                            
                    #TODO better alloc
                    # geo = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_u = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_v = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    # laps = set_matrices!(
                    #     num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
                    #     op.opC_pL, op.opC_uL, op.opC_vL,
                    #     periodic_x, periodic_y)

                    if electrolysis_advection

                    # if imposed_velocity == "zero" && num.current_iter ==16
                        #    num.ϵwall = 0.0
                        #   print("\n changed eps wall")
                    # end
                        #printstyled(color=:cyan, @sprintf "\n epsilon %.2e %.2e \n" num.ϵ num.ϵwall )

                        update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y)
                    end

                    
                    # printstyled(color=:red, @sprintf "\n return before debug mem\n")
                    # return

                    # if num.nLS ==1 #TODO dev multiple levelsets


                      #PDI (IO)      
        
                    if num.io_pdi>0

                        #or permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1)) (3, 1, 2)
                        try                            
                            # num.iLSpdi = 1 # TODO all grid_p.LS                
                            PDI_status = @ccall "libpdi".PDI_multi_expose("write_capacities"::Cstring,                    
                            # "dcap"::Cstring, permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1))::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_1"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_2"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_3"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_4"::Cstring, grid_p.LS[num.iLSpdi].geoL.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint,
                            C_NULL::Ptr{Cvoid})::Cint                           
                        catch error
                            printstyled(color=:red, @sprintf "\n PDI error \n")
                            print(error)
                            printstyled(color=:red, @sprintf "\n PDI error \n")
                        end
                    end #if io_pdi


                   
                    # # #TODO better alloc
                    # geo = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_u = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_v = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    # laps = set_matrices!(
                    #     num, grid_p, geo, grid_u, geo_u, grid_v, geo_v,
                    #     op.opC_pL, op.opC_uL, op.opC_vL,
                    #     periodic_x, periodic_y)

                    #interface term in rhs comes from op.opC_TL
                    #TODO variable CL
                    
                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_iso"::Cstring,
                    "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                    "levelset_iso"::Cstring, grid_p.LS[num.iLSpdi].iso::Ptr{Cdouble}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint

                    scalar_transport!(num, grid_p, grid_u, grid_v,
                    op.opC_TL, #op
                    op.opL, #op_conv
                    phL, 
                    BC_trans_scal,
                    BC_int,                     
                    Ascal,
                    Bscal,
                    all_CUTCT,
                    rhs_scal,
                    tmp_vec_p, #used to store a0 for rhs of interfacial value
                    tmp_vec_u,
                    tmp_vec_v,
                    mass_transfer_rate,
                    periodic_x, 
                    periodic_y, 
                    electrolysis_convection,                 
                    ls_advection)

                    print("\n num.stop_simulation after scalar ",num.stop_simulation)


                    PDI_status = @ccall "libpdi".PDI_multi_expose("check_concentrations"::Cstring,
                        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,            
                        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint

                    print("\n BC_trans_scal ",BC_trans_scal)

                    # PDI_status = @ccall "libpdi".PDI_multi_expose("check_concentrations"::Cstring,
                
                    # scalar_transport!(BC_trans_scal, num, grid_p, , grid_p.LS[1].geoL, phL, num.concentration0,
                    # grid_p.LS[1].MIXED, grid_p.LS[1].geoL.projection, op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    # periodic_x, periodic_y, electrolysis_convection, ls_advection, BC_int, num.diffusion_coeff,Ascal,Bscal,all_CUTCT,rhs_scal)

                    # else
                    #     printstyled(color=:red, @sprintf "\n TODO multiple LS \n" )

                    # end

                    # scalar_transport_2!(BC_trans_scal, num, grid_p, op.opC_TL, grid_p.LS[1].geoL, phL, num.concentration0,
                    # grid_p.LS[1].MIXED, grid_p.LS[1].geoL.projection, op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    # periodic_x, periodic_y, electrolysis_convection, true, BC_int, num.diffusion_coeff)

                    if ((num.current_iter-1)%show_every == 0) 
                        # printstyled(color=:cyan, @sprintf "\n after scalar transport \n")
                        # print_electrolysis_statistics(num,grid_p,phL)
                        
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                          
                        C_NULL::Ptr{Cvoid})::Cint


                        # scal_error=0.0
                        # for iscal in 1:num.nb_transported_scalars
                        #     # printstyled(color=:cyan, @sprintf "\n after scalar transport %.2i %.2e \n" iscal maximum(abs.(phL.trans_scal[:,:,iscal])))
                        #     # printstyled(color=:cyan, @sprintf "\n after scalar transport %.2i %.2e \n" iscal maximum(abs.(phL.trans_scalD[:,iscal])))
                        #     # scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-num.concentration0[iscal])./num.concentration0[iscal])
                        #     # scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid_p).-num.concentration0[iscal])./num.concentration0[iscal])
                        #     # scal_error = max(scal_error_bulk,scal_error_border,scal_error)
                        # end

                        # printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e num.CFL %.2e\n" scal_error num.v_inlet*num.timestep_n/grid_p.dx[1,1])


                    end

                    if imposed_velocity != "none" && num.debug== "scalar_testing"
                        scal_error=0.0
                        for iscal in 1:num.nb_transported_scalars

                            # print("\n maximum ",maximum(phL.trans_scalD[:,iscal]), )
                            # printstyled(color=:cyan, @sprintf "\n error after scalar transport max %.2e min %.2e\n" maximum(phL.trans_scalD[:,iscal]) minimum(phL.trans_scalD[:,iscal]))

                            scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-num.concentration0[iscal])./num.concentration0[iscal])
                            scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid_p).-num.concentration0[iscal])./num.concentration0[iscal])
                            scal_error = max(scal_error_bulk,scal_error_border,scal_error)

                        end

                        printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e num.CFL %.2e\n" scal_error num.v_inlet*num.timestep_n/grid_p.dx[1,1])

                        # printstyled(color=:red, @sprintf "\n Poiseuille \n")

                        # # Check the velocity field before the scalar transport
                        # test_Poiseuille(num,phL.vD,grid_v)

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

                        # compute_grad_p!(num,grid_p, grid_u, grid_v, phL.pD, op.opC_pL, op.opC_uL, op.opC_vL)

                    end

                end #num.nb_transported_scalars>0

                #endregion Scalar transport

                # if imposed_velocity =="none"
                #     printstyled(color=:red, @sprintf "\n after scalar transport \n")

                #     # Check the velocity field before the scalar transport
                #     test_Poiseuille(num,phL,grid_v)
                    
                # end
                
   
            end #if electrolysis_liquid_phase        
        end #if electrolysis
        #endregion Electrolysis



        
        #PDI (IO)      
        if num.io_pdi>0

            try
                # printstyled(color=:red, @sprintf "\n PDI test \n" )
        
            
                # phi_array=phL.phi_ele #do not transpose since python row major
                
                # Compute electrical current, interpolate velocity on scalar grid_p
                #                     if num.electrical_potential>0 compute_grad_phi_ele!(num, grid_p, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
        
                # if electrolysis && num.nb_transported_scalars>1
                #     if heat 
                #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
                #     else
                #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
                #     end
                # else 
                #     elec_cond = ones(grid_p)
                #     printstyled(color=:green, @sprintf "\n conductivity one")

                # end 

                # phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude


                #store in us, vs instead of Eus, Evs
                # interpolate_grid_liquid!(grid_p,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

                # printstyled(color=:red, @sprintf "\n test i current\n" )
                # print("\n phi ",phL.phi_ele[64,127]," ",phL.phi_ele[64,128]," ",(phL.phi_ele[64,128]-phL.phi_ele[64,127])/grid_p.dx[64,128]," ",(phL.phi_ele[64,128]-phL.phi_ele[64,127])/grid_p.dx[64,127]*elec_cond[64,127], " cond ", elec_cond[64,127]," \n")


                # tmp_vec_p .*= elec_cond
                # tmp_vec_p0 .*= elec_cond


                # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # C_NULL::Ptr{Cvoid})::Cint
            
                interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
                    
                # num.iLSpdi = 1 # TODO all grid_p.LS

                # Exposing data to PDI for IO    
                # if writing "D" array (bulk, interface, border), add "_1D" to the name

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
                "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
                "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,            
                "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                
                "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                C_NULL::Ptr{Cvoid})::Cint
        
                
                #TODO debug with volume fraction

                # A = zeros(gv.ny, gv.nx+1)
                # for jplot in 1:gv.ny
                #     for iplot in 1:gv.nx+1
                #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                #     pII = lexicographic(II, grid_p.ny + 1)
                #     A[jplot,iplot] =  1 ./ op.opC_vL.iMx.diag[pII]
                #     end
                # end
            

                # print("\n after write \n ")
                    
                # printstyled(color=:red, @sprintf "\n PDI test end\n" )
        
            catch error
                printstyled(color=:red, @sprintf "\n PDI error \n")
                print(error)
                printstyled(color=:red, @sprintf "\n PDI error \n")
            end
        end #if io_pdi

        #region compute mass transfer
        if num.solve_Navier_Stokes_liquid_phase == 1 && num.phase_change_method >0
            # num.status = 
            compute_mass_transfer_rate_main!(num, grid_p, grid_u, grid_v, op, phL, phS, BC_int, electrolysis, electrolysis_phase_change_case, 
            periodic_x, periodic_y, λ, Vmean, num.iLSpdi, mode_2d, show_every, 
            mass_transfer_rate,mass_transfer_rate_vec1,
            mass_transfer_rate_vecb,mass_transfer_rate_veci, mass_transfer_rate_redistributed, tmp_vec_p, tmp_vec_p0, tmp_vec_p1,    
            nb_gaz_acceptors, volume_fraction, interface_length)
        #endregion compute mass transfer
        end

        #region Navier-Stokes
        if num.solve_Navier_Stokes_liquid_phase == 1
            if num.time < num.nucleation_time
                printstyled(color=:red, @sprintf "\n Navier-Stokes not solved")

                navier_stokes = false
            else
                printstyled(color=:red, @sprintf "\n Navier-Stokes solved")
                navier_stokes = true
            end

            if navier_stokes


                # Adapt cell volume W for gradients 
                # cf 4/3 factor for Laplacian
                if num.laplacian == 1
                    AvLcopy = copy(AvL) #does not work if reset A after operations in compute_divergence!
                    Lvm1_L,bc_Lvm1_L,bc_Lvm1_b_L=compute_divergence!(num, 
                    # grid_p, 
                    # grid_u, 
                    grid_v, 
                    op.opC_vL,
                    AvLcopy, 
                    # rhs_scal,
                    # tmp_vec_p, #a0
                    tmp_vec_1D_v0,
                    tmp_vec_1D_v,
                    Lvm1_L, 
                    bc_Lvm1_L, 
                    bc_Lvm1_b_L
                    # tmp_vec_u0,
                    # tmp_vec_v0,
                    # tmp_vec_1D,
                    # ls_advection
                    )
                    print("\nbefore pressure")
                    II = CartesianIndex(div(grid_v.ny,2),1)
                    pII = lexicographic(II,grid_v.ny)
                    print("\nLv[pII,:] ",Lvm1_L[pII,:])
                    print("\bc_Lvm1_b[pII,:] ",bc_Lvm1_b_L[pII,:])
                
                    ApLcopy = copy(Ascal)
                    Lpm1_L,bc_Lpm1_L,bc_Lpm1_b_L=compute_divergence!(num, 
                    # grid_p, 
                    # grid_u, 
                    grid_p, 
                    op.opC_pL,
                    ApLcopy, 
                    # rhs_scal,
                    # tmp_vec_p, #a0
                    tmp_vec_1D_p0,
                    tmp_vec_1D_p,
                    Lpm1_L, 
                    bc_Lpm1_L, 
                    bc_Lpm1_b_L
                    # tmp_vec_u0,
                    # tmp_vec_v0,
                    # tmp_vec_1D,
                    # ls_advection
                    )

                end  

                if num.pressure_velocity_coupling == 3

                    # num.iLSpdi = 1 # TODO all grid_p.LS                
                    PDI_status = @ccall "libpdi".PDI_multi_expose("print_capacities"::Cstring,                    
                    # "dcap"::Cstring, permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1))::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_1"::Cstring, grid_u.LS[num.iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_2"::Cstring, grid_u.LS[num.iLSpdi].geoL.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_3"::Cstring, grid_u.LS[num.iLSpdi].geoL.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_4"::Cstring, grid_u.LS[num.iLSpdi].geoL.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint                           

                    PDI_status = @ccall "libpdi".PDI_multi_expose("print_capacities"::Cstring,                    
                    # "dcap"::Cstring, permutedims(grid_p.LS[num.iLSpdi].geoL.dcap, (3, 2, 1))::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_1"::Cstring, grid_v.LS[num.iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_2"::Cstring, grid_v.LS[num.iLSpdi].geoL.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_3"::Cstring, grid_v.LS[num.iLSpdi].geoL.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "dcap_4"::Cstring, grid_v.LS[num.iLSpdi].geoL.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint   

                    # print("\n cap_1 ",grid_u.LS[num.iLSpdi].geoL.dcap[1,:,1])
                    # print("\n cap_1 ",grid_u.LS[num.iLSpdi].geoL.dcap[2,:,1])

                    # print("\n cap_1 ",grid_v.LS[num.iLSpdi].geoL.dcap[1,:,1])
                    # print("\n cap_1 ",grid_v.LS[num.iLSpdi].geoL.dcap[2,:,1])

                    # print("\n cap_2 ",grid_v.LS[num.iLSpdi].geoL.dcap[1,:,2])
                    # print("\n cap_2 ",grid_v.LS[num.iLSpdi].geoL.dcap[2,:,2])

                    # print("\n cap_3 ",grid_v.LS[num.iLSpdi].geoL.dcap[1,:,3])
                    # print("\n cap_3 ",grid_v.LS[num.iLSpdi].geoL.dcap[2,:,3])

                    # for u grid_p
                    # cap 1 same 
                    # cap_1 [0.0, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5]


        

                end


                # if !advection
                #     @time no_slip_condition!(num, grid_p, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
                #     # grid_u.V .= num.Δ / (1 * num.timestep_n)
                #     # grid_v.V .= 0.0
                # end

                # Pressure-velocity coupling

                if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) 

                    # interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)

                    tmp_vec_p  .= 0.0
                    tmp_vec_p0 .= 0.0

                    for j = 1:grid_p.ny
                    for i = 1:grid_p.nx
                        tmp_vec_p[j,i] =(phL.u[j,i]+phL.u[j,i+1])/2
                        tmp_vec_p0[j,i]=(phL.v[j,i]+phL.v[j+1,i])/2
                    end
                    end


                    update_one_fluid_density_viscosity(num,grid_p,grid_u,grid_v,volume_fraction,levelset_one_fluid,rho_one_fluid,
                                                        rho_one_fluid_u,rho_one_fluid_v,mu_one_fluid,tmp_vec_p0)

                    total_interface_length = compute_interface_length!(num, grid_p, 1, interface_length)
                    # MIXED =

                    nb_levelsets = num.nLS
                    num.nLS = 0 #1 #deactivate cut-cell for one-fluid model
                    
                    #region deactivate LS
                    empty_capacities = vcat(zeros(7), zeros(4))
                    full_capacities = vcat(ones(7), 0.5.*ones(4))

                    for grid_iter in [grid_p,grid_u,grid_v]

                        grid_iter.LS[1].u .= 1.0
                        grid_iter.LS[end].u .= 1.0


                        # for j in 1:grid_iter.ny
                        #     for i in 1:grid_iter.nx
                        #         II = CartesianIndex(j,i)
                        #         # grid_iter.LS[end].geoL.cap[II,:] .= full_capacities
                        #         # grid_iter.LS[end].geoS.cap[II,:] .= empty_capacities
                        #     end
                        # end
                    end

                    # grid_u.LS[end].geoL.cap[:,:,:] .= full_capacities
                    # grid_u.LS[end].geoS.cap[:,:,:] .= empty_capacities

                    # grid_v.LS[end].geoL.cap[:,:,:] .= full_capacities
                    # grid_v.LS[end].geoS.cap[:,:,:] .= empty_capacities

                    #endregion  deactivate LS

                    NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y,true,true) 


                    geoL = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                    geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    #region reset centroids 
                    for grid_iter in [grid_p,grid_u,grid_v]
                        for II in grid_iter.ind.inside
                            grid_iter.LS[1].geoL.centroid[II] = Point(0.0,0.0)
                        end
                    end
                    #endregion reset centroids 
                    #TODO reactivate centroids ?
                    # TODO update density

                    # if num.pressure_velocity_coupling == 0 

                    #region update LS

                    #endregion update LS






                    #pbm mass_transfer_rateL
                    PDI_status = @ccall "libpdi".PDI_multi_expose("check_mass_transfer_rate_NS"::Cstring,
                    # "conservation"::Cstring, conservation::Ref{Cdouble}, PDI_OUT::Cint,
                    "mass_transfer_rate"::Cstring, mass_transfer_rate::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint



                    #region update LS for convection (bool=true)+ one fluid

                    # At every iteration, update_all_ls_data is called twice, once inside run.jl 
                    # and another one (if there's advection of the levelset) inside set_heat!. 
                    # The difference between both is a flag as last argument, inside run.jl is implicitly defined 
                    # as true and inside set_heat! is false. If you're calling your version of set_heat! several times, 
                    # then you're calling the version with the flag set to false, but for the convective term it has to be set to true.
                    # The flag=true, the capacities are set for the convection, the flag=false they are set for the other operators

                    if advection
                        # update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y, true) 
                        update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y, true,true) 
                        if num.convection == 0
                            # op.opL is op_conv
                            set_convection_preallocated!(num, grid_p, geoL[end], grid_u, grid_u.LS, grid_v, grid_v.LS, phL.u, phL.v, op.opL,
                            phL, BC_uL, BC_vL,op.opC_pL, op.opC_uL, op.opC_vL,
                            velocity_and_BC_convection_u_x ,
                            velocity_and_BC_convection_u_y ,
                            velocity_and_BC_convection_v_x ,
                            velocity_and_BC_convection_v_y)

                            # Mm1_L .= op.opC_pL.M
                            # Mum1_L .= op.opC_uL.M
                            # Mvm1_L .= op.opC_vL.M
                        else


                        end

                        
                    end
                    #endregion update LS for convection (bool=true)+ one fluid


                    #region update LS for other operators than convection (bool=false)+ one fluid
                    # update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)
                    update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y, false,true)

                    #endregion



                    if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) 
                        if num.activate_interface == 0
                            volumic_surface_tension_u .= 0.0
                            volumic_surface_tension_v .= 0.0
                        else
                            if num.surface_tension == 0
                                compute_surface_tension_VOF!(num,grid_p, grid_u, grid_v, op.opC_pL, op.opC_uL, op.opC_vL, 
                                volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)
                            elseif num.surface_tension == 1
                                # compute_surface_tension_LS!(num,grid_p, grid_u, grid_v, opC_p, opC_u, opC_v, 
                                # volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0)
                                compute_surface_tension_LS!(num,grid_p, grid_u, grid_v, op.opC_pL, op.opC_uL, op.opC_vL, 
                                volume_fraction,levelset_one_fluid,volumic_surface_tension_u,volumic_surface_tension_v,tmp_vec_p,tmp_vec_p0,
                                levelset_1D, levelset_heavyside_2D, 
                                tmp_vec_u0,tmp_vec_v0, #normal_and_dirac_u, normal_and_dirac_v,
                                tmp_vec_u1,tmp_vec_v1,#normal_u, normal_v, 
                                tmp_vec_u,tmp_vec_v,#curvature_u, curvature_v
                                )
                            # elseif num. TODO HF                               
                            end

                        end
                    end
                    
                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_one_fluid_surface_tension_concise"::Cstring,
                    "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,             
                    "volumic_surface_tension_u"::Cstring, volumic_surface_tension_u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "volumic_surface_tension_v"::Cstring, volumic_surface_tension_v::Ptr{Cdouble}, PDI_OUT::Cint,              
                    C_NULL::Ptr{Cvoid})::Cint


                    #BC TODO

                    # @error("\n check temporal terms")

                    # Mum1_L is put in B matrix that multiplies v
                    # print("\n advection ", ns_advection, " adv ",advection)
                    one_fluid_NS_ls_advection = true #update matrix


                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L,
                    Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = solve_one_fluid_NS!(
                    time_scheme, BC_int,
                    num, grid_p, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    op.opC_pL, op.opC_uL, op.opC_vL, op.opL,
                    # op.opC_TL,
                    AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,rhs_uv,
                    Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                    Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                    periodic_x, periodic_y, ns_advection, one_fluid_NS_ls_advection, num.current_iter, Ra, navier,
                    volume_fraction,
                    levelset_one_fluid,
                    rho_one_fluid,
                    mu_one_fluid,
                    rho_one_fluid_u,
                    # mu_one_fluid_u,
                    rho_one_fluid_v,
                    # mu_one_fluid_v,
                    volumic_surface_tension_u,
                    volumic_surface_tension_v,
                    convection_u,convection_v,
                    viscosity_coeff_for_du_dx ,
                    viscosity_coeff_for_du_dy ,
                    viscosity_coeff_for_dv_dx ,
                    viscosity_coeff_for_dv_dy,
                    # velocity_and_BC_convection_u_x ,
                    # velocity_and_BC_convection_u_y ,
                    # velocity_and_BC_convection_v_x ,
                    # velocity_and_BC_convection_v_y ,
                    tmp_vec_p,
                    tmp_vec_p0,            
                    rhs_phi,
                    pres_free_surfaceL,jump_mass_transfer_rateL,mass_transfer_rate )  

                    #region ciprianoMulticomponentDropletEvaporation2024
                    if extend_liquid_velocity

                        # num.phase_change_currently_activated == 1 
                        phase_change_currently_activated = 0

                        # no need to recompute Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L,
                        # Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L

                        # Mm1_L_ext_vel, Mum1_L_ext_vel, Mvm1_L_ext_vel same as Mm1_L, Mum1_L, Mvm1_L
                        # because :
                        # Mm1_L = copy(op.opC_pL.M)
                        # Mum1_L = copy(op.opC_uL.M)
                        # Mvm1_L = copy(op.opC_vL.M)

                        Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L,
                        Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, 
                        Mm1_L, Mum1_L, Mvm1_L, 
                        Cum1L_ext_vel, Cvm1L_ext_vel = solve_one_fluid_NS_no_phase!(
                        time_scheme, BC_int,
                        num, grid_p, geoL, grid_u, geo_uL, grid_v, geo_vL, 
                        # phL,
                        p_ext_vel, pD_ext_vel, phi_ext_vel, u_ext_vel, v_ext_vel, u_predictionD_ext_vel, v_predictionD_ext_vel, uD_ext_vel, vD_ext_vel, u_prediction_ext_vel, v_prediction_ext_vel, uT_ext_vel,
                        pres_grad_x, pres_grad_y,                    
                        phase_change_currently_activated,
                        BC_uL, BC_vL, BC_pL,
                        op.opC_pL, op.opC_uL, op.opC_vL, op.opL,
                        # op.opC_TL,
                        AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,rhs_uv,
                        Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, 
                        Lum1_L, bc_Lum1_L, bc_Lum1_b_L, 
                        Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                        Cum1L_ext_vel, Cvm1L_ext_vel, 
                        Mum1_L, Mvm1_L,
                        periodic_x, periodic_y, ns_advection, advection, num.current_iter, Ra, navier,
                        volume_fraction,
                        levelset_one_fluid,
                        rho_one_fluid,
                        mu_one_fluid,
                        rho_one_fluid_u,
                        rho_one_fluid_v,
                        volumic_surface_tension_u,
                        volumic_surface_tension_v,
                        convection_u,convection_v,
                        viscosity_coeff_for_du_dx ,
                        viscosity_coeff_for_du_dy ,
                        viscosity_coeff_for_dv_dx ,
                        viscosity_coeff_for_dv_dy,             
                        tmp_vec_p,
                        tmp_vec_p0,            
                        rhs_phi,
                        pres_free_surfaceL,jump_mass_transfer_rateL,mass_transfer_rate )  

                        PDI_status = @ccall "libpdi".PDI_multi_expose("velocity_extension"::Cstring,                    
                        "u_ext_vel"::Cstring, u_ext_vel::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_ext_vel"::Cstring, v_ext_vel::Ptr{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint

                    end
                    #endregion ciprianoMulticomponentDropletEvaporation2024


                    print("\n num.stop_simulation after NS ",num.stop_simulation)


                    #region reactivate cut-cell for one-fluid model
                    num.nLS = nb_levelsets 

                    grid_p.LS[1].u .= levelset_one_fluid
                    grid_p.LS[end].u .= levelset_one_fluid

                    NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y) 

                    # TODO check after activation NS after nucleation

                    # geoL = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    # reactivate centroids
                    # printstyled(color=:red, @sprintf "\n reactivate centroids \n")

                    # x_centroid = grid_p.x .+ getproperty.(grid_p.LS[1].geoL.centroid, :x) .* grid_p.dx #geoS
                    # display(x_centroid)

                    # y_centroid = grid_p.y .+ getproperty.(grid_p.LS[1].geoL.centroid, :y) .* grid_p.dy #geoS
                    # display(y_centroid)

                    # display(levelset_one_fluid)
                    #endregion reactivate cut-cell for one-fluid model


                else

                    if num.pressure_velocity_coupling == 0 
                        if ns_solid_phase
                            geoS = [grid_p.LS[iLS].geoS for iLS in 1:num._nLS]
                            geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:num._nLS]
                            geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:num._nLS]
                            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                                time_scheme, BC_int,
                                num, grid_p, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                                BC_uS, BC_vS, BC_pS,
                                op.opC_pS, op.opC_uS, op.opC_vS, op.opS,
                                AuS, BuS, AvS, BvS, AϕS, AuvS, BuvS,
                                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                                Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                                periodic_x, periodic_y, ns_advection, advection, num.current_iter, Ra, navier,pres_free_surfaceS,jump_mass_transfer_rateS,mass_transfer_rateS
                            )
                        end
                        if ns_liquid_phase
                            geoL = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                            # Mum1_L is put in B matrix that multiplies v

                            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                                time_scheme, BC_int,
                                num, grid_p, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                                BC_uL, BC_vL, BC_pL,
                                op.opC_pL, op.opC_uL, op.opC_vL, op.opL,
                                AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,
                                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                                Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                                periodic_x, periodic_y, ns_advection, advection, num.current_iter, Ra, navier,pres_free_surfaceL,jump_mass_transfer_rateL,mass_transfer_rate
                            )
                            # if num.current_iter == 1
                            #     phL.u .= -0.5 .* grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy
                            #     phL.v .= 0.5 .* grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
                            #     phL.u[grid_u.LS[1].SOLID] .= 0.0
                            #     phL.v[grid_v.LS[1].SOLID] .= 0.0
                            # end
                            # linear_advection!(
                            #     num, grid_p, grid_p.LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
                            #     BC_uL, BC_vL, op.opL
                            # )
                        end

                    elseif num.pressure_velocity_coupling > 1

                        if ns_liquid_phase
                            geoL = [grid_p.LS[iLS].geoL for iLS in 1:num._nLS]
                            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                            # Mum1_L is put in B matrix that multiplies v

                            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = coupled_pressure_velocity!(
                                time_scheme, BC_int,
                                num, grid_p, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                                BC_uL, BC_vL, BC_pL,
                                op.opC_pL, op.opC_uL, op.opC_vL, op.opL,
                                AuL, BuL, AvL, BvL, AϕL, AuvL, BuvL,rhs_uv,
                                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L,
                                Cum1L, Cvm1L, Mum1_L, Mvm1_L,
                                periodic_x, periodic_y, ns_advection, advection, num.current_iter, Ra, navier,pres_free_surfaceL,jump_mass_transfer_rateL,mass_transfer_rate
                            )
                            # if num.current_iter == 1
                            #     phL.u .= -0.5 .* grid_u.y .+ getproperty.(grid_u.LS[1].geoL.centroid, :y) .* grid_u.dy
                            #     phL.v .= 0.5 .* grid_v.x .+ getproperty.(grid_v.LS[1].geoL.centroid, :x) .* grid_v.dx
                            #     phL.u[grid_u.LS[1].SOLID] .= 0.0
                            #     phL.v[grid_v.LS[1].SOLID] .= 0.0
                            # end
                            # linear_advection!(
                            #     num, grid_p, grid_p.LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
                            #     BC_uL, BC_vL, op.opL
                            # )
                        end

                    end #if num.pressure_velocity_coupling

                end # if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1)

            end # if navier_stokes


            #region conservation, divergence checks
            # TODO check global mass conservation and divergence free
            not_divergence_free = true
            #need one-fluid capa
            conservation = compute_conservation_mass(num,phL, grid_p ,grid_u, grid_v, rho_one_fluid)

            # Compute divergence of velocity
            Duv = op.opC_pL.AxT * vec1(phL.uD,grid_u) .+ op.opC_pL.Gx_b * vecb(phL.uD,grid_u) .+
            op.opC_pL.AyT * vec1(phL.vD,grid_v) .+ op.opC_pL.Gy_b * vecb(phL.vD,grid_v)
            for iLS in 1:num.nLS
                if !is_navier(BC_int[iLS]) && !is_navier_cl(BC_int[iLS]) #otherwise normal velocity null if no blowing
                    Duv .+= op.opC_pL.Gx[iLS] * veci(phL.uD,grid_u,iLS+1) .+ 
                            op.opC_pL.Gy[iLS] * veci(phL.vD,grid_v,iLS+1)
                end
            end
            #TODO divergence when Navier ?

            PDI_status = @ccall "libpdi".PDI_multi_expose("print_conservation"::Cstring,
            "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
            "conservation"::Cstring, conservation::Ref{Cdouble}, PDI_OUT::Cint,
            "velocity_divergence"::Cstring, Duv::Ptr{Cdouble}, PDI_OUT::Cint,
            # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

            if maximum(Duv)< num.epsilon_divergence
                not_divergence_free = false
            end

            if abs(conservation) > num.epsilon_conservation || not_divergence_free
                
                PDI_status = @ccall "libpdi".PDI_multi_expose("conservation_error"::Cstring,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint
            end
            #endregion conservation, divergence checks

            # PDI_status = @ccall "libpdi".PDI_multi_expose("check_pressure_velocity"::Cstring,
            # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
            # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
            # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
            # C_NULL::Ptr{Cvoid})::Cint

    end
        #endregion Navier-Stokes 


        #region old block phase change
        #cf old_phase_change_block.jl 
        #endregion old block phase change

        # printstyled(color=:red, @sprintf "\n after phase change radius: %.2e \n" num.current_radius)

        if verbose && adaptative_t
            println("num.timestep_n = $num.timestep_n")
        end

        printstyled(color=:red, @sprintf "\n advection")
        print("\n num.advection_LS_mode ",num.advection_LS_mode," advection ",advection)
        printstyled(color=:red, @sprintf "\n advection")


        #region Advection 
        if advection || electrolysis_advection

            if num.io_pdi>0

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_before_LS_adv"::Cstring,                
                "normal_velocity_intfc"::Cstring, grid_p.V::Ptr{Cdouble}, PDI_OUT::Cint,                   
                C_NULL::Ptr{Cvoid})::Cint

            end # if num.io_pdi>0

            printstyled(color=:red, @sprintf "\n select advection")


            select_advection!(num, grid_p, BC_int, BC_u, grid_u, grid_v, CFL_sc, periodic_x, periodic_y, 
                θ_out, rhs_LS, utmp, electrolysis_phase_change_case, mass_transfer_rate, levelset_1D,
                volume_fraction,
                tmp_vec_p,tmp_vec_p0,
                tmp_vec_u,tmp_vec_v,tmp_vec_u0,tmp_vec_v0,tmp_vec_u1,tmp_vec_v1,
                op,
                phL, u_ext_vel, v_ext_vel)



            # printstyled(color=:red, @sprintf "\n after advection_LS_mode radius: %.2e \n" num.current_radius)

            #region reinitialize Levelset

            #TODO document
            if num.levelset_reinitialize == 0


                if analytical
                    u[grid_p.ind.b_top[1]] .= sqrt.(grid_p.x[grid_p.ind.b_top[1]] .^ 2 + grid_p.y[grid_p.ind.b_top[1]] .^ 2) .- (num.R + speed*num.current_iter*num.timestep_n);
                    u[grid_p.ind.b_bottom[1]] .= sqrt.(grid_p.x[grid_p.ind.b_bottom[1]] .^ 2 + grid_p.y[grid_p.ind.b_bottom[1]] .^ 2) .- (num.R + speed*num.current_iter*num.timestep_n);
                    u[grid_p.ind.b_left[1]] .= sqrt.(grid_p.x[grid_p.ind.b_left[1]] .^ 2 + grid_p.y[grid_p.ind.b_left[1]] .^ 2) .- (num.R + speed*num.current_iter*num.timestep_n);
                    u[grid_p.ind.b_right[1]] .= sqrt.(grid_p.x[grid_p.ind.b_right[1]] .^ 2 + grid_p.y[grid_p.ind.b_right[1]] .^ 2) .- (num.R + speed*num.current_iter*num.timestep_n);
                elseif num.nb_reinit > 0
                    if auto_reinit == 1 && (num.current_iter-1)%num.reinit_every == 0
                        for iLS in 1:num.nLS
                            if !is_wall(BC_int[iLS])
                                ls_rg, rl_rg_v = rg(num, grid_p, grid_p.LS[iLS].u, periodic_x, periodic_y, BC_int)
                                println("$(ls_rg)")
                                printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
                                if ls_rg >= num.δreinit || num.current_iter == 1
                                    print("(ls_rg >= num.δreinit || num.current_iter == 1): yes")
                                    # println("yes")
                                    RK2_reinit!(ls_scheme, grid_p, grid_p.ind, iLS, grid_p.LS[iLS].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                                    
                                    ls_rg, rl_rg_v = rg(num, grid_p, grid_p.LS[iLS].u, periodic_x, periodic_y, BC_int)
                                    println("$(ls_rg) ")
                                    printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
                                end
                            end
                        end
                    elseif (num.current_iter-1)%num.reinit_every == 0
                        for iLS in 1:num.nLS
                            if !is_wall(BC_int[iLS])
                                RK2_reinit!(ls_scheme, grid_p, grid_p.ind, iLS, grid_p.LS[iLS].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                            end
                        end
                    # elseif num.nLS > 1
                    #     for iLS in 1:num.nLS
                    #         if !is_wall(BC_int[iLS])
                    #             RK2_reinit!(ls_scheme, grid_p, grid_p.ind, iLS, grid_p.LS[iLS].u, 2num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int, true)
                    #         end
                    #     end
                            end
                        end

                # Numerical breakup
                if free_surface && breakup ==1
                    count, id_break = breakup_n(grid_p.LS[1].u, grid_p.nx, grid_p.ny, grid_p.dx, grid_p.dy, periodic_x, periodic_y, NB_indices, 5e-2)
                    println(count)
                    if count > count_limit_breakup
                        println("BREAK UP!!") 
                        breakup_f(grid_p, grid_p.LS[1].u, id_break)
                        RK2_reinit!(ls_scheme, grid_p, grid_p.ind, 1, grid_p.LS[1].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                    end
                end

            elseif num.levelset_reinitialize == -1
                printstyled(color=:red, @sprintf "\n No reinit")


            end #if num.levelset_reinitialize

        end

        #endregion reinitialize Levelset
            
        # printstyled(color=:red, @sprintf "\n after reinit radius: %.2e \n" num.current_radius)

        if verbose
            if (num.current_iter-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) | t = %.2e \n" (num.current_iter) 100*(num.current_iter)/num.max_iterations num.time)
                
                #TODO CFL
                # printstyled(color=:green, @sprintf "\n num.CFL : %.2e num.CFL : %.2e num.timestep_n : %.2e\n" num.CFL max(abs.(grid_p.V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)*num.timestep_n/num.Δ num.timestep_n)
                
                if heat && length(grid_p.LS[end].MIXED) != 0
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" mean(grid_p.V[grid_p.LS[1].MIXED]) findmax(grid_p.V[grid_p.LS[1].MIXED])[1] findmin(grid_p.V[grid_p.LS[1].MIXED])[1])
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(grid_p.LS[1].κ[grid_p.LS[1].MIXED]) findmax(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] findmin(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1])
                elseif advection && length(grid_p.LS[end].MIXED) != 0
                    V_mean = mean([mean(grid_u.V[grid_p.LS[1].MIXED]), mean(grid_v.V[grid_p.LS[1].MIXED])])
                    V_max = max(findmax(grid_u.V[grid_p.LS[1].MIXED])[1], findmax(grid_v.V[grid_p.LS[1].MIXED])[1])
                    V_min = min(findmin(grid_u.V[grid_p.LS[1].MIXED])[1], findmin(grid_v.V[grid_p.LS[1].MIXED])[1])
                    # print(@sprintf "Vol_ratio = %.3f%%\n" (volume(grid_p.LS[end].geoL) / V0L * 100))
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(grid_p.LS[1].κ[grid_p.LS[1].MIXED]) findmax(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] findmin(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1])
                end
                if navier_stokes
                    if ns_solid_phase
                        normuS = norm(phS.u)
                        normvS = norm(phS.v)
                        normpS = norm(phS.p.*num.timestep_n)
                        print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                    end
                    if ns_liquid_phase
                        # normuL = norm(phL.u)
                        # normvL = norm(phL.v)
                        # normpL = norm(phL.p.*num.timestep_n)
                        # print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                        if electrolysis
                            # print_electrolysis_statistics(num,grid_p,phL) 
                            PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                            "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                            "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                            "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                            "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                            "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                            
                            C_NULL::Ptr{Cvoid})::Cint
                        end 
                    end
                end
            end
        end

        PDI_status = @ccall "libpdi".PDI_multi_expose("write_after_advection"::Cstring,
        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,      
        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,      
        C_NULL::Ptr{Cvoid})::Cint

        # if levelset && (advection || num.current_iter<2 || electrolysis_advection)
        if levelset && (advection || num.current_iter<2)
            try
                NB_indices = update_all_ls_data(num, grid_p, grid_u, grid_v, BC_int, periodic_x, periodic_y) 
            catch errorLS
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                printstyled(color=:red, @sprintf "\n grid_p.LS not updated \n")
                print(errorLS)
                print(errorLS.task.exception)
                return
            end
            # printstyled(color=:red, @sprintf "\n levelset 4:\n")
            # println(grid_p.LS[1].geoL.dcap[1,1,:])

            grid_p.LS[end].geoL.fresh .= false
            grid_u.LS[end].geoL.fresh .= false
            grid_v.LS[end].geoL.fresh .= false

            get_fresh_cells!(grid_p, grid_p.LS[end].geoL, Mm1_L, grid_p.ind.all_indices)
            get_fresh_cells!(grid_u, grid_u.LS[end].geoL, Mum1_L, grid_u.ind.all_indices)
            get_fresh_cells!(grid_v, grid_v.LS[end].geoL, Mvm1_L, grid_v.ind.all_indices)

            FRESH_L_u = findall(grid_u.LS[end].geoL.fresh)
            FRESH_L_v = findall(grid_v.LS[end].geoL.fresh)

            if navier_stokes

                init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1),
                    grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1),
                    grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                if num.nLS>1 #not monophasic
                    init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1),
                        grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                    init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1),
                        grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                end
            end

            if num.solve_solid == 1 
                grid_p.LS[end].geoS.fresh .= false
                grid_u.LS[end].geoS.fresh .= false
                grid_v.LS[end].geoS.fresh .= false

                get_fresh_cells!(grid_p, grid_p.LS[end].geoS, Mm1_S, grid_p.ind.all_indices)
                get_fresh_cells!(grid_u, grid_u.LS[end].geoS, Mum1_S, grid_u.ind.all_indices)
                get_fresh_cells!(grid_v, grid_v.LS[end].geoS, Mvm1_S, grid_v.ind.all_indices)

                FRESH_S_u = findall(grid_u.LS[end].geoS.fresh)
                FRESH_S_v = findall(grid_v.LS[end].geoS.fresh)
                
                if navier_stokes
                    init_fresh_cells!(grid_u, veci(phL.uD,grid_u,1), veci(phL.uD,grid_u,1),
                        grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                    init_fresh_cells!(grid_v, veci(phL.vD,grid_v,1), veci(phL.vD,grid_v,1),
                        grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                    if num.nLS>1 #not monophasic
                        init_fresh_cells!(grid_u, veci(phL.uD,grid_u,2), veci(phL.uD,grid_u,1),
                            grid_u.LS[end].geoL.projection, FRESH_L_u, periodic_x, periodic_y)
                        init_fresh_cells!(grid_v, veci(phL.vD,grid_v,2), veci(phL.vD,grid_v,1),
                            grid_v.LS[end].geoL.projection, FRESH_L_v, periodic_x, periodic_y)
                    end
                end
            end

            if iszero(num.current_iter%num.save_every) || num.current_iter==num.max_iterations
                snap = num.current_iter÷num.save_every+1
                if save_radius
                    radius[snap] = find_radius(grid_p, grid_p.LS[1])
                end
                if hill
                    a = zeros(length(grid_p.LS[1].MIXED))
                    for i in eachindex(grid_p.LS[1].MIXED)
                        a[i] = grid_p.LS[1].geoL.projection[grid_p.LS[1].MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                # if save_length
                #     fwd.length[snap] = arc_length2(grid_p.LS[1].geoS.projection, grid_p.LS[1].MIXED)
                # end
            end
        end


        if num.nLS>1 #not monophasic
            intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num, intfc_seg_num = convert_interfacial_D_to_segments(num,grid_p,phL.TD,1,2)
    
            barycenter_x_coord = mean(intfc_vtx_x)

            PDI_status = @ccall "libpdi".PDI_multi_expose("update_levelset"::Cstring,
            "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,      
            "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,        
            "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
            "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
            "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
            "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
            "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
            "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
            "barycenter_x_coord"::Cstring, barycenter_x_coord::Ref{Cdouble}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint
        end

        #endregion Advection 

        #region old Navier-Stokes block
        # cf old_Navier_Stokes_block.jl
        #endregion old Navier-Stokes block


        #region end loop post-processing

        if (num.one_fluid_model == 1 &&  num.solve_Navier_Stokes_liquid_phase == 1) 

            # interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)

            tmp_vec_p  .= 0.0
            tmp_vec_p0 .= 0.0

            for j = 1:grid_p.ny
            for i = 1:grid_p.nx
                tmp_vec_p[j,i] =(phL.u[j,i]+phL.u[j,i+1])/2
                tmp_vec_p0[j,i]=(phL.v[j,i]+phL.v[j+1,i])/2
            end
            end

            #TODO compute once, write once
            update_one_fluid_density_viscosity(num,grid_p,grid_u,grid_v,volume_fraction,levelset_one_fluid,rho_one_fluid,
                                                rho_one_fluid_u,rho_one_fluid_v,mu_one_fluid,tmp_vec_p0)

        end

        #endregion end loop post-processing

        # cD, cL, D, L = force_coefficients!(num, grid_p, grid_u, grid_v, op.opL, fwd, phL; step = num.current_iter+1, saveCoeffs = false)

        # if iszero(num.current_iter%num.save_every) || num.current_iter==num.max_iterations
        #     snap = num.current_iter÷num.save_every+1
        #     if num.current_iter==num.max_iterations
        #         snap = size(fwd.T,1)
        #     end
        #     fwd.t[snap] = num.time
        #     @views fwd.V[snap,:,:] .= grid_p.V
        #     if advection
        #         fwdS.Vratio[snap] = volume(grid_p.LS[end].geoS) / V0S
        #         fwdL.Vratio[snap] = volume(grid_p.LS[end].geoL) / V0L
        #     end
        # end
        # @views fwd.Cd[num.current_iter+1] = cD
        # @views fwd.Cl[num.current_iter+1] = cL
        # # @views fwd.radius[num.current_iter+1] = num.current_radius

        # PDI (IO)
        if electrolysis
            if num.io_pdi>0

                try
                    # printstyled(color=:red, @sprintf "\n PDI test \n" )
            
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    # Compute electrical current, interpolate velocity on scalar grid_p
                    #                     if num.electrical_potential>0 compute_grad_phi_ele!(num, grid_p, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current

                    if num.electrical_potential>0
                        compute_grad_phi_ele!(num, grid_p, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL,
                        op.opC_pL, elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                    end

                    # #store in us, vs instead of Eus, Evs
                    # interpolate_grid_liquid!(grid_p,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

                    # #TODO i_current_mag need cond

                    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # C_NULL::Ptr{Cvoid})::Cint

                    # interpolate_grid_liquid!(grid_p,grid_u,grid_v,phL.Eu, phL.Ev,Eus,Evs)
                    
                    interpolate_staggered_u_v_to_scalar_grid_one_fluid_or_one_phase!(num,grid_p,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
                        
                    # num.iLSpdi = 1 # TODO all grid_p.LS

                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                    "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                    "timestep"::Cstring, num.timestep_n::Ref{Cdouble}, PDI_OUT::Cint,  
                    "nx"::Cstring, grid_p.nx::Ref{Clonglong}, PDI_OUT::Cint,
                    "ny"::Cstring, grid_p.ny::Ref{Clonglong}, PDI_OUT::Cint,    
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,                   
                    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                    C_NULL::Ptr{Cvoid})::Cint
            
                        
            
                catch error
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                    print(error)
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                end

               
            end #if io_pdi



            if num.status == 1 #due to num.nH2<0...
                return
            end


            #region test NaN
            if num.solve_solid == 1

                if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                    any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                    any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                    norm(phL.u) > 1e8 || norm(phL.T) > 1e8 || 
                    norm(phS.u) > 1e8 || norm(phS.T) > 1e8 || 
                    # norm(phL.trans_scal) > 1e8 || 
                    norm(phL.phi_ele) > 1e8 
                    # ||
                    # any(phL.trans_scal .<0)
                    )
                    println(@sprintf "\n CRASHED start \n")
                    
                    print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) ,
                    "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                    "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                    "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , 
                    "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 ,
                    "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8,"\n any(phL.trans_scal .<0): ", any(phL.trans_scal .<0))
    
                    num.status =  1
    
                end
            else 
    
                if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                    any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                    norm(phL.u) > 1e8 || norm(phL.T) > 1e8 || 
                    # norm(phL.trans_scal) > 1e8 || 
                    norm(phL.phi_ele) > 1e8 
                    # ||
                    # any(phL.trans_scal .<0)
                    )
                    println(@sprintf "\n CRASHED start \n")
    
                    # println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                    
                    print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) ,
                    "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                    "\n phL.u: ",norm(phL.u) > 1e8, "\n phL.T: ",norm(phL.T) > 1e8 , 
                    "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 ,
                    "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8,"\n any(phL.trans_scal .<0): ", any(phL.trans_scal .<0))
    
                    num.status =  1
    
                end
    
            end

            if num.status == 1
                
                PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                           
                        C_NULL::Ptr{Cvoid})::Cint
               
                return num.current_iter
            end


          
        else #not electrolysis

            if num.solve_solid == 1

                if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                    any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                    norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
                    println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                
                    num.status =  1
                    return
                    
                end

            else

                if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || 
                    norm(phL.u) > 1e8 || norm(phL.T) > 1e8 )
                    println(@sprintf "\n CRASHED after %d iterations \n" num.current_iter)
                
                    num.status =  1
                    return
                    
                end

            end

        end #if electrolysis
        
        #endregion test NaN


        #region update iter number and time

        # num.current_iter += 1
      
        if num.time + num.timestep_n > num.end_time
            print("\n num.time + num.timestep_n > num.end_time, break")
            break
        end

       
        #endregion update iter number and time

    end 
    #endregion time loop

    #region print end
    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) | t = %.2e \n" (num.current_iter-1) 100*(num.current_iter-1)/num.max_iterations num.time)
            print("\n num.time ",num.time," stop_simulation ",num.stop_simulation)
            print("\n num.time ",num.time," num.end_time",num.end_time,"num.timestep_n",num.timestep_n)
            print("\n")
            if stefan && advection
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(grid_p.V[grid_p.LS[1].MIXED]) findmax(grid_p.V[grid_p.LS[1].MIXED])[1] findmin(grid_p.V[grid_p.LS[1].MIXED])[1] std(grid_p.V[grid_p.LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(grid_p.LS[1].κ[grid_p.LS[1].MIXED]) findmax(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] findmin(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] std(grid_p.LS[1].κ[grid_p.LS[1].MIXED]))
            end
            if free_surface && advection
                # print(@sprintf "Vol_ratio = %.3f%%\n" (volume(grid_p.LS[end].geoL) / V0L * 100))
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(grid_p.V[grid_p.LS[1].MIXED]) findmax(grid_p.V[grid_p.LS[1].MIXED])[1] findmin(grid_p.V[grid_p.LS[1].MIXED])[1] std(grid_p.V[grid_p.LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(grid_p.LS[1].κ[grid_p.LS[1].MIXED]) findmax(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] findmin(grid_p.LS[1].κ[grid_p.LS[1].MIXED])[1] std(grid_p.LS[1].κ[grid_p.LS[1].MIXED]))
            end
            if navier_stokes
                if ns_solid_phase
                    normuS = norm(phS.u)
                    normvS = norm(phS.v)
                    normpS = norm(phS.p.*num.timestep_n)
                    print("$(@sprintf("norm(uS) %.6e", normuS))\t$(@sprintf("norm(vS) %.6e", normvS))\t$(@sprintf("norm(pS) %.6e", normpS))\n")
                end
                if ns_liquid_phase
                    normuL = norm(phL.u)
                    normvL = norm(phL.v)
                    normpL = norm(phL.p.*num.timestep_n)
                    print("$(@sprintf("norm(uL) %.6e", normuL))\t$(@sprintf("norm(vL) %.6e", normvL))\t$(@sprintf("norm(pL) %.6e", normpL))\n")
                    if electrolysis
                        print_electrolysis_statistics(num,grid_p,phL)
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, num.time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid_p.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[num.iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,                        
                        C_NULL::Ptr{Cvoid})::Cint

                    end 
                end
            end
            print("\n\n")
        catch
            @show (length(grid_p.LS[end].MIXED))
        end
    end #if verbose
    #endregion print end

    if levelset && (save_radius || hill)
        #TODO save radius
        return # radius
    # elseif flapping
    #     return xc, yc
    else
        return
    end
end

