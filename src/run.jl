
#TODO test with run_forward! or alloc grid inside or ; or no unpack
#TODO remove duplicate variables: cap (already dcap), emptied (already iso),...
#TODO floating-point comparison
#TODO Crank-Nicolson with moving interfaces: simulations explode
#TODO flow rate, check conservation
#TODO i_current_x ... multiply by elec_cond
#TODO no concentration prefactor


# About PDI: the the Parallel Data Interface https://pdi.dev/1.8/

"""

Main function of Flower.jl code to run a simulation

"""
function run_forward!(
    num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    op::DiscreteOperators{Float64, Int64}, 
    phS::Phase{Float64}, 
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
    electrolysis_reaction::String = "nothing",
    imposed_velocity::String = "none",
    adapt_timestep_mode::Int64 = 0,
    non_dimensionalize::Int64=1,
    mode_2d::Int64=0,
    test_laplacian::Bool = false,
    )


    # Initialize simulation parameters

    # λ::Float64 = 1,
    λ = 1 #for Stefan velocity
    # speed::Float64 = 0.0,
    speed = 0.0

    #index of bubble interface
    iLSbubble = 1
    
    #TODO Re
    #TODO homogenize set_poisson... functions which use -b - in the original implementation and not +b + for the Robin BC


    if num.epsilon_mode == 1 || num.epsilon_mode ==2
        num.epsilon_dist = eps(0.01) * num.Δ
        num.epsilon_vol = (eps(0.01)*num.Δ)^2
        #TODO kill dead cells
        #TODO 1e-...
    end


    if length(BC_int) != num.nLS
        @error ("You have to specify $(num.nLS) boundary conditions.")
        return nothing
    end

    crashed=false
    num.current_i=1
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

    if electrolysis
        electrolysis_advection = true
   
    end

    if free_surface && stefan
        @error ("Cannot advect the levelset using both free-surface and stefan condition.")
        return nothing
    elseif free_surface || stefan #|| electrolysis_phase_change_case !="none"
        advection = true
    else
        advection = false
    end

    # The count threshold shouldn't be smaller than 2
    count_limit_breakup = 6
    
    if num.verbosity > 0
        printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e\n" num.CFL num.τ)
    end

    if adapt_timestep_mode !=0
        num.τ = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
        # print("after adapt_timestep!")
        printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e num.τ : %.2e\n" num.CFL num.τ num.τ)
    end

    pres_free_surfaceS = 0.0
    # pres_free_surfaceL = 0.0
    pres_free_surfaceL = num.pres0

    if occursin("levelset",electrolysis_phase_change_case)
        jump_mass_fluxS = false 
        jump_mass_fluxL = true
    else
        jump_mass_fluxS = false 
        jump_mass_fluxL = false
    end

    mass_fluxS = 0.0
    mass_fluxL = 0.0
    

    iRe = 1.0 / num.Re

    # CFL_sc is not a CFL, it is supposed to be the timestep divided by the cell volume
    CFL_sc = num.τ / num.Δ^2

    irho = 1.0

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
    if electrolysis
        
        if electrolysis_phase_change_case != "none"
            num.current_radius = num.R

            printstyled(color=:green, @sprintf "\n radius: %.2e \n" num.current_radius)

            p_liq= num.pres0 + mean(veci(phL.pD,grid,2)) #TODO here one bubble
            p_g=p_liq + 2 * num.σ / num.current_radius

            #TODO using num.temperature0
            if mode_2d==0
                nH2 = p_g * 4.0 / 3.0 * pi * num.current_radius ^ 3 / (num.temperature0 * num.Ru) 
            elseif mode_2d == 1 #reference thickness for a cylinder
                nH2 = p_g * pi * num.current_radius ^ 2 * num.ref_thickness_2d / (num.temperature0 * num.Ru) 
            elseif mode_2d==2 #mol/meter
                nH2=num.concentration0[num.index_phase_change]* pi * num.current_radius ^ 2
            elseif mode_2d==3 #mol/meter half circle
                nH2=1.0/2.0*num.concentration0[num.index_phase_change]* pi * num.current_radius ^ 2
            end
            # nH2 = 4.0/3.0 * pi * num.current_radius^3 * num.rho2 / num.MWH2

            printstyled(color=:green, @sprintf "\n Mole: %.2e \n" nH2)

            printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" num.concentration0[num.index_phase_change]*4.0/3.0*pi*num.current_radius^3 p_g*4.0/3.0*pi*num.current_radius^3/(num.temperature0*num.Ru))

        end

    end # if electrolysis    
    

    local NB_indices;


    #region Allocations
    # Allocations

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
    utmp = copy(grid.LS[1].u)
    rhs_LS = fzeros(grid)

    #vectors used reset at start of functions to limit allocations
    tmp_vec_u = zeros(grid_u) 
    tmp_vec_v = zeros(grid_v) 
    tmp_vec_u0 = zeros(grid_u) 
    tmp_vec_v0 = zeros(grid_v)
    tmp_vec_p = zeros(grid) 
    tmp_vec_p0 = zeros(grid) 
    tmp_vec_p1 = zeros(grid) 

    tmp_vec_1D = fnzeros(grid,num)
    tmp_vec_1D_2 = fnzeros(grid,num)

    # tmp_vec_1D_u = fnzeros(grid,num)
    tmp_vec_1D_v = fnzeros(grid_v,num)
    tmp_vec_1D_v0 = fnzeros(grid_v,num)


    if electrolysis
        if num.nb_transported_scalars>1
            elec_cond = zeros(grid)
            elec_condD = fnzeros(grid,num)

            if heat 
                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
                elec_cond .= reshape(vec1(elec_condD,grid),grid)
                # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
                # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
            else
                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
                elec_cond .= reshape(vec1(elec_condD,grid),grid)
                # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
                # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
            end
        else 
            elec_cond = ones(grid)
            elec_condD = fnones(grid,num)
            printstyled(color=:green, @sprintf "\n conductivity one")
        end 

        i_butler = zeros(grid.ny) #left wall
    end #electrolysis

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
            n_snaps = iszero(num.max_iterations%num.save_every) ? num.max_iterations÷num.save_every+1 : num.max_iterations÷num.save_every+2
            local radius = zeros(n_snaps)
            radius[1] = find_radius(grid, grid.LS[1])
        end
        if hill
            local radius = zeros(num.max_iterations+1)
            a = zeros(length(grid.LS[1].MIXED))
            for i in eachindex(grid.LS[1].MIXED)
                a[i] = grid.LS[1].geoL.projection[grid.LS[1].MIXED[i]].pos.y
            end
            radius[1] = mean(a)
        end
    elseif !levelset
        grid.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_u.LS[1].MIXED = [CartesianIndex(-1,-1)]
        grid_v.LS[1].MIXED = [CartesianIndex(-1,-1)]
    end

    # if save_length
    #     fwd.length[1] = arc_length2(grid.LS[1].geoS.projection, grid.LS[1].MIXED)
    # end

    #endregion


    #region Initialisation of bulk and interfacial values
    # Initialisation of bulk and interfacial values
    
    #TODO which grid.LS
    presintfc = 0.0
    # presintfc = pres0 + p_lapl ? #TODO init pressure
    #TODO perio, intfc, ... check init_fields_2!

    #No scalar in "solid" phase
    # @views init_fields_multiple_levelsets!(num,phS.trans_scalD[:,iscal],phS.trans_scal[:,:,iscal],HS,BC_trans_scal[iscal],grid,num.concentration0[iscal])
    # @views phS.trans_scal[:,:,iscal] .= num.concentration0[iscal]

    # TODO reset zero
    tmp_vec_u .= 0.0

    #tmp_vec_u
    #from LS 1 to centroid grid_u.LS[end].geoS
    get_height!(grid_u.LS[1],grid_u.ind,grid_u.dx,grid_u.dy,grid_u.LS[end].geoS,tmp_vec_u) #here tmp_vec_u solid

    init_fields_multiple_levelsets!(num,phS.uD,phS.u,tmp_vec_u,BC_uS,grid_u,num.uD,"uS")
    # init_fields_multiple_levelsets!(num,phS.ucorrD,phS.u,HSu,BC_uS,grid_u,num.uD)

    get_height!(grid_u.LS[1],grid_u.ind,grid_u.dx,grid_u.dy,grid_u.LS[end].geoL,tmp_vec_u)  #here tmp_vec_u liquid

    init_fields_multiple_levelsets!(num,phL.uD,phL.u,tmp_vec_u,BC_uL,grid_u,num.uD,"uL")
    # init_fields_multiple_levelsets!(num,phL.ucorrD,phL.u,HLu,BC_uL,grid_u,num.uD)


    # TODO reset zero
    tmp_vec_v .= 0.0

    get_height!(grid_v.LS[1],grid_v.ind,grid_v.dx,grid_v.dy,grid_v.LS[end].geoS,tmp_vec_v) 

    init_fields_multiple_levelsets!(num,phS.vD,phS.v,tmp_vec_v,BC_vS,grid_v,num.vD,"vS")
    # init_fields_multiple_levelsets!(num,phS.vcorrD,phS.v,HSv,BC_vS,grid_v,num.vD)

    get_height!(grid_v.LS[1],grid_v.ind,grid_v.dx,grid_v.dy,grid_v.LS[end].geoL,tmp_vec_v)

    init_fields_multiple_levelsets!(num,phL.vD,phL.v,tmp_vec_v,BC_vL,grid_v,num.vD,"vL")
    # init_fields_multiple_levelsets!(num,phL.vcorrD,phL.v,HLv,BC_vL,grid_v,num.vD)


    # TODO reset zero
    tmp_vec_p .= 0.0

    get_height!(grid.LS[1],grid.ind,grid.dx,grid.dy,grid.LS[end].geoS,tmp_vec_p) #here tmp_vec_p solid

    init_fields_multiple_levelsets!(num,phS.pD,phS.p,tmp_vec_p,BC_pS,grid,presintfc,"pS")

    if heat
        init_fields_multiple_levelsets!(num,phS.TD,phS.T,tmp_vec_p,BC_TS,grid,num.θd,"TS")
    end

    #Electrolysis
    if electrolysis && num.electrical_potential > 0
        init_fields_multiple_levelsets!(num,phS.phi_eleD,phS.phi_ele,tmp_vec_p,BC_phi_ele,grid,num.phi_ele0,"phiS")
    end

    get_height!(grid.LS[1],grid.ind,grid.dx,grid.dy,grid.LS[end].geoL,tmp_vec_p) #here tmp_vec_p liquid

    init_fields_multiple_levelsets!(num,phL.pD,phL.p,tmp_vec_p,BC_pL,grid,presintfc,"pL")

    if heat
        init_fields_multiple_levelsets!(num,phL.TD,phL.T,tmp_vec_p,BC_TL,grid,num.θd,"TL")
    end

    # Electrolysis
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check %s %s %s %s %.2e %.2e %2i\n" heat heat_convection electrolysis electrolysis_convection num.τ num.θd num.nb_transported_scalars)

        for iscal=1:num.nb_transported_scalars
            @views phL.trans_scal[:,:,iscal] .= num.concentration0[iscal]
            @views init_fields_multiple_levelsets!(num,phL.trans_scalD[:,iscal],phL.trans_scal[:,:,iscal],
            tmp_vec_p,BC_trans_scal[iscal],grid,num.concentration0[iscal],"scalL")
            # tmp_vec_p,BC_trans_scal[iscal],grid,num.concentration0[iscal],"testscalL")

        end

        if num.electrical_potential > 0 #TODO init phi =0 or with Neumann for BC of concentration? for now not done since "phiL" given in arg
            init_fields_multiple_levelsets!(num,phL.phi_eleD,phL.phi_ele,tmp_vec_p,BC_phi_ele,grid,num.phi_ele0,"phiL")
            if num.electrical_potential == 3
                vecb(phL.phi_eleD,grid) .= 0.0
            end
        end
    end  
    
    if electrolysis

        # current_t = 0.
        # num.current_i =0
        nstep = 0
        time = 0.
        iLSpdi = 1


        printstyled(color=:green, @sprintf "\n Check init_fields_2!\n")
        # print_electrolysis_statistics(num,grid,phL)

        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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

        interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)

        PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
            "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
            "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
            "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
            # "normal_x"::Cstring, normal_x::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "normal_y"::Cstring, normal_y::Ptr{Cdouble}, PDI_OUT::Cint,  
            # grid_u.LS[iLS].α
            # "normal_angle"::Cstring, grid.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
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
    
    kill_dead_cells!(phS.T, grid, grid.LS[end].geoS)
    kill_dead_cells!(phL.T, grid, grid.LS[end].geoL)

    #TODO check timestep coefficients num.n-1 

    #Electrolysis
    # TODO kill_dead_cells! for [:,:,iscal]
    if electrolysis
        for iscal=1:num.nb_transported_scalars
            # @views kill_dead_cells!(phS.trans_scal[:,:,iscal], grid, grid.LS[end].geoS) #TODO
            # @views kill_dead_cells!(phL.trans_scal[:,:,iscal], grid, grid.LS[end].geoL) 
            # @views kill_dead_cells_val!(phS.trans_scal[:,:,iscal], grid, grid.LS[end].geoS) #TODO
            # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, grid.LS[end].geoL,num.concentration0[iscal]) 
            @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, grid.LS[end].geoL,0.0) 
            @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

        end
    end #if electrolysis    

    #endregion
    

    # print("\n vecb_L(elec_condD, grid) after kill \n ", vecb_L(phL.trans_scalD[:,2], grid) )


    # if num.nb_transported_scalars>0
    #     printstyled(color=:green, @sprintf "\n after kill \n")
    #     print_electrolysis_statistics(num,phL)
    #     printstyled(color=:green, @sprintf "\n average T %s\n" average!(phL.T, grid, grid.LS[1].geoL, num))
    # end

    # Set matrices/operators
    if is_Forward_Euler(time_scheme) || is_Crank_Nicolson(time_scheme)
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y, false)

        # printstyled(color=:red, @sprintf "\n levelset 2:\n")
        # println(grid.LS[1].geoL.dcap[1,1,:])

        if navier_stokes || heat || electrolysis
            geoS = [grid.LS[iLS].geoS for iLS in 1:num._nLS]
            geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:num._nLS]
            geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:num._nLS]
            Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S = set_matrices!(
                num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS,
                op.opC_pS, op.opC_uS, op.opC_vS, periodic_x, periodic_y
            )

            geoL = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
            geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
            geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]
            Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L = set_matrices!(
                num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL,
                op.opC_pL, op.opC_uL, op.opC_vL, periodic_x, periodic_y
            )
        end

        Mm1_L = copy(op.opC_pL.M)
        Mm1_S = copy(op.opC_pS.M)
        Mum1_L = copy(op.opC_uL.M)
        Mum1_S = copy(op.opC_uS.M)
        Mvm1_L = copy(op.opC_vL.M)
        Mvm1_S = copy(op.opC_vS.M)

        if navier_stokes || heat || electrolysis

            #Allocations
            #TODO pre-allocate at start to save up allocations
            #TODO optimize allocations

            #Preallocate for mass flux computations
            mass_flux_vec1 = fzeros(grid)
            mass_flux_vecb = fzeros(grid)
            mass_flux_veci = fzeros(grid)
            mass_flux = zeros(grid)


            #Allocations for scalar grid
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb

            a1_p = spdiagm(ni,ni,zeros(ni))

            Ascal = spzeros(nt, nt)
            Bscal = spzeros(nt, nt)
            rhs_scal = fnzeros(grid, num)

            all_CUTCT = zeros(grid.ny * grid.nx, num.nb_transported_scalars)

            # us=zeros(grid)
            # vs=zeros(grid)

            #save up memory: do not allocate Eus
            # Eus=zeros(grid) 
            # Evs=zeros(grid)

            AϕS = spzeros(nt, nt)
            AϕL = spzeros(nt, nt)

            if electrolysis 
                # Aphi_eleL = spzeros(nt, nt)

                # coeffDu = zeros(grid_u)
                # coeffDv = zeros(grid_v)
                coeffDx_interface = zeros(grid_u)
                coeffDy_interface = zeros(grid_v)
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
            nt = (num.nLS - num.nNavier + 1) * ni + num.nNavier * grid.nx * grid.ny + nb
            AuvS = spzeros(nt, nt)
            AuvL = spzeros(nt, nt)
            BuvS = spzeros(nt, nt)
            BuvL = spzeros(nt, nt)



            # Adapt cell volume W for gradients 
            # cf 4/3 factor for Laplacian
            if num.laplacian == 1
                compute_divergence!(num, 
                # grid, 
                # grid_u, 
                grid_v, 
                op,
                AvL, 
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
            end  
           
         
            #TODO why this call without interface initialization ?
            if !navier
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
            else
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    op.opC_pS, op.opC_uS, op.opC_vS,
                    AuvS, BuvS,
                    iRe.*Lum1_S, iRe.*bc_Lum1_S, iRe.*bc_Lum1_b_S, Mum1_S, BC_uS,
                    iRe.*Lvm1_S, iRe.*bc_Lvm1_S, iRe.*bc_Lvm1_b_S, Mvm1_S, BC_vS,
                    true
                )
            end
            #TODO remove alloc a0_p
            a0_p = []
            for i in 1:num.nLS
                push!(a0_p, zeros(grid))
            end
            
            
            if !advection
            #call to set_heat! is there to set up the matrices for the heat equation. 
            #If the level-set is not advected, then after this call there is no need to update these matrices anymore

                _ = set_poisson(
                    BC_int, num, grid, a0_p, op.opC_pS, op.opC_uS, op.opC_vS,
                    AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                    true
                )

                set_heat!(
                    BC_int[1], num, grid, op.opC_TS, grid.LS[1].geoS, phS, num.θd, BC_TS, grid.LS[1].MIXED, grid.LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    op.opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, true, BC_int
                )
            end

            if test_laplacian
                num.exact_laplacian = test_laplacian_pressure(num,grid_v,phL.vD,op.opC_pL, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L)
                return
            end

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
                # Coupled resolution if navier BC is activated
                _ = FE_set_momentum_coupled(
                    BC_int, num, grid, grid_u, grid_v,
                    op.opC_pL, op.opC_uL, op.opC_vL,
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
                BC_int, num, grid, a0_p, op.opC_pL, op.opC_uL, op.opC_vL,
                AϕL, Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, BC_pL,
                true
            )

            #call to set_heat! is there to set up the matrices for the heat equation. 
            #If the level-set is not advected, then after this call there is no need to update these matrices anymore
            set_heat!(
                BC_int[1], num, grid, op.opC_TL, grid.LS[1].geoL, phL, num.θd, BC_TL, grid.LS[1].MIXED, grid.LS[1].geoL.projection,
                ATL, BTL,rhs_scal,
                op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                periodic_x, periodic_y, heat_convection, true, BC_int
            )
         
        end
    else
        error("Unknown time scheme. Available options are ForwardEuler and CrankNicolson")
    end

    if heat_convection || electrolysis_convection
        NB_indices = update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
        # println(grid.LS[1].geoL.dcap[1,1,:])
    end

    # V0S = volume(grid.LS[end].geoS)
    # V0L = volume(grid.LS[end].geoL)

    if num.debug =="allocations_start"
        print("\n STOP allocations")
        return
    end


    


    


    current_t = 0.
    num.current_i =1

    #Time loop
    while num.current_i < num.max_iterations + 1

        
        #region Adapt timestep
        # printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e\n" num.CFL num.τ)
        if adapt_timestep_mode !=0
            num.τ = adapt_timestep!(num, phL, phS, grid_u, grid_v,adapt_timestep_mode)
            # print("after adapt_timestep!")
            printstyled(color=:green, @sprintf "\n num.CFL : %.2e dt : %.2e num.τ : %.2e\n" num.CFL num.τ num.τ)
        end
        #endregion

        printstyled(color=:red, @sprintf "\n iter: %5i\n" num.current_i)
        println("\n grid.LS[1].geoL.dcap[1,1,:]",grid.LS[1].geoL.dcap[1,1,:])

        # PDI_multi_expose()

        # PDI_status = @ccall "libpdi".PDI_multi_expose("check_levelset"::Cstring,
        # "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
        # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        # # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
        # # "normal_x"::Cstring, normal_x::Ptr{Cdouble}, PDI_OUT::Cint,   
        # # "normal_y"::Cstring, normal_y::Ptr{Cdouble}, PDI_OUT::Cint,  
        # # grid_u.LS[iLS].α
        # # "normal_angle"::Cstring, grid.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
        # "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
        # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
        # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
        # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
        # C_NULL::Ptr{Cvoid})::Cint

        # PDI (IO)
        if electrolysis

            #TODO not necessary to expose everything now for ex only grid.LS ? and the rest later
    
            printstyled(color=:red, @sprintf "\n test segments\n" )
            # print("\n type of elec_cond ", typeof(elec_cond)," \n")

    
            intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num, intfc_seg_num = convert_interfacial_D_to_segments(num,grid,phL.TD,1,2)
            # print("\n number of interface points intfc_vtx_num ", intfc_vtx_num)
            # print("\n intfc_vtx_connectivities ",intfc_vtx_connectivities)
            # print("\n len ", size(intfc_vtx_connectivities),intfc_seg_num)
    
            # print("\n intfc_vtx_x ",intfc_vtx_x)
            # print("\n intfc_vtx_x ",intfc_vtx_y)
    

            #TODO when to write elec dat, ...
    
            if num.io_pdi>0
    
                try
                    # printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                    time = current_t #Cdouble
                    nstep = num.current_i
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    
                    if num.electrical_potential > 0

                        # print("\n type of elec_cond ", typeof(elec_cond)," \n")
                        try
                            compute_grad_phi_ele!(num, grid, grid_u, grid_v,grid_u.LS[end], grid_v.LS[end], phL, phS, 
                            op.opC_pL, op.opC_pS, elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                        catch
                            @error("compute_grad_phi_ele!")
                            break
                        end
                    end
                    # #store in us, vs instead of Eus, Evs
                    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)
    
                    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # C_NULL::Ptr{Cvoid})::Cint
    
                    interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
            
                    # print("\n before write \n ")
            
                    iLSpdi = 1 # all grid.LS iLS = 1 # or all grid.LS ?
    
                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name
                    
                    printstyled(color=:magenta, @sprintf "\n PDI write_data_start_loop %.5i \n" num.current_i)
    
                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data_start_loop"::Cstring,
                    "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "normal_x"::Cstring, normal_x::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "normal_y"::Cstring, normal_y::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # grid_u.LS[iLS].α
                    # "normal_angle"::Cstring, grid.LS[iLSpdi].α::Ptr{Cdouble}, PDI_OUT::Cint,
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
            grid.V .= speed #*ones(grid.ny, grid.nx)
        end

        #region Heat equation
        # Solve heat equation
        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, grid.LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                set_heat!(
                    BC_int[1], num, grid, op.opC_TS, grid.LS[1].geoS, phS, num.θd, BC_TS, grid.LS[1].MIXED, grid.LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    op.opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs_scal
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, grid.LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                set_heat!(
                    BC_int[1], num, grid, op.opC_TL, grid.LS[1].geoL, phL, num.θd, BC_TL, grid.LS[1].MIXED, grid.LS[1].geoL.projection,
                    ATL, BTL,rhs_scal,
                    op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs_scal
                phL.T .= reshape(veci(phL.TD,grid,1), grid)

            end
        end # #if heat
        #endregion heat equation

        
        #region Electrolysis 
        if electrolysis
            if electrolysis_liquid_phase

                #region Electrolysis: Poisson  
                if num.electrical_potential > 0
                    
                    # Electroneutrality assumption
                    a0_p = [] 
                    for i in 1:num.nLS
                        push!(a0_p, zeros(grid))
                    end

                    # Constant electrical conductivity assumption
                    #TODO electrical conductivity depends on concentration
                    #iKOH index of KOH 
                    # kappa_ele=2*num.Faraday^2*num.concentration0[iKOH]*num.diffusion_coeff[iKOH]/(num.Ru*T)
                    # elec_cond=2*num.Faraday^2*trans_scal[iKOH]*num.diffusion_coeff[iKOH]/(num.Ru*T)

                    # #print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")

                    if electrolysis && num.nb_transported_scalars>1
                        if heat 
                            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
                            elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
                            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
                        else
                            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
                            elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
                            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
                            
                            if num.bulk_conductivity == 3
                                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                                elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            end
                       
                        end
                    else
                        elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                        elec_cond .= reshape(vec1(elec_condD,grid),grid)
                    end

         

                    # printstyled(color=:red, @sprintf "\n test conductivity")

                    # # elec_condD .= 2*num.Faraday^2 .*num.concentration0[2].*num.diffusion_coeff[2]./(num.Ru.*num.temperature0)
                    # test_filter_concentration!(num,grid,phL.trans_scalD[:,2],num.concentration0[2])

                    # elec_condD = 2*num.Faraday^2 .*phL.trans_scalD[:,2].*num.diffusion_coeff[2]./(num.Ru.*num.temperature0)
 
                    # # TODO icurrent mag and replace huge val scal
                    
                    # print("\n test i ",-2*num.Faraday^2*num.concentration0[2]*num.diffusion_coeff[2]./(num.Ru.*num.temperature0))*(num.phi0-num.phi1)/(1e-4)
                    # print("\n test i ",-2*(num.Faraday^2)*num.concentration0[2]*num.diffusion_coeff[2]./(num.Ru.*num.temperature0))*(num.phi0-num.phi1)/(1e-4)

                    # printstyled(color=:red, @sprintf "\n test conductivity")





                    #Update Butler-Volmer Boundary Condition with new potential 

                    if occursin("Butler",electrolysis_reaction) && num.nLS == 1

                        printstyled(color=:red, @sprintf "\n Recomputing Butler \n" )

                        #region Update current
                        if electrolysis_reaction == "Butler_no_concentration"                
                            update_electrical_current_from_Butler_Volmer!(num,grid,heat,phL.phi_eleD,i_butler;phL.T)
                        end
                        #endregion Update current

                        update_BC_electrical_potential!(num,grid,BC_phi_ele,elec_cond,elec_condD,i_butler)


                        # if heat
                        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,phL.T)./elec_cond[:,1]
                        # else
                        #     BC_phi_ele.left.val = -butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./elec_cond[:,1]
                            
                        #     # for iscal=1:num.nb_transported_scalars
                        #     #     BC_trans_scal[iscal].left.val = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./(2*num.Faraday*num.diffusion_coeff[iscal])
                        #     #     if iscal==1 || iscal==2
                        #     #         BC_trans_scal[iscal].left.val .*=-1 #H2O
                        #     #     end
                        #     # end
                        # end    

                    # elseif electrolysis_reaction == ""
                    #     # BC_phi_ele.left.val = -butler_volmer_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,phL.phi_ele[:,1],num.phi_ele1,num.Ru,num.temperature0)./elec_cond
                    
                        
                
                        # TODO 
                        #Remove Nan when dividing by conductivity which may be null

                        # TODO bug 1                           

                        for iLS in 1:num.nLS
                            # kill_dead_bc_left_wall!(vecb(elec_condD,grid), grid, iLS,1.0)
                            for i = 1:grid.ny
                                # print("vecb cap",vecb_L(grid.LS[iLS].geoL.cap[:,5],grid))
                               
                                # II = CartesianIndex(i,1)
                                # II = grid.ind.b_left[1][i]
                                # opC.χ_b[i, i] = geo.dcap[II,1]
                                # TODO not cleat why zero: grid.LS[iLS].geoL.cap[II,1]
                                #TODO cf update LS convection not convection where something is overwritten
                                # wall_liquid_height = grid.LS[iLS].geoL.cap[II,1]
                                wall_liquid_height = op.opC_pL.χ_b[i, i]
                                if wall_liquid_height < 1e-12
                                    BC_phi_ele.left.val[i] = 1.0
                                    print("\n bug BC_phi_ele.left.val[i] ",II," ",grid.LS[iLS].geoL.cap[II,:])
                                    # print("\n opC.χ_b[i, i] ",op.opC_pL.χ_b[i, i])
                                end
                            end
                        end



                    end #if occursin("Butler",electrolysis_reaction)

          
                    #TODO nLS
                    #TODO kill_dead_cells! ?
                    kill_dead_cells!(phL.phi_ele, grid, grid.LS[1].geoL)
                    veci(phL.phi_eleD,grid,1) .= vec(phL.phi_ele)
                    

            

                    # iterate (non-linear BC with Butler) 
                    for poisson_iter=1:num.electrical_potential_max_iter

                        printstyled(color=:orange, @sprintf "\n poisson iter %.2i \n" poisson_iter)

                        compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
                        elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current

                    
                        residual_electrical_potential = maximum(abs.(-tmp_vec_p[div(grid.ny,2),:].+butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                                    num.phi_ele1,num.Ru,num.temperature0)))

                        printstyled(color=:orange, @sprintf "\n Residual %.3e criterion %.3e\n" residual_electrical_potential num.electrical_potential_residual)


                        @ccall "libpdi".PDI_multi_expose("check_electrical_potential"::Cstring,
                        "poisson_iter"::Cstring, poisson_iter ::Ref{Clonglong}, PDI_OUT::Cint,
                        "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                        "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                        "i_current_mag"::Cstring, tmp_vec_p1::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                        "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,  
                        "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "residual_electrical_potential"::Cstring, residual_electrical_potential ::Ref{Cdouble}, PDI_OUT::Cint,
                        # "grad_phi_ele_u"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,  
                        C_NULL::Ptr{Cvoid})::Cint

                        if residual_electrical_potential<num.electrical_potential_residual
                            printstyled(color=:orange, @sprintf "\n End Poisson loop \n")
                            break
                        end
                        
                        # printstyled(color=:orange, @sprintf "\n grad poisson iter %.2i \n" poisson_iter)

                        # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),:]," \n")

                        # @ccall "libpdi".PDI_multi_expose("solve_poisson"::Cstring,
                        # # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                        # # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,  
                        # "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
                        # # "grad_phi_ele_u"::Cstring, tmp_vec_u::Ptr{Cdouble}, PDI_OUT::Cint,  
                        # C_NULL::Ptr{Cvoid})::Cint

                    
                        if electrolysis_reaction == "Butler_no_concentration"

                            # if num.poisson_newton ==1
                            #     vecb_L(phL.phi_eleD, grid) = vecb_L(phL.phi_eleD, grid) - (partial...+)/deriv
                            # end


                            #TODO dev multiple levelsets
                            if heat
                                i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                                num.phi_ele1,num.Ru,phL.T)
                            else
                                if num.nLS == 1
                                    i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                                    num.phi_ele1,num.Ru,num.temperature0)
                                # else
                                    #imposed by LS 2
                                    # iLS_elec = 2
                                    # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid,iLS_elec+1),
                                    # num.phi_ele1,num.Ru,num.temperature0)
                                end
                                    
                            end   

                            if poisson_iter>1
                                if num.bulk_conductivity == 0
                                    BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)
        
                                elseif num.bulk_conductivity == 1
                                    # Recommended as long as cell merging not implemented:
                                    # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
                                    # which produces higher concentration in front of the contact line
                                    BC_phi_ele.left.val .= i_butler./elec_cond[:,1]
        
        
                                elseif num.bulk_conductivity == 2 || num.bulk_conductivity == 3
                                    BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)
        
                                    iLS = 1 #TODO end ? if several grid.LS ?
                                    for j in 1:grid.ny
                                        II = CartesianIndex(j,1)
                                        if grid.LS[iLS].geoL.cap[II,5] < num.ϵ
                                            BC_phi_ele.left.val[j] = i_butler[j]/elec_cond[j,1] 
                                        end
                                    end
                                    
                                # if num.bulk_conductivity == 3
                                #     elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                                #     elec_cond .= reshape(vec1(elec_condD,grid),grid)
                                # end
                                
                                end
                            end 
                            # print("\n BC_phi_ele",BC_phi_ele,"\n")

                            
                        end

                        print("\n BC_phi_ele ",BC_phi_ele)

                        solve_poisson_variable_coeff!(num, 
                        grid, 
                        grid_u, 
                        grid_v, 
                        op.opC_pL,
                        Ascal, 
                        rhs_scal,
                        tmp_vec_p, #a0
                        a1_p,
                        BC_phi_ele,
                        phL,    
                        elec_cond,                    
                        elec_condD,
                        tmp_vec_u,
                        tmp_vec_v,
                        # tmp_vec_u0,
                        # tmp_vec_v0,
                        i_butler,
                        ls_advection,
                        heat)


                        @ccall "libpdi".PDI_multi_expose("solve_poisson"::Cstring,
                        # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                        # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                        "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                        "elec_cond_1D"::Cstring, elec_condD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "BC_phi_ele_left"::Cstring, BC_phi_ele.left.val::Ptr{Cdouble}, PDI_OUT::Cint,  
                        C_NULL::Ptr{Cvoid})::Cint

                        #TODO or linearize 

                        #TODO compute grad

                        if num.electrical_potential>0
                            compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
                            elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                            
                            # printstyled(color=:orange, @sprintf "\n grad poisson iter %.2i \n" poisson_iter)

                            # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),:]," \n")
                            
                            # print("\n grad ", tmp_vec_u[div(grid_u.ny,2),1]," \n")
                            # print("\n BC_phi_ele ", BC_phi_ele.left.val[div(grid_u.ny,2)]," \n")
                            # print("\n i_butler ", i_butler," \n")

                        end

                    end #for loop Poisson

                    # printstyled(color=:cyan, @sprintf "\n after solve_poisson_variable_coeff! \n")
                    # print_electrolysis_statistics(num,grid,phL)

                    PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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

                    if any(isnan, phL.phi_eleD)
                        print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                        "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                        "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
            
                        print("\n phL.phi_eleD: ",any(isnan, phL.phi_eleD),"\n phL.phi_ele: ",any(isnan, phL.phi_ele),"\n")

                        print("\n Ascal: ",any(isnan, Ascal),"\n rhs_scal: ",any(isnan, rhs_scal),"\n")

                        print("\n \n vecb_L",vecb_L(phL.phi_eleD[:,1], grid))

                    end

                
                    # TODO compute magnitude of exchange current
                    # gradient!(::Neumann, Ox, Oy, Bx, By, HNx, HNy, Divx, Divy, dcap, num.n, BC, all_indices, b_left_u, b_bottom_v, b_right_u, b_top_v, b_left_p, b_bottom_p, b_right_p, b_top_p)
                    # TODO add post-treatment variables

                    #TODO update BC concentration

                   

                    if electrolysis && num.nb_transported_scalars>1
                        if heat 
                            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.TD, phL.trans_scalD[:,num.index_electrolyte])
                            elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, phL.T, phL.trans_scal)
                            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
                        else
                            elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scalD[:,num.index_electrolyte])
                            elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            # elec_cond .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, phL.trans_scal)
                            # elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
                       
                            if num.bulk_conductivity == 3
                                elec_condD .= compute_ele_cond.(num.Faraday,num.diffusion_coeff[num.index_electrolyte],num.Ru, num.temperature0, num.concentration0[num.index_electrolyte])
                                elec_cond .= reshape(vec1(elec_condD,grid),grid)
                            end
                        end
                    end

                    # if electrolysis && num.nb_transported_scalars>1
                    #     if heat 
                    #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
                    #     else
                    #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
                    #     end
                    # else 
                    #     elec_cond = ones(grid)
                    #     printstyled(color=:green, @sprintf "\n conductivity one")

                    # end 
                    if num.electrical_potential>0
                        compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS, op.opC_pL, op.opC_pS, 
                        elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                    end

                    # scal_magnitude

                    # phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude

                    printstyled(color=:green, @sprintf "\n test grad")

                    # compute_grad_p!(num,grid, grid_u, grid_v, phL.phi_eleD, op.opC_pL, op.opC_uL, op.opC_vL)


                    # #store in us, vs instead of Eus, Evs
                    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

                    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # C_NULL::Ptr{Cvoid})::Cint

                end #electrical_potential

                #endregion Poisson


               
                #TODO check BC not overwritten by different scalars
                #TODO check ls advection true when num.n scalars

                # print("\n vecb_L",vecb_L(phL.phi_eleD, grid))
                # print("\n phL.phi_ele[:,1]",phL.phi_ele[:,1])

                #TODO phL.phi_ele[:,1] or vecb_L(phL.phi_eleD, grid)
                # we suppose phi(x=0)=... cf Khalighi
                #but here BC
                
                # New start scalar loop

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

                        if ((num.current_i-1)%show_every == 0) 
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
                    # vecb...
                    # ....
                    
                    if ((num.current_i-1)%show_every == 0) 
                        # printstyled(color=:red, @sprintf "\n Imposed velocity v min %.2e max %.2e\n" minimum(phL.vD) maximum(phL.vD))
                        # printstyled(color=:red, @sprintf "\n Imposed velocity u min %.2e max %.2e\n" minimum(phL.uD) maximum(phL.uD))       
                        
                        # printstyled(color=:cyan, @sprintf "\n before scalar transport 0\n")

                        # print_electrolysis_statistics(num,grid,phL)

                        # PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        # "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        # # "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # # "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                        # # "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                        # # "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
                        # # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
                        # # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
                        # # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
                        # # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
                        # # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
                        # # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
                        # C_NULL::Ptr{Cvoid})::Cint

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
                    
                
                end #imposed_velocity
                #endregion Impose velocity

                #region Update current
                if electrolysis_reaction == "Butler_no_concentration"
                    
                    update_electrical_current_from_Butler_Volmer!(num,grid,heat,phL.phi_eleD,i_butler;phL.T)

                    # #TODO dev multiple levelsets
                    # if heat
                    #     i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                    #     num.phi_ele1,num.Ru,phL.T)
                    # else
                    #     if num.nLS == 1
                    #         i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,vecb_L(phL.phi_eleD, grid),
                    #         num.phi_ele1,num.Ru,num.temperature0)
                    #     # else
                    #         #imposed by LS 2
                    #         # iLS_elec = 2
                    #         # i_butler = butler_volmer_no_concentration.(num.alpha_a,num.alpha_c,num.Faraday,num.i0,veci(phL.phi_eleD, grid,iLS_elec+1),
                    #         # num.phi_ele1,num.Ru,num.temperature0)
                    #     end
                    # end   
                end
                #endregion Update current


                #region Scalar transport: update boundary conditions
                if num.nb_transported_scalars>0

                    # printstyled(color=:magenta, @sprintf "\n num.nb_transported_scalars %.5i " num.nb_transported_scalars)

                    for iscal=1:num.nb_transported_scalars
                        @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, grid.LS[1].geoL,0.0) 
                        # @views kill_dead_cells_val!(phL.trans_scal[:,:,iscal], grid, grid.LS[1].geoL,num.concentration0[iscal]) 

                        @views veci(phL.trans_scalD[:,iscal],grid,1) .= vec(phL.trans_scal[:,:,iscal])

                        if electrolysis_reaction == "Butler_no_concentration" && num.nLS == 1

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

                            # for testn in 1:grid.ny
                            #     printstyled(color=:green, @sprintf "\n jtmp : %.5i j : %.5i border %.5e\n" testn grid.ny-testn+1 vecb_L(phL.trans_scalD[:,iscal], grid)[testn])
                            # end

                        end
                    end

                    #TODO convection_Cdivu BC divergence
                    #TODO check num.nb_transported_scalars>1

                    if ((num.current_i-1)%show_every == 0) 
                        # printstyled(color=:cyan, @sprintf "\n before scalar transport \n")
                        # print_electrolysis_statistics(num,grid,phL)
                        
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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
                    end
                
                    # printstyled(color=:red, @sprintf "\n levelset: before scalar_transport\n")
                    # println(grid.LS[1].geoL.dcap[1,1,:])

                            
                    #TODO better alloc
                    # geo = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_u = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_v = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    # laps = set_matrices!(
                    #     num, grid, geo, grid_u, geo_u, grid_v, geo_v,
                    #     op.opC_pL, op.opC_uL, op.opC_vL,
                    #     periodic_x, periodic_y)

                    if electrolysis_advection

                    # if imposed_velocity == "zero" && num.current_i ==16
                        #    num.ϵwall = 0.0
                        #   print("\n changed eps wall")
                    # end
                        #printstyled(color=:cyan, @sprintf "\n epsilon %.2e %.2e \n" num.ϵ num.ϵwall )

                        update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
                    end

                    
                    # printstyled(color=:red, @sprintf "\n return before debug mem\n")
                    # return

                    # if num.nLS ==1 #TODO dev multiple levelsets


                      #PDI (IO)      
        
                    if num.io_pdi>0

                        #or permutedims(grid.LS[iLSpdi].geoL.dcap, (3, 2, 1)) (3, 1, 2)
                        try                            
                            iLSpdi = 1 # TODO all grid.LS                
                            PDI_status = @ccall "libpdi".PDI_multi_expose("write_capacities"::Cstring,                    
                            # "dcap"::Cstring, permutedims(grid.LS[iLSpdi].geoL.dcap, (3, 2, 1))::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_1"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_2"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_3"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                            "dcap_4"::Cstring, grid.LS[iLSpdi].geoL.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint,

                            C_NULL::Ptr{Cvoid})::Cint                           
                        catch error
                            printstyled(color=:red, @sprintf "\n PDI error \n")
                            print(error)
                            printstyled(color=:red, @sprintf "\n PDI error \n")
                        end
                    end #if io_pdi


                   
                    # # #TODO better alloc
                    # geo = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_u = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                    # geo_v = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]

                    # laps = set_matrices!(
                    #     num, grid, geo, grid_u, geo_u, grid_v, geo_v,
                    #     op.opC_pL, op.opC_uL, op.opC_vL,
                    #     periodic_x, periodic_y)

                    #interface term in rhs comes from op.opC_TL
                    #TODO variable CL

                    scalar_transport!(num, grid, grid_u, grid_v,
                    op.opC_TL,
                    op.opL,
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
                    periodic_x, 
                    periodic_y, 
                    electrolysis_convection, 
                    ls_advection)

                    # PDI_status = @ccall "libpdi".PDI_multi_expose("check_concentrations"::Cstring,
                    # "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    # C_NULL::Ptr{Cvoid})::Cint

                    # scalar_transport!(BC_trans_scal, num, grid, , grid.LS[1].geoL, phL, num.concentration0,
                    # grid.LS[1].MIXED, grid.LS[1].geoL.projection, op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    # periodic_x, periodic_y, electrolysis_convection, ls_advection, BC_int, num.diffusion_coeff,Ascal,Bscal,all_CUTCT,rhs_scal)

                    # else
                    #     printstyled(color=:red, @sprintf "\n TODO multiple LS \n" )

                    # end

                    # scalar_transport_2!(BC_trans_scal, num, grid, op.opC_TL, grid.LS[1].geoL, phL, num.concentration0,
                    # grid.LS[1].MIXED, grid.LS[1].geoL.projection, op.opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    # periodic_x, periodic_y, electrolysis_convection, true, BC_int, num.diffusion_coeff)

                    if ((num.current_i-1)%show_every == 0) 
                        # printstyled(color=:cyan, @sprintf "\n after scalar transport \n")
                        # print_electrolysis_statistics(num,grid,phL)
                        
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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


                        # scal_error=0.0
                        # for iscal in 1:num.nb_transported_scalars
                        #     # printstyled(color=:cyan, @sprintf "\n after scalar transport %.2i %.2e \n" iscal maximum(abs.(phL.trans_scal[:,:,iscal])))
                        #     # printstyled(color=:cyan, @sprintf "\n after scalar transport %.2i %.2e \n" iscal maximum(abs.(phL.trans_scalD[:,iscal])))
                        #     # scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-num.concentration0[iscal])./num.concentration0[iscal])
                        #     # scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid).-num.concentration0[iscal])./num.concentration0[iscal])
                        #     # scal_error = max(scal_error_bulk,scal_error_border,scal_error)
                        # end

                        # printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e num.CFL %.2e\n" scal_error num.v_inlet*num.dt0/grid.dx[1,1])


                    end

                    if imposed_velocity != "none" && num.debug== "scalar_testing"
                        scal_error=0.0
                        for iscal in 1:num.nb_transported_scalars

                            # print("\n maximum ",maximum(phL.trans_scalD[:,iscal]), )
                            # printstyled(color=:cyan, @sprintf "\n error after scalar transport max %.2e min %.2e\n" maximum(phL.trans_scalD[:,iscal]) minimum(phL.trans_scalD[:,iscal]))

                            scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-num.concentration0[iscal])./num.concentration0[iscal])
                            scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid).-num.concentration0[iscal])./num.concentration0[iscal])
                            scal_error = max(scal_error_bulk,scal_error_border,scal_error)

                        end

                        printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e num.CFL %.2e\n" scal_error num.v_inlet*num.dt0/grid.dx[1,1])

                        # printstyled(color=:red, @sprintf "\n Poiseuille \n")

                        # # Check the velocity field before the scalar transport
                        # test_Poiseuille(num,phL.vD,grid_v)

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

                        # compute_grad_p!(num,grid, grid_u, grid_v, phL.pD, op.opC_pL, op.opC_uL, op.opC_vL)

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
        
                time = current_t #Cdouble
                nstep = num.current_i
            
                # phi_array=phL.phi_ele #do not transpose since python row major
                
                # Compute electrical current, interpolate velocity on scalar grid
                #                     if num.electrical_potential>0 compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
        
                # if electrolysis && num.nb_transported_scalars>1
                #     if heat 
                #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru.*phL.T) #phL.T
                #     else
                #         elec_cond = 2*num.Faraday^2 .*phL.trans_scal[:,:,2].*num.diffusion_coeff[2]./(num.Ru*num.temperature0) 
                #     end
                # else 
                #     elec_cond = ones(grid)
                #     printstyled(color=:green, @sprintf "\n conductivity one")

                # end 

                # phL.i_current_mag .*= elec_cond # i=-κ∇ϕ here magnitude


                #store in us, vs instead of Eus, Evs
                # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

                # printstyled(color=:red, @sprintf "\n test i current\n" )
                # print("\n phi ",phL.phi_ele[64,127]," ",phL.phi_ele[64,128]," ",(phL.phi_ele[64,128]-phL.phi_ele[64,127])/grid.dx[64,128]," ",(phL.phi_ele[64,128]-phL.phi_ele[64,127])/grid.dx[64,127]*elec_cond[64,127], " cond ", elec_cond[64,127]," \n")


                # tmp_vec_p .*= elec_cond
                # tmp_vec_p0 .*= elec_cond


                # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # C_NULL::Ptr{Cvoid})::Cint
            
                interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
                    
                iLSpdi = 1 # TODO all grid.LS

                # Exposing data to PDI for IO    
                # if writing "D" array (bulk, interface, border), add "_1D" to the name

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                C_NULL::Ptr{Cvoid})::Cint
        
                
                #TODO debug with volume fraction

                # A = zeros(gv.ny, gv.nx+1)
                # for jplot in 1:gv.ny
                #     for iplot in 1:gv.nx+1
                #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                #     pII = lexicographic(II, grid.ny + 1)
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

        
        #region Phase change

        print("\n electrolysis_phase_change_case ",electrolysis_phase_change_case)

        if electrolysis && electrolysis_phase_change_case != "None"
            printstyled(color=:magenta, @sprintf "\n integrate_mass_flux_over_interface\n")

            # @views integrate_mass_flux_over_interface(num,grid,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,mass_flux)

            # @views integrate_mass_flux_over_interface_2(num,grid,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,mass_flux)

            @views integrate_mass_flux_over_interface(num,grid,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,
            mass_flux_vecb,mass_flux_veci, tmp_vec_p, tmp_vec_p0, tmp_vec_p1, mass_flux,num.index_phase_change) #1

            print("\n sum mass flux all levelsets (walls and interfaces alike) ", sum(mass_flux),"\n ")
        end
        
        #    grid.LS[i].α  which is the angle of the outward point normal with respect to the horizontal axis

        for iLS in 1:num.nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, grid.LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
            elseif is_fs(BC_int[iLS]) || (occursin("levelset",electrolysis_phase_change_case) && iLS == iLSbubble)
                printstyled(color=:green, @sprintf "\n grid.V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(grid.V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))

                if electrolysis_phase_change_case!="none"    
                    if occursin("levelset",electrolysis_phase_change_case)

                        printstyled(color=:magenta, @sprintf "\n phase-change for LS %.2i " iLS)

                        # plot_electrolysis_velocity!(num, grid, grid.LS, grid.V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)

                        # TODO send to PDI points and velocity for phase change like in plot_electrolysis_velocity!
                        
                        # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 

                       
                            
                        flower_status = update_free_surface_velocity_electrolysis!(num, grid, grid_u, grid_v, iLS, phL.uD, phL.vD, 
                        periodic_x, periodic_y, num.average_velocity, phL.trans_scalD[:,num.index_phase_change],phL.trans_scal[:,:,num.index_phase_change],
                        num.diffusion_coeff[num.index_phase_change],num.concentration0[num.index_phase_change],electrolysis_phase_change_case,mass_flux)

                        if flower_status !=0
                            printstyled(color=:red, @sprintf "\n Stopping simulation %.3i " flower_status)
                            return
                        end
                            
                            # # iLS = 1
                            # # intfc_length = 0.0
                            # # @inbounds @threads for II in grid.LS[iLS].MIXED
                            # #     intfc_length += 
                            # # end


                            # printstyled(color=:green, @sprintf "\n pi*R %.2e len : %.2e \n" π*num.R intfc_length)

                            # #TODO u-vphase change

                            # #TODO check velocity
                            # @inbounds @threads for II in grid.LS[iLS].MIXED
                            #     grid.V[II] = sum(mass_flux) * num.diffusion_coeff[num.index_phase_change] *(1.0/num.rho2-1.0/num.rho1).*num.diffusion_coeff[num.index_phase_change].*num.MWH2
                            # end


                        if num.mass_flux == 0
                            varnH2 = num.sum_mass_flux * num.diffusion_coeff[num.index_phase_change] 

                            new_nH2 = nH2 + varnH2 * num.τ

                            print("\n varn ",varnH2 ," dt ", num.τ," dn ",varnH2 * num.τ, " sum ", num.sum_mass_flux)
                            printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new nH2 %.2e \n" num.current_i nH2 varnH2*num.τ new_nH2)

                            if varnH2 < 0.0 
                                # print(@sprintf "error nH2 %.2e dnH2 %.2e new nH2 %.2e\n" nH2-varnH2*num.τ varnH2*num.τ nH2 )
                                @error ("error nH2")
                                crashed = true
                                new_nH2 = nH2
                                print("wrong nH2 ")
                                # println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                                return
                            else
                                nH2 = new_nH2
                            end

                         end #num.mass_flux == 0


                    end
                    # update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
                    printstyled(color=:green, @sprintf "\n grid.V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(grid.V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))
                    
                    printstyled(color=:green, @sprintf "\n grid.V %.2e dx : %.2e CFL %.2e\n" maximum(abs.(grid.V)) grid.dx[1,1] maximum(abs.(grid.V))*num.τ/grid.dx[1,1])


                else
                    update_free_surface_velocity(num, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y)
                end

            
            elseif (electrolysis && occursin("Khalighi",electrolysis_phase_change_case))

                if ((num.current_i-1)%show_every == 0) 
                    # print_electrolysis_statistics(num,grid,phL)
                    PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                    "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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
                end

                previous_radius = num.current_radius

                # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 
                varnH2 =  sum(mass_flux) * num.diffusion_coeff[num.index_phase_change] 

                #TODO mode_2d==0 flux corresponds to cylinder of length 1
                #2D cylinder reference length
                if mode_2d == 1
                    varnH2 .*= num.ref_thickness_2d
                end

                #Pliquid is the average value of p over the bubble interface plus the ambient operating pressure (P).
                p_liq= num.pres0 + mean(veci(phL.pD,grid,2)) #TODO here one bubble
                # p_g=p_liq + 2 * num.σ / num.current_radius #3D
                p_g=p_liq + num.σ / num.current_radius #2D

                new_nH2 = nH2 + varnH2 * num.τ

               
                printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new nH2 %.2e \n" num.current_i nH2 varnH2*num.τ new_nH2)

                if varnH2 < 0.0 
                    # print(@sprintf "error nH2 %.2e dnH2 %.2e new nH2 %.2e\n" nH2-varnH2*num.τ varnH2*num.τ nH2 )
                    @error ("error nH2")
                    crashed = true
                    new_nH2 = nH2
                    print("wrong nH2 ")
                    # println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                    return
                end

                if occursin("Khalighi_no_update",electrolysis_phase_change_case)

                else
                    nH2 = new_nH2
                end
                
                #TODO using num.temperature0
                if mode_2d == 0
                    num.current_radius = cbrt(3.0 * nH2 * num.Ru * num.temperature0/( 4.0 * pi * p_g) )
                elseif mode_2d == 1
                    num.current_radius = sqrt(nH2 * num.Ru * num.temperature0/( pi * p_g * num.ref_thickness_2d) )
                elseif mode_2d == 2
                    num.current_radius = sqrt(nH2/(num.concentration0[num.index_phase_change] * pi))
                elseif mode_2d == 3
                    num.current_radius = sqrt(2*nH2/(num.concentration0[num.index_phase_change] * pi))
                elseif mode_2d == 4 #TODO
                    num.current_radius = sqrt(nH2 * num.Ru * num.temperature0/( pi * p_g) )
                end

                printstyled(color=:green, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid.nx))

                if (num.current_radius-previous_radius)/(num.L0/grid.nx) > 0.5
                    printstyled(color=:red, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid.nx))
                    @error ("num.CFL radius")
                    crashed = true
                    return
                end

               
                printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e %.2e\n" grid.nx num.τ num.L0/grid.nx (num.current_radius-previous_radius)/(num.L0/grid.nx) sum(mass_flux))
                
                printstyled(color=:green, @sprintf "\n num.n(H2): %.2e added %.2e old R %.2e new R %.2e \n" nH2 varnH2*num.τ previous_radius num.current_radius)
                printstyled(color=:green, @sprintf "\n p0: %.2e p_liq %.2e p_lapl %.2e \n" num.pres0 p_liq p_g)

                if mode_2d == 3
                    grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)                  
                else
                    grid.LS[1].u .= sqrt.((grid.x .- num.xcoord .- num.current_radius .+ num.R ).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)
                end
                # init_franck!(grid, TL, R, num.T_inf, 0)
                # u

            elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius")

                #num.CFL 0.5
                num.current_radius = num.current_radius + grid.dx[1,1]/2

                grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)                  

            elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius4")

                #num.CFL 0.5
                num.current_radius = num.current_radius + grid.dx[1,1]/4

                grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)                  


            end #phase change

        end #iLS

        #endregion

        if verbose && adaptative_t
            println("num.τ = $num.τ")
        end


        #region Advection 
        if advection

            if num.io_pdi>0

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_before_LS_adv"::Cstring,
                # "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                # "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                # "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                # "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                # "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                "normal_velocity_intfc"::Cstring, grid.V::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "velocity_x_intfc"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "velocity_y_intfc"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                # "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                C_NULL::Ptr{Cvoid})::Cint

            end # if num.io_pdi>0


            #In the current implementations, the first cell is not solved: levelset not advected
            for (iLS, bc) in enumerate(BC_int)
                if is_stefan(bc)
                    IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                    grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u)), grid)
                    # u .= sqrt.((grid.x .- num.current_i*num.Δ/1).^ 2 + grid.y .^ 2) - (0.5) * ones(grid.nx, grid.ny);
                elseif is_fs(bc) || (occursin("levelset",electrolysis_phase_change_case) && iLS == iLSbubble)

                    if num.advection_LS_mode == 0
                        rhs_LS .= 0.0
                        grid.LS[iLS].A.nzval .= 0.0
                        grid.LS[iLS].B.nzval .= 0.0
                        IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        utmp .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                        rhs_LS .= 0.0
                        S2IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, utmp, grid.LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                    elseif num.advection_LS_mode == 1

                        # Project velocities to the normal and use advection scheme for advection just
                        # in the normal direction
                        tmpVx = zeros(grid)
                        tmpVy = zeros(grid)
                        grid.V .= 0.0
                        @inbounds @threads for II in grid.LS[iLS].MIXED
                            cap1 = grid_u.LS[iLS].geoL.cap[II,5]
                            cap3 = grid_u.LS[iLS].geoL.cap[δx⁺(II),5]
                            tmpVx[II] = (grid_u.V[II] * cap1 + grid_u.V[δx⁺(II)] * cap3) / (cap1 + cap3 + eps(0.01))

                            cap2 = grid_v.LS[iLS].geoL.cap[II,5]
                            cap4 = grid_v.LS[iLS].geoL.cap[δy⁺(II),5]
                            tmpVy[II] = (grid_v.V[II] * cap2 + grid_v.V[δy⁺(II)] * cap4) / (cap2 + cap4 + eps(0.01))

                            tmpV = sqrt(tmpVx[II]^2 + tmpVy[II]^2)
                            β = atan(tmpVy[II], tmpVx[II])
                            if grid.LS[iLS].α[II] > 0.0 && β < 0.0
                                β += 2π
                            end
                            if grid.LS[iLS].α[II] < 0.0 && β > 0.0
                                β -= 2π
                            end

                            grid.V[II] = tmpV * cos(β - grid.LS[iLS].α[II])
                        end

                        i_ext, l_ext, b_ext, r_ext, t_ext = indices_extension(grid, grid.LS[iLS], grid.ind.inside, periodic_x, periodic_y)
                        field_extension!(grid, grid.LS[iLS].u, grid.V, i_ext, l_ext, b_ext, r_ext, t_ext, num.NB, periodic_x, periodic_y)

                        rhs_LS .= 0.0
                        IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                        BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        BC_LS_interior!(num, grid, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                        # Impose contact angle if a wall is present
                        rhs_LS .= 0.0
                        grid.LS[iLS].A.nzval .= 0.0
                        grid.LS[iLS].B.nzval .= 0.0
                        for II in grid.ind.all_indices
                            pII = lexicographic(II, grid.ny)
                            grid.LS[iLS].A[pII,pII] = 1.0
                            grid.LS[iLS].B[pII,pII] = 1.0
                        end
                        BC_LS_interior!(num, grid, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                    elseif num.advection_LS_mode == 2
                        print("\n num.advection_LS_mode == 2 iLS", iLS)


                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)


                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                        IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                      
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u)), grid)



                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)


                    elseif num.advection_LS_mode == 3

                        # from scalar grid, normal to u and v
                        #TODO 
                        print("\n setting 0 u and v velocities \n")
                        grid_u.V .= 0
                        grid_v.V .= 0

                        interpolate_scalar!(grid, grid_u, grid_v, grid.V, grid_u.V, grid_v.V)

                        normalx = cos.(grid_u.LS[iLS].α)
                        normaly = sin.(grid_v.LS[iLS].α)
                    
                        grid_u.V .*= normalx
                        grid_v.V .*= normaly

                        print("\n num.advection_LS_mode == 3 iLS", iLS)


                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)


                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))


                        rhs_LS .= 0.0
                        grid.LS[iLS].A.nzval .= 0.0
                        grid.LS[iLS].B.nzval .= 0.0
                        IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        utmp .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                        rhs_LS .= 0.0
                        S2IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, utmp, grid.LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                        BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)


                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)

                    elseif num.advection_LS_mode == 4

                        previous_radius = num.current_radius

                        # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 
                        varnH2 = sign_mass_flux * sum(mass_flux) * num.diffusion_coeff[num.index_phase_change] 

                        #TODO mode_2d==0 flux corresponds to cylinder of length 1
                        #2D cylinder reference length
                        if mode_2d == 1
                            varnH2 .*= num.ref_thickness_2d
                        end

                      
                        nH2 = nH2 + varnH2 * num.τ

                        # printstyled(color=:green, @sprintf "\n it %.5i Mole: %.2e dn %.2e new nH2 %.2e \n" num.current_i nH2 varnH2*num.τ new_nH2)

                        # if varnH2 < 0.0 
                        #     # print(@sprintf "error nH2 %.2e dnH2 %.2e new nH2 %.2e\n" nH2-varnH2*num.τ varnH2*num.τ nH2 )
                        #     @error ("error nH2")
                        #     crashed = true
                        #     new_nH2 = nH2
                        #     print("wrong nH2 ")
                        #     # println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                        #     return
                        # end
        
                        # if occursin("Khalighi_no_update",electrolysis_phase_change_case)
        
                        # else
                        #     nH2 = new_nH2
                        # end
                        
                        #TODO using num.temperature0
                        if mode_2d == 0
                            num.current_radius = cbrt(3.0 * nH2 * num.Ru * num.temperature0/( 4.0 * pi * p_g) )
                        elseif mode_2d == 1
                            num.current_radius = sqrt(nH2 * num.Ru * num.temperature0/( pi * p_g * num.ref_thickness_2d) )
                        elseif mode_2d == 2
                            num.current_radius = sqrt(nH2/(num.concentration0[num.index_phase_change] * pi))
                        elseif mode_2d == 3
                            num.current_radius = sqrt(2*nH2/(num.concentration0[num.index_phase_change] * pi))
                        end
        
                        printstyled(color=:green, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid.nx))
        
                        if (num.current_radius-previous_radius)/(num.L0/grid.nx) > 0.5
                            printstyled(color=:red, @sprintf "\n radius num.CFL: %.2e \n" (num.current_radius-previous_radius)/(num.L0/grid.nx))
                            @error ("num.CFL radius")
                            crashed = true
                            return
                        end
        
                       
                        printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e \n" grid.nx num.τ num.L0/grid.nx (num.current_radius-previous_radius)/(num.L0/grid.nx)) 
                        # sum(mass_flux))
                        
                        printstyled(color=:green, @sprintf "\n num.n(H2): %.2e added %.2e old R %.2e new R %.2e \n" nH2 varnH2*num.τ previous_radius num.current_radius)
                        printstyled(color=:green, @sprintf "\n p0: %.2e p_liq %.2e p_lapl %.2e \n" num.pres0 p_liq p_g)
        
                        if mode_2d == 3
                            grid.LS[1].u .= sqrt.((grid.x.- num.xcoord).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)                  
                        else
                            grid.LS[1].u .= sqrt.((grid.x .- num.xcoord .- num.current_radius .+ num.R ).^ 2 + (grid.y .- num.ycoord) .^ 2) - (num.current_radius) * ones(grid.ny, grid.nx)
                        end


                    elseif num.advection_LS_mode == 5
                        print("\n num.advection_LS_mode == 5 iLS", iLS)


                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)

                        grid.V .=0.25*grid.dx[1,1]/num.τ  

                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                        IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u)), grid)



                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)

                    elseif (num.advection_LS_mode == 6) || (num.advection_LS_mode == 7)

                        rhs_LS .= 0.0
                        # grid.LS[iLS].A.nzval .= 0.0
                        # grid.LS[iLS].B.nzval .= 0.0

                        print("\n num.advection_LS_mode == 2 iLS", iLS)


                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)


                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                        IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)

                        # IIOE_normal_indices!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)
                        # grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u)), grid)


                        # rhs_LS .= 0.0
                        # grid.LS[iLS].A.nzval .= 0.0
                        # grid.LS[iLS].B.nzval .= 0.0
                        # IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                        # BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        # BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        # utmp .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                        # rhs_LS .= 0.0
                        # S2IIOE!(grid, grid_u, grid_v, grid.LS[iLS].A, grid.LS[iLS].B, utmp, grid.LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                        # BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        # BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                        # grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                       

                        if num.advection_LS_mode == 6

                            # IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                            BC_LS!(grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)
                            BC_LS_interior!(num, grid, grid_u, grid_v, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        
                        end
                        
                        grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)

                        # # Impose contact angle if a wall is present
                        # rhs_LS .= 0.0
                        # grid.LS[iLS].A.nzval .= 0.0
                        # grid.LS[iLS].B.nzval .= 0.0
                        # for II in grid.ind.all_indices
                        #     pII = lexicographic(II, grid.ny)
                        #     grid.LS[iLS].A[pII,pII] = 1.0
                        #     grid.LS[iLS].B[pII,pII] = 1.0
                        # end
                        # BC_LS_interior!(num, grid, iLS, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                        # grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)





                        print_CL_length(num,grid, grid.LS[iLS].u, grid.LS[iLS].A, grid.LS[iLS].B, rhs_LS, BC_u)


                    elseif num.advection_LS_mode == 8
                        print("\n num.advection_LS_mode == 8 iLS", iLS)

                        nghost = 1

                        Aghost, Bghost = allocate_ghost_matrices(grid.nx,grid.ny,nghost)

                        print("\n periodic ",periodic_x," y ",periodic_y)


                        print_CL_length(num,grid, grid.LS[iLS].u, Aghost, Bghost, rhs_LS, BC_u)


                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                        # IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                        # IIOE_normal_indices!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # print("\n sizes LS A ",size(grid.LS[iLS].A), " B ", size(grid.LS[iLS].B)," LS ",size(grid.LS[iLS].u)," V ",size(grid.V),"\n")

                        #recopy value in ghost cell
                        LSghost = init_ghost_neumann(grid.LS[iLS].u,grid.nx,grid.ny,nghost)
                        
                        Vghost = init_ghost_neumann(grid.V,grid.nx,grid.ny,nghost)

                        # print("\n LSghost \n")
                        # print("\n ",LSghost[0,:])
                        # print("\n ",LSghost[1,:])
                        # print("\n ",LSghost[:,0])
                        # print("\n ",LSghost[:,1])

                        IIOE_normal_indices!(grid, Aghost, Bghost, grid.LS[iLS].u, LSghost, 
                        grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # IIOE_normal_indices!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, LSghost, 
                        # grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # IIOE_normal_indices!(grid, Aghost, Bghost, LSghost, Vghost, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # print("\n sizes LS A ",size(grid.LS[iLS].A), " B ", size(grid.LS[iLS].B)," LS ",size(grid.LS[iLS].u)," V ",size(grid.V),"\n")
                        # print("\n sizes LS A ",size(OffsetArrays.no_offset_view(Aghost)), " B ", size(OffsetArrays.no_offset_view(Bghost))," LS ",size(LSghost)," V ",size(grid.V),"\n")

                        OffsetArrays.no_offset_view(LSghost) .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(OffsetArrays.no_offset_view(LSghost))), (grid.ny+2,grid.nx+2))

                        grid.LS[iLS].u .= LSghost[1:grid.ny,1:grid.nx]

                        # grid.LS[iLS].u .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(grid.LS[iLS].u)), grid)



                        print_CL_length(num,grid, grid.LS[iLS].u, Aghost, Bghost, rhs_LS, BC_u)

                    elseif ((num.advection_LS_mode == 9) || (num.advection_LS_mode == 10))
                        print("\n num.advection_LS_mode == 8 iLS", iLS)

                        if num.advection_LS_mode == 10
                            grid.V .=0.25*grid.dx[1,1]/num.τ  
                        end

                        nghost = 1

                        Aghost, Bghost = allocate_ghost_matrices_2(grid.nx,grid.ny,nghost)


                        print_CL_length(num,grid, grid.LS[iLS].u, Aghost, Bghost, rhs_LS, BC_u)


                        printstyled(color=:green, @sprintf "\n grid p u v max : %.2e %.2e %.2e\n" maximum(abs.(grid.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_u.V[grid.LS[iLS].MIXED])) maximum(abs.(grid_v.V[grid_v.LS[iLS].MIXED])))

                        # IIOE_normal!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y)
                        # IIOE_normal_indices!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # print("\n sizes LS A ",size(grid.LS[iLS].A), " B ", size(grid.LS[iLS].B)," LS ",size(grid.LS[iLS].u)," V ",size(grid.V),"\n")


                        LSghost = init_ghost_neumann_2(grid.LS[iLS].u,grid.nx,grid.ny,nghost)

                        # print("\n LSghost \n")
                        # print("\n ",LSghost[0,:])
                        # print("\n ",LSghost[1,:])
                        # print("\n ",LSghost[:,0])
                        # print("\n ",LSghost[:,1])
                                               
                        Vghost = init_ghost_neumann_2(grid.V,grid.nx,grid.ny,nghost)

                        IIOE_normal_indices_2!(grid, Aghost, Bghost, grid.LS[iLS].u, LSghost, 
                        Vghost, CFL_sc, periodic_x, periodic_y,nghost)

                        # IIOE_normal_indices!(grid, grid.LS[iLS].A, grid.LS[iLS].B, grid.LS[iLS].u, LSghost, 
                        # grid.V, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # IIOE_normal_indices!(grid, Aghost, Bghost, LSghost, Vghost, CFL_sc, periodic_x, periodic_y,grid.ind.all_indices)

                        # print("\n sizes LS A ",size(grid.LS[iLS].A), " B ", size(grid.LS[iLS].B)," LS ",size(grid.LS[iLS].u)," V ",size(grid.V),"\n")
                        # print("\n sizes LS A ",size(OffsetArrays.no_offset_view(Aghost)), " B ", size(OffsetArrays.no_offset_view(Bghost))," LS ",size(LSghost)," V ",size(grid.V),"\n")

                        LSghost .= reshape(gmres(Aghost, Bghost * vec(LSghost)), (grid.ny+2*nghost,grid.nx+2*nghost))


                        # grid.LS[iLS].u .= LSghost[1:grid.ny,1:grid.nx]

                        #Store result of LS advection without ghost cells
                        for j=1:grid.ny
                            for i=1:grid.nx
                                grid.LS[iLS].u[j,i] = LSghost[j+1,i+1]
                            end
                        end

                        # grid.LS[iLS].u .= reshape(gmres(OffsetArrays.no_offset_view(Aghost), OffsetArrays.no_offset_view(Bghost) * vec(grid.LS[iLS].u)), grid)



                        print_CL_length(num,grid, grid.LS[iLS].u, Aghost, Bghost, rhs_LS, BC_u)



                    end #num.advection_LS_mode == 
                end
            end
            if analytical
                u[grid.ind.b_top[1]] .= sqrt.(grid.x[grid.ind.b_top[1]] .^ 2 + grid.y[grid.ind.b_top[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[grid.ind.b_bottom[1]] .= sqrt.(grid.x[grid.ind.b_bottom[1]] .^ 2 + grid.y[grid.ind.b_bottom[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[grid.ind.b_left[1]] .= sqrt.(grid.x[grid.ind.b_left[1]] .^ 2 + grid.y[grid.ind.b_left[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
                u[grid.ind.b_right[1]] .= sqrt.(grid.x[grid.ind.b_right[1]] .^ 2 + grid.y[grid.ind.b_right[1]] .^ 2) .- (num.R + speed*num.current_i*num.τ);
            elseif num.nb_reinit > 0
                if auto_reinit == 1 && (num.current_i-1)%num.reinit_every == 0
                    for iLS in 1:num.nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg, rl_rg_v = rg(num, grid, grid.LS[iLS].u, periodic_x, periodic_y, BC_int)
                            println("$(ls_rg)")
                            printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
                            if ls_rg >= num.δreinit || num.current_i == 1
                                print("(ls_rg >= num.δreinit || num.current_i == 1): yes")
                                # println("yes")
                                RK2_reinit!(ls_scheme, grid, grid.ind, iLS, grid.LS[iLS].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                                
                                ls_rg, rl_rg_v = rg(num, grid, grid.LS[iLS].u, periodic_x, periodic_y, BC_int)
                                println("$(ls_rg) ")
                                printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
                            end
                        end
                    end
                elseif (num.current_i-1)%num.reinit_every == 0
                    for iLS in 1:num.nLS
                        if !is_wall(BC_int[iLS])
                            RK2_reinit!(ls_scheme, grid, grid.ind, iLS, grid.LS[iLS].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                        end
                    end
                # elseif num.nLS > 1
                #     for iLS in 1:num.nLS
                #         if !is_wall(BC_int[iLS])
                #             RK2_reinit!(ls_scheme, grid, grid.ind, iLS, grid.LS[iLS].u, 2num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int, true)
                #         end
                #     end
                        end
                    end

            # Numerical breakup
            if free_surface && breakup ==1
                count, id_break = breakup_n(grid.LS[1].u, grid.nx, grid.ny, grid.dx, grid.dy, periodic_x, periodic_y, NB_indices, 5e-2)
                println(count)
                if count > count_limit_breakup
                    println("BREAK UP!!") 
                    breakup_f(grid, grid.LS[1].u, id_break)
                    RK2_reinit!(ls_scheme, grid, grid.ind, 1, grid.LS[1].u, num.nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                end
            end
        end

        if verbose
            if (num.current_i-1)%show_every == 0
                printstyled(color=:green, @sprintf "\n Current iteration : %d (%d%%) | t = %.2e \n" (num.current_i-1) 100*(num.current_i-1)/num.max_iterations current_t)
                printstyled(color=:green, @sprintf "\n num.CFL : %.2e num.CFL : %.2e num.τ : %.2e\n" num.CFL max(abs.(grid.V)..., abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)*num.τ/num.Δ num.τ)
                if heat && length(grid.LS[end].MIXED) != 0
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" mean(grid.V[grid.LS[1].MIXED]) findmax(grid.V[grid.LS[1].MIXED])[1] findmin(grid.V[grid.LS[1].MIXED])[1])
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(grid.LS[1].κ[grid.LS[1].MIXED]) findmax(grid.LS[1].κ[grid.LS[1].MIXED])[1] findmin(grid.LS[1].κ[grid.LS[1].MIXED])[1])
                elseif advection && length(grid.LS[end].MIXED) != 0
                    V_mean = mean([mean(grid_u.V[grid.LS[1].MIXED]), mean(grid_v.V[grid.LS[1].MIXED])])
                    V_max = max(findmax(grid_u.V[grid.LS[1].MIXED])[1], findmax(grid_v.V[grid.LS[1].MIXED])[1])
                    V_min = min(findmin(grid_u.V[grid.LS[1].MIXED])[1], findmin(grid_v.V[grid.LS[1].MIXED])[1])
                    # print(@sprintf "Vol_ratio = %.3f%%\n" (volume(grid.LS[end].geoL) / V0L * 100))
                    print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e\n" V_mean V_max V_min)
                    print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e\n" mean(grid.LS[1].κ[grid.LS[1].MIXED]) findmax(grid.LS[1].κ[grid.LS[1].MIXED])[1] findmin(grid.LS[1].κ[grid.LS[1].MIXED])[1])
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
                            # print_electrolysis_statistics(num,grid,phL) 
                            PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                            "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                            "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                            "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                            "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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
                printstyled(color=:red, @sprintf "\n grid.LS not updated \n")
                print(errorLS)
                return
            end
            # printstyled(color=:red, @sprintf "\n levelset 4:\n")
            # println(grid.LS[1].geoL.dcap[1,1,:])

            grid.LS[end].geoL.fresh .= false
            grid.LS[end].geoS.fresh .= false
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

            if iszero(num.current_i%num.save_every) || num.current_i==num.max_iterations
                snap = num.current_i÷num.save_every+1
                if save_radius
                    radius[snap] = find_radius(grid, grid.LS[1])
                end
                if hill
                    a = zeros(length(grid.LS[1].MIXED))
                    for i in eachindex(grid.LS[1].MIXED)
                        a[i] = grid.LS[1].geoL.projection[grid.LS[1].MIXED[i]].pos.y
                    end
                    radius[snap] = mean(a)
                end
                # if save_length
                #     fwd.length[snap] = arc_length2(grid.LS[1].geoS.projection, grid.LS[1].MIXED)
                # end
            end
        end



        intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num, intfc_seg_num = convert_interfacial_D_to_segments(num,grid,phL.TD,1,2)
        # print("\n number of interface points intfc_vtx_num ", intfc_vtx_num)
        # print("\n intfc_vtx_connectivities ",intfc_vtx_connectivities)
        # print("\n len ", size(intfc_vtx_connectivities),intfc_seg_num)

        # print("\n intfc_vtx_x ",intfc_vtx_x)
        # print("\n intfc_vtx_x ",intfc_vtx_y)

        barycenter_x_coord = mean(intfc_vtx_x)

        PDI_status = @ccall "libpdi".PDI_multi_expose("update_levelset"::Cstring,
        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        # "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        # "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
        # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
        # "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
        "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
        "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
        "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
        "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
        "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
        "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
        "barycenter_x_coord"::Cstring, barycenter_x_coord::Ref{Cdouble}, PDI_OUT::Cint,
        C_NULL::Ptr{Cvoid})::Cint

        #endregion Advection 


        #region Navier-Stokes

        if navier_stokes
            # if !advection
            #     @time no_slip_condition!(num, grid, grid_u, grid_u.LS[1], grid_v, grid_v.LS[1], periodic_x, periodic_y)
            #     # grid_u.V .= num.Δ / (1 * num.τ)
            #     # grid_v.V .= 0.0
            # end

            # Pressure-velocity coupling

            if ns_solid_phase
                geoS = [grid.LS[iLS].geoS for iLS in 1:num._nLS]
                geo_uS = [grid_u.LS[iLS].geoS for iLS in 1:num._nLS]
                geo_vS = [grid_v.LS[iLS].geoS for iLS in 1:num._nLS]
                Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,Mm1_S, Mum1_S, Mvm1_S, Cum1S, Cvm1S = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoS, grid_u, geo_uS, grid_v, geo_vS, phS,
                    BC_uS, BC_vS, BC_pS,
                    op.opC_pS, op.opC_uS, op.opC_vS, op.opS,
                    AuS, BuS, AvS, BvS, AϕS, AuvS, BuvS,
                    Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, Lum1_S, bc_Lum1_S, bc_Lum1_b_S, Lvm1_S, bc_Lvm1_S, bc_Lvm1_b_S,
                    Cum1S, Cvm1S, Mum1_S, Mvm1_S,
                    periodic_x, periodic_y, ns_advection, advection, num.current_i, Ra, navier,pres_free_surfaceS,jump_mass_fluxS,mass_fluxS
                )
            end
            if ns_liquid_phase
                geoL = [grid.LS[iLS].geoL for iLS in 1:num._nLS]
                geo_uL = [grid_u.LS[iLS].geoL for iLS in 1:num._nLS]
                geo_vL = [grid_v.LS[iLS].geoL for iLS in 1:num._nLS]
                Lpm1_L, bc_Lpm1_L, bc_Lpm1_b_L, Lum1_L, bc_Lum1_L, bc_Lum1_b_L, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L, Mm1_L, Mum1_L, Mvm1_L, Cum1L, Cvm1L = pressure_projection!(
                    time_scheme, BC_int,
                    num, grid, geoL, grid_u, geo_uL, grid_v, geo_vL, phL,
                    BC_uL, BC_vL, BC_pL,
                    op.opC_pL, op.opC_uL, op.opC_vL, op.opL,
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
                #     num, grid, grid.LS[1].geoL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL, phL,
                #     BC_uL, BC_vL, op.opL
                # )
            end
        end # if navier_stokes

        #endregion Navier-Stokes 


        # cD, cL, D, L = force_coefficients!(num, grid, grid_u, grid_v, op.opL, fwd, phL; step = num.current_i+1, saveCoeffs = false)

        # if iszero(num.current_i%num.save_every) || num.current_i==num.max_iterations
        #     snap = num.current_i÷num.save_every+1
        #     if num.current_i==num.max_iterations
        #         snap = size(fwd.T,1)
        #     end
        #     fwd.t[snap] = current_t
        #     @views fwd.V[snap,:,:] .= grid.V
        #     if advection
        #         fwdS.Vratio[snap] = volume(grid.LS[end].geoS) / V0S
        #         fwdL.Vratio[snap] = volume(grid.LS[end].geoL) / V0L
        #     end
        # end
        # @views fwd.Cd[num.current_i+1] = cD
        # @views fwd.Cl[num.current_i+1] = cL
        # # @views fwd.radius[num.current_i+1] = num.current_radius

        # PDI (IO)
        if electrolysis
            if num.io_pdi>0

                try
                    # printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                    time = current_t #Cdouble
                    nstep = num.current_i
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    # Compute electrical current, interpolate velocity on scalar grid
                    #                     if num.electrical_potential>0 compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current

                    if num.electrical_potential>0
                        compute_grad_phi_ele!(num, grid, grid_u, grid_v, grid_u.LS[end], grid_v.LS[end], phL, phS,
                        op.opC_pL, op.opC_pS, elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current
                    end

                    # #store in us, vs instead of Eus, Evs
                    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,tmp_vec_p,tmp_vec_p0)

                    # #TODO i_current_mag need cond

                    # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,  
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # C_NULL::Ptr{Cvoid})::Cint

                    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,Eus,Evs)
            
                    interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,tmp_vec_p,tmp_vec_p0)
                        
                    iLSpdi = 1 # TODO all grid.LS

                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                    "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    "velocity_x"::Cstring, tmp_vec_p::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, tmp_vec_p0::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, num.current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                    C_NULL::Ptr{Cvoid})::Cint
            
                    # print("\n after write \n ")
                        
                    # printstyled(color=:red, @sprintf "\n PDI test end\n" )
            
                catch error
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                    print(error)
                    printstyled(color=:red, @sprintf "\n PDI error \n")
                end

               
            end #if io_pdi

            ####################################################################################################


            if crashed #due to nH2<0...
                return
            end

        
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                any(isnan, phL.trans_scalD) || any(isnan, phL.phi_eleD) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8 || norm(phL.trans_scal) > 1e8 || norm(phL.phi_ele) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
                
                print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
    

                # print_electrolysis_statistics(num,grid,phL)

                 PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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

                crashed=true
                return

            end
        else
            if (any(isnan, phL.uD) || any(isnan, phL.vD) || any(isnan, phL.TD) || any(isnan, phS.uD) || any(isnan, phS.vD) || any(isnan, phS.TD) ||
                norm(phL.u) > 1e8 || norm(phS.u) > 1e8 || norm(phL.T) > 1e8 || norm(phS.T) > 1e8)
                println(@sprintf "\n CRASHED after %d iterations \n" num.current_i)
               
                crashed=true
                return
                
            end

        end

        #update iter number and time

        num.current_i += 1
        #update time
        current_t += num.τ

        # if adaptative_t
           
        # elseif adaptative_t
        #     num.τ = min(num.CFL*num.Δ^2*Re, num.CFL*num.Δ/max(
        #         abs.(grid.V)..., abs.(grid_u.V)..., abs.(grid_v.V)..., 
        #         abs.(phL.u)..., abs.(phL.v)..., abs.(phS.u)..., abs.(phS.v)...)
        #     )
        # end
    end

    if verbose
        try
            printstyled(color=:blue, @sprintf "\n Final iteration : %d (%d%%) | t = %.2e \n" (num.current_i-1) 100*(num.current_i-1)/num.max_iterations current_t)
            if stefan && advection
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(grid.V[grid.LS[1].MIXED]) findmax(grid.V[grid.LS[1].MIXED])[1] findmin(grid.V[grid.LS[1].MIXED])[1] std(grid.V[grid.LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(grid.LS[1].κ[grid.LS[1].MIXED]) findmax(grid.LS[1].κ[grid.LS[1].MIXED])[1] findmin(grid.LS[1].κ[grid.LS[1].MIXED])[1] std(grid.LS[1].κ[grid.LS[1].MIXED]))
            end
            if free_surface && advection
                # print(@sprintf "Vol_ratio = %.3f%%\n" (volume(grid.LS[end].geoL) / V0L * 100))
                print(@sprintf "V_mean = %.2e  V_max = %.2e  V_min = %.2e  V_stdev = %.5f\n" mean(grid.V[grid.LS[1].MIXED]) findmax(grid.V[grid.LS[1].MIXED])[1] findmin(grid.V[grid.LS[1].MIXED])[1] std(grid.V[grid.LS[1].MIXED]))
                print(@sprintf "κ_mean = %.2e  κ_max = %.2e  κ_min = %.2e  κ_stdev = %.5f\n" mean(grid.LS[1].κ[grid.LS[1].MIXED]) findmax(grid.LS[1].κ[grid.LS[1].MIXED])[1] findmin(grid.LS[1].κ[grid.LS[1].MIXED])[1] std(grid.LS[1].κ[grid.LS[1].MIXED]))
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
                        print_electrolysis_statistics(num,grid,phL)
                        PDI_status = @ccall "libpdi".PDI_multi_expose("print_variables"::Cstring,
                        "nstep"::Cstring, nstep ::Ref{Clonglong}, PDI_OUT::Cint,
                        "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_p"::Cstring, grid.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                        "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
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

                    end 
                end
            end
            print("\n\n")
        catch
            @show (length(grid.LS[end].MIXED))
        end
    end

    if levelset && (save_radius || hill)
        #TODO save radius
        return # radius
    # elseif flapping
    #     return xc, yc
    else
        return
    end
end

