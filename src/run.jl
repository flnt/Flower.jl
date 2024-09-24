
#TODO test with run_forward! or alloc grid inside or ; or no unpack


#precompile(Tuple{typeof(Core.kwcall), 
#NamedTuple{(:periodic_x, :periodic_y, :BC_uL, :BC_uS, :BC_vL, :BC_vS, :BC_pL, :BC_pS, :BC_u, :BC_int, 
# :BC_trans_scal, :BC_phi_ele, :auto_reinit, :time_scheme, :electrolysis,
# :navier_stokes, :ns_advection, :ns_liquid_phase, :verbose, :show_every,
# :electrolysis_convection, :electrolysis_liquid_phase, :electrolysis_phase_change_case, 
# :electrolysis_reaction, :imposed_velocity, :adapt_timestep_mode, :non_dimensionalize, :mode_2d, :breakup), 
# Tuple{Bool, Bool, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries,
# Flower.Boundaries, Tuple{}, Array{Flower.WallNoSlip{Float64, Float64}, 1}, Tuple{Flower.BoundariesInt, Flower.BoundariesInt, 
# Flower.BoundariesInt}, Flower.BoundariesInt, Int64, Flower.ForwardEuler, 
# Bool, Bool, Bool, Bool, Bool, Int64, Bool, Bool, String, String, String, Vararg{Int64, 4}}}, 
# typeof(Flower.run_forward), Flower.Numerical{Float64, Int64}, Flower.Mesh{Flower.GridCC, Float64, Int64}, 
# Flower.Mesh{Flower.GridFCx, Float64, Int64}, Flower.Mesh{Flower.GridFCy, Float64, Int64}, 
# Flower.DiscreteOperators{Float64, Int64}, Flower.Phase{Float64}, Flower.Phase{Float64}})

#precompile(Tuple{typeof(Core.kwcall), 
#NamedTuple{(:periodic_x, :periodic_y, :BC_uL, :BC_uS, :BC_vL, :BC_vS, :BC_pL, :BC_pS, :BC_u, :BC_int, 
# :BC_trans_scal, :BC_phi_ele, :auto_reinit, :time_scheme, :electrolysis,
# :navier_stokes, :ns_advection, :ns_liquid_phase, :verbose, :show_every,
# :electrolysis_convection, :electrolysis_liquid_phase, :electrolysis_phase_change_case, 
# :electrolysis_reaction, :imposed_velocity, :adapt_timestep_mode, :non_dimensionalize, :mode_2d, :breakup), 
# Tuple{Bool, Bool, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries, Flower.Boundaries,
# Flower.Boundaries, Tuple{}, Array{Flower.WallNoSlip{Float64, Float64}, 1}, Tuple{Flower.BoundariesInt, Flower.BoundariesInt, 
# Flower.BoundariesInt}, Flower.BoundariesInt, Int64, Flower.ForwardEuler, 
# Bool, Bool, Bool, Bool, Bool, Int64, Bool, Bool, String, String, String, Vararg{Int64, 4}}}, 
# typeof(Flower.run_forward), ,  
# , , 
# , 

function run_forward!(
    num::Numerical{Float64, Int64},
    grid::Mesh{Flower.GridCC, Float64, Int64},
    grid_u::Mesh{Flower.GridFCx, Float64, Int64},
    grid_v::Mesh{Flower.GridFCy, Float64, Int64},
    op::DiscreteOperators{Float64, Int64}, 
    phS::Phase{Float64}, 
    phL::Phase{Float64};
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
    auto_reinit = 0,
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
    save_radius = false,
    adaptative_t = false,
    breakup = 0,
    Ra = 0.0,
    λ = 1,
    electrolysis = false,
    electrolysis_convection = false,  
    electrolysis_liquid_phase = false,
    electrolysis_solid_phase = false,
    electrolysis_phase_change_case = "Khalighi",
    electrolysis_reaction = "nothing",
    imposed_velocity = "none",
    adapt_timestep_mode = 0,
    non_dimensionalize=1,
    mode_2d=0,
    test_laplacian = false,
    )

# function run_forward(
#     num, grid, grid_u, grid_v, op, phS, phL;
#     periodic_x = false,
#     periodic_y = false,
#     BC_TS = Boundaries(),
#     BC_TL = Boundaries(),
#     BC_pS = Boundaries(),
#     BC_pL = Boundaries(),
#     BC_uS = Boundaries(),
#     BC_uL = Boundaries(),
#     BC_vS = Boundaries(),
#     BC_vL = Boundaries(),
#     BC_u = Boundaries(),
#     BC_trans_scal = Vector{BoundariesInt}(),
#     BC_phi_ele = BoundariesInt(),
#     BC_int = [WallNoSlip()],
#     time_scheme = CN,
#     ls_scheme = weno5,
#     auto_reinit = 0,
#     heat = false,
#     heat_convection = false,
#     heat_liquid_phase = false,
#     heat_solid_phase = false,
#     navier_stokes = false,
#     ns_advection = false,
#     ns_liquid_phase = false,
#     ns_solid_phase = false,
#     hill = false,
#     Vmean = false,
#     levelset = true,
#     speed = 0,
#     analytical = false,
#     verbose = false,
#     show_every = 100,
#     save_radius = false,
#     adaptative_t = false,
#     breakup = 0,
#     Ra = 0.0,
#     λ = 1,
#     electrolysis = false,
#     electrolysis_convection = false,  
#     electrolysis_liquid_phase = false,
#     electrolysis_solid_phase = false,
#     electrolysis_phase_change_case = "Khalighi",
#     electrolysis_reaction = "nothing",
#     imposed_velocity = "none",
#     adapt_timestep_mode = 0,
#     non_dimensionalize=1,
#     mode_2d=0,
#     test_laplacian = false,
#     )

    # print("\n before unpack \n")

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re, max_iterations, save_every, reinit_every, nb_reinit, δreinit, ϵ, m, θ₀, aniso, nLS, _nLS, nNavier,
            concentration0, diffusion_coeff, nb_transported_scalars, temperature0, i0, phi_ele0, phi_ele1, alpha_c,
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
    num.current_i=1
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
        electrolysis_advection = true
   
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

    if occursin("levelset",electrolysis_phase_change_case)
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
        
        if electrolysis_phase_change_case != "none"
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


    
    # Initialisation
    
    #TODO which LS
    presintfc = 0.0
    # presintfc = pres0 + p_lapl ? #TODO init pressure
    #TODO perio, intfc, ... check init_fields_2!

    #No scalar in "solid" phase
    # @views init_fields_2!(phS.trans_scalD[:,iscal],phS.trans_scal[:,:,iscal],HS,BC_trans_scal[iscal],grid,concentration0[iscal])
    # @views phS.trans_scal[:,:,iscal] .= concentration0[iscal]

    Hu = zeros(grid_u) 
    get_height!(grid_u,ind,dx,dy,LS[end].geoS,Hu) #here Hu solid

    init_fields_2!(phS.uD,phS.u,Hu,BC_uS,grid_u,num.uD)
    # init_fields_2!(phS.ucorrD,phS.u,HSu,BC_uS,grid_u,num.uD)

    get_height!(grid_u,ind,dx,dy,LS[end].geoL,Hu)  #here Hu liquid

    init_fields_2!(phL.uD,phL.u,Hu,BC_uL,grid_u,num.uD)
    # init_fields_2!(phL.ucorrD,phL.u,HLu,BC_uL,grid_u,num.uD)


    Hv = zeros(grid_v) 

    get_height!(grid_v,ind,dx,dy,LS[end].geoS,Hv) 

    init_fields_2!(phS.vD,phS.v,Hv,BC_vS,grid_v,num.vD)
    # init_fields_2!(phS.vcorrD,phS.v,HSv,BC_vS,grid_v,num.vD)

    get_height!(grid_v,ind,dx,dy,LS[end].geoL,Hv)

    init_fields_2!(phL.vD,phL.v,Hv,BC_vL,grid_v,num.vD)
    # init_fields_2!(phL.vcorrD,phL.v,HLv,BC_vL,grid_v,num.vD)


    Hsc = zeros(grid) 

    get_height!(grid,ind,dx,dy,LS[end].geoS,Hsc) #here Hsc solid

    init_fields_2!(phS.pD,phS.p,Hsc,BC_pS,grid,presintfc)

    if heat
        init_fields_2!(phS.TD,phS.T,Hsc,BC_TS,grid,θd)
    end

    #Electrolysis
    if electrolysis && num.electric_potential == 1
        init_fields_2!(phS.phi_eleD,phS.phi_ele,Hsc,BC_phi_ele,grid,phi_ele0) 
    end

    get_height!(grid,ind,dx,dy,LS[end].geoL,Hsc) #here Hsc liquid

    init_fields_2!(phL.pD,phL.p,Hsc,BC_pL,grid,presintfc)

    if heat
        init_fields_2!(phL.TD,phL.T,Hsc,BC_TL,grid,θd)
    end


    #Electrolysis
    
    if electrolysis

        printstyled(color=:green, @sprintf "\n Check %s %s %s %s %.2e %.2e %2i\n" heat heat_convection electrolysis electrolysis_convection num.τ θd nb_transported_scalars)

        for iscal=1:nb_transported_scalars
            @views phL.trans_scal[:,:,iscal] .= concentration0[iscal]
            @views init_fields_2!(phL.trans_scalD[:,iscal],phL.trans_scal[:,:,iscal],Hsc,BC_trans_scal[iscal],grid,concentration0[iscal])
        end

        if num.electric_potential == 1
            init_fields_2!(phL.phi_eleD,phL.phi_ele,Hsc,BC_phi_ele,grid,phi_ele0)
        end
    end  
    
    if electrolysis
        printstyled(color=:green, @sprintf "\n Check init_fields_2!\n")
        print_electrolysis_statistics(nb_transported_scalars,grid,phL)

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


            #TODO pre-allocate at start to save up allocations

            #Preallocate for mass flux computations
            mass_flux_vec1 = fzeros(grid)
            mass_flux_vecb = fzeros(grid)
            mass_flux_veci = fzeros(grid)
            mass_flux = zeros(grid)

           
            if advection

            #Allocations for scalar grid
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb

            Ascal = spzeros(nt, nt)
            Bscal = spzeros(nt, nt)
            rhs_scal = fnzeros(grid, num)


            all_CUTCT = zeros(grid.ny * grid.nx, nb_transported_scalars)

            us=zeros(grid)
            vs=zeros(grid)

            #save up memory: do not allocate Eus
            # Eus=zeros(grid) 
            # Evs=zeros(grid)

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


            else

            #Allocations for scalar grid
            ni = grid.nx * grid.ny
            nb = 2 * grid.nx + 2 * grid.ny
            nt = (num.nLS + 1) * ni + nb

            Ascal = spzeros(nt, nt)
            Bscal = spzeros(nt, nt)
            rhs_scal = fnzeros(grid, num)

            all_CUTCT = zeros(grid.ny * grid.nx, nb_transported_scalars)

            us=zeros(grid)
            vs=zeros(grid)


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

            end

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
            
            
            if !advection
            #call to set_heat! is there to set up the matrices for the heat equation. 
            #If the level-set is not advected, then after this call there is no need to update these matrices anymore

                _ = set_poisson(
                    BC_int, num, grid, a0_p, opC_pS, opC_uS, opC_vS,
                    AϕS, Lpm1_S, bc_Lpm1_S, bc_Lpm1_b_S, BC_pS,
                    true
                )

                set_heat!(
                    BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, true, BC_int
                )
            end

            if test_laplacian
                num.exact_laplacian = test_laplacian_pressure(num,grid_v,phL,opC_pL, Lvm1_L, bc_Lvm1_L, bc_Lvm1_b_L)
                return
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

            #call to set_heat! is there to set up the matrices for the heat equation. 
            #If the level-set is not advected, then after this call there is no need to update these matrices anymore
            set_heat!(
                BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                ATL, BTL,rhs_scal,
                opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
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

    V0S = volume(LS[end].geoS)
    V0L = volume(LS[end].geoL)


    if electrolysis

        ####################################################################################################
        #PDI (IO)
        ####################################################################################################
        
        #TODO not necessary to expose everything now for ex only LS ? and the rest later


        printstyled(color=:red, @sprintf "\n test segments\n" )


        intfc_vtx_x,intfc_vtx_y,intfc_vtx_field,intfc_vtx_connectivities,intfc_vtx_num, intfc_seg_num = convert_interfacial_D_to_segments(num,grid,phL.T,1)
        # print("\n number of interface points intfc_vtx_num ", intfc_vtx_num)
        # print("\n intfc_vtx_connectivities ",intfc_vtx_connectivities)
        # print("\n len ", size(intfc_vtx_connectivities),intfc_seg_num)

        # print("\n intfc_vtx_x ",intfc_vtx_x)
        # print("\n intfc_vtx_x ",intfc_vtx_y)


        current_t = 0.

        if num.io_pdi>0

            try
                # printstyled(color=:red, @sprintf "\n PDI test \n" )
        
                time = current_t #Cdouble
                nstep = num.current_i
           
                # phi_array=phL.phi_ele #do not transpose since python row major
                
                compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
        
                #store in us, vs instead of Eus, Evs
                interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,us,vs)

                @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                "i_current_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                "i_current_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,  
                "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                C_NULL::Ptr{Cvoid})::Cint

                interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,us,vs)
        
                # print("\n before write \n ")
        
                iLSpdi = 1 # all LS iLS = 1 # or all LS ?

                # Exposing data to PDI for IO    
                # if writing "D" array (bulk, interface, border), add "_1D" to the name
                
                printstyled(color=:magenta, @sprintf "\n PDI write_data_start_loop %.5i \n" num.current_i)

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_data_start_loop"::Cstring,
                "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
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

        ####################################################################################################

    end #if electrolysis

    if num.debug =="allocations_start"
        print("\n STOP allocations")
        return
    end


    


    


    current_t = 0.
    num.current_i =1

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
        println("\n grid.LS[1].geoL.dcap[1,1,:]",grid.LS[1].geoL.dcap[1,1,:])

        


        if !stefan
            V .= speed*ones(ny, nx)
        end

        # Solve heat equation
        if heat
            if heat_solid_phase
                kill_dead_cells!(phS.T, grid, LS[1].geoS)
                veci(phS.TD,grid,1) .= vec(phS.T)
                set_heat!(
                    BC_int[1], num, grid, opC_TS, LS[1].geoS, phS, θd, BC_TS, LS[1].MIXED, LS[1].geoS.projection,
                    ATS, BTS,rhs_scal,
                    opS, grid_u, grid_u.LS[1].geoS, grid_v, grid_v.LS[1].geoS,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTS, phS.TD, 1.0, 1.0)

                phS.TD .= ATS \ rhs_scal
                phS.T .= reshape(veci(phS.TD,grid,1), grid)
            end
            if heat_liquid_phase
                kill_dead_cells!(phL.T, grid, LS[1].geoL)
                veci(phL.TD,grid,1) .= vec(phL.T)
                set_heat!(
                    BC_int[1], num, grid, opC_TL, LS[1].geoL, phL, θd, BC_TL, LS[1].MIXED, LS[1].geoL.projection,
                    ATL, BTL,rhs_scal,
                    opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, heat_convection, advection, BC_int
                )
                mul!(rhs_scal, BTL, phL.TD, 1.0, 1.0)

                phL.TD .= ATL \ rhs_scal
                phL.T .= reshape(veci(phL.TD,grid,1), grid)

            end
        end

        ####################################################################################################
        #Electrolysis
        ####################################################################################################  
        if electrolysis
 
            if electrolysis_liquid_phase

                # printstyled(color=:green, @sprintf "\n Before transport\n")
                # print_electrolysis_statistics(nb_transported_scalars,grid,phL)

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

                
                if nb_transported_scalars>0

                    printstyled(color=:magenta, @sprintf "\n nb_transported_scalars %.5i " nb_transported_scalars)


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

                    # if imposed_velocity == "zero" && num.current_i ==16
                        #    num.ϵwall = 0.0
                        #   print("\n changed eps wall")
                    # end
                        #printstyled(color=:cyan, @sprintf "\n epsilon %.2e %.2e \n" num.ϵ num.ϵwall )

                        update_all_ls_data(num, grid, grid_u, grid_v, BC_int, periodic_x, periodic_y)
                    end

                    
                    # printstyled(color=:red, @sprintf "\n return before debug mem\n")
                    # return

                    scalar_transport!(BC_trans_scal, num, grid, opC_TL, LS[1].geoL, phL, concentration0,
                    LS[1].MIXED, LS[1].geoL.projection, opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff,Ascal,Bscal,all_CUTCT,rhs_scal)

                

                    # scalar_transport_2!(BC_trans_scal, num, grid, opC_TL, LS[1].geoL, phL, concentration0,
                    # LS[1].MIXED, LS[1].geoL.projection, opL, grid_u, grid_u.LS[1].geoL, grid_v, grid_v.LS[1].geoL,
                    # periodic_x, periodic_y, electrolysis_convection, true, BC_int, diffusion_coeff)

                    if ((num.current_i-1)%show_every == 0) 
                        printstyled(color=:cyan, @sprintf "\n after scalar transport \n")
                        print_electrolysis_statistics(nb_transported_scalars,grid,phL)
                    end

                    if imposed_velocity != "none" && num.debug== "scalar_testing"
                        scal_error=0.0
                        for iscal in 1:nb_transported_scalars

                            # print("\n maximum ",maximum(phL.trans_scalD[:,iscal]), )
                            # printstyled(color=:cyan, @sprintf "\n error after scalar transport max %.2e min %.2e\n" maximum(phL.trans_scalD[:,iscal]) minimum(phL.trans_scalD[:,iscal]))

                            scal_error_bulk = maximum(abs.(phL.trans_scal[:,:,iscal].-concentration0[iscal])./concentration0[iscal])
                            scal_error_border = maximum(abs.(vecb(phL.trans_scalD[:,iscal],grid).-concentration0[iscal])./concentration0[iscal])
                            scal_error = max(scal_error_bulk,scal_error_border,scal_error)

                        end

                        printstyled(color=:cyan, @sprintf "\n error after scalar transport %.2e CFL %.2e\n" scal_error num.v_inlet*num.dt0/grid.dx[1,1])

                        # printstyled(color=:red, @sprintf "\n Poiseuille \n")

                        # # Check the velocity field before the scalar transport
                        # test_Poiseuille(num,phL.vD,grid_v)

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

                        # printstyled(color=:cyan, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

                        # compute_grad_p!(num,grid, grid_u, grid_v, phL.pD, opC_pL, opC_uL, opC_vL)

                    end

                end #nb_transported_scalars>0

                # if imposed_velocity =="none"
                #     printstyled(color=:red, @sprintf "\n after scalar transport \n")

                #     # Check the velocity field before the scalar transport
                #     test_Poiseuille(num,phL,grid_v)
                    
                # end
                
    
                #Electrolysis: Poisson  
                if num.electric_potential == 1
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

                    if electrolysis && nb_transported_scalars>1 && num.bulk_conductivity != 0
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
                        # TODO -(-i/kappa) in Flower ? so i_butler not -i_butler
                        # For small cells
                        if num.bulk_conductivity == 0
                            BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)
                        elseif num.bulk_conductivity == 1
                            # Recommended as long as cell merging not implemented:
                            # Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
                            # which produces higher concentration in front of the contact line
                            BC_phi_ele.left.val .= i_butler./elec_cond[:,1]
                        elseif num.bulk_conductivity == 2
                            BC_phi_ele.left.val .= i_butler./vecb_L(elec_condD, grid)

                            iLS = 1 #TODO end ? if several LS ?
                            for j in 1:grid.ny
                                II = CartesianIndex(j,1)
                                if grid.LS[iLS].geoL.cap[II,5] < ϵ
                                    BC_phi_ele.left.val[j] = i_butler[j]/elec_cond[j,1] 
                                end
                            end
                        end

                        # print("\n before ",BC_phi_ele.left.val) 

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

                        # TODO
                        #Remove Nan when dividing by conductivity which may be null
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

                    #    print(BC_phi_ele.left.val) 

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
                     set_poisson_variable_coeff!(
                        [BC_phi_ele.int], num, grid, grid_u, grid_v, a0_p, opC_pL, opC_uL, opC_vL,
                        Aphi_eleL, 
                        # elec_Lpm1_L, elec_bc_Lpm1_L, elec_bc_Lpm1_b_L, 
                        BC_phi_ele,rhs_scal,
                        true,elec_condD
                    )

                    # print("\n Aphi_eleL: ",any(isnan, Aphi_eleL),"\n rhs_scal: ",any(isnan, rhs_scal),"\n")

                    # print("\n veci rhs_scal 2: ",any(isnan, veci(rhs_scal,grid,1)),"\n veci rhs_scal 1 : ",any(isnan, veci(rhs_scal,grid,1)),"\n")
            

                    # print("\n \n elec_condD",vecb_L(elec_condD, grid))
                    # print("\n \n BC_phi_ele.left.val",BC_phi_ele.left.val)

                    

                    # print("\n \n vecb_L",vecb_L(rhs_scal, grid))

                    # print("\n \n vecb_R",vecb_R(rhs_scal, grid))

                    # print("\n \n vecb_T",vecb_T(rhs_scal, grid))

                    # print("\n \n vecb_B",vecb_B(rhs_scal, grid))

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


                    veci(rhs_scal,grid,1) .+= op.opC_pL.M * vec(b_phi_ele)
                
                    res_phi_ele = zeros(size(rhs_scal))
                
                    @time @inbounds @threads for i in 1:Aphi_eleL.m
                        @inbounds Aphi_eleL[i,i] += 1e-10
                    end
                    
                    # @time res_phi_ele .= Aphi_eleL \ rhs_scal

                    @time res_phi_ele .= Aphi_eleL \ rhs_scal

                    #TODO or use mul!(rhs_scal, BTL, phL.TD, 1.0, 1.0) like in :

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


                    # phL.phi_eleD .= Aphi_eleL \ rhs_scal
                    phL.phi_eleD .= res_phi_ele

                    phL.phi_ele .= reshape(veci(phL.phi_eleD,grid,1), grid)


                    if any(isnan, phL.phi_eleD)
                        print("\n phL.uD: ",any(isnan, phL.uD) , "\n phL.vD: ",any(isnan, phL.vD) , "\n phL.TD: ",any(isnan, phL.TD) , "\n phS.uD: ",any(isnan, phS.uD) , "\n phS.vD: ",any(isnan, phS.vD) , "\n phS.TD: ",any(isnan, phS.TD) ,
                        "\n phL.trans_scalD: ",any(isnan, phL.trans_scalD) , "\n phL.phi_eleD: ",any(isnan, phL.phi_eleD) ,
                        "\n phL.u: ",norm(phL.u) > 1e8 , "\n phS.u: ",norm(phS.u) > 1e8 , "\n phL.T: ",norm(phL.T) > 1e8 , "\n phS.T: ",norm(phS.T) > 1e8 , "\n phL.trans_scal: ",norm(phL.trans_scal) > 1e8 , "\n phL.phi_ele: ",norm(phL.phi_ele) > 1e8)
            
                        print("\n phL.phi_eleD: ",any(isnan, phL.phi_eleD),"\n phL.phi_ele: ",any(isnan, phL.phi_ele),"\n")

                        print("\n Aphi_eleL: ",any(isnan, Aphi_eleL),"\n rhs_scal: ",any(isnan, rhs_scal),"\n")

                        print("\n \n vecb_L",vecb_L(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_R",vecb_R(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_T",vecb_T(phL.phi_eleD[:,1], grid))

                        # print("\n \n vecb_B",vecb_B(phL.phi_eleD[:,1], grid))
                
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

                        # #############################################################################################################
                        # #Print values on line
                        # id_x = 1
                        # for iplot in 1:ny
                        #     printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n" iplot grid.y[iplot]/num.plot_xscale phL.saved_scal[iplot,id_x,2] phL.saved_scal[iplot,id_x,3] phL.saved_scal[iplot,id_x,4] phL.trans_scal[iplot,id_x,1] grid.LS[1].u[iplot,id_x] phL.trans_scalD[iplot,id_x,1])
                        # end
                        # #############################################################################################################


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

                end #electric_potential

            end
        end
        



        
        #PDI (IO)      
        if num.io_pdi>0

            try
                # printstyled(color=:red, @sprintf "\n PDI test \n" )
        
                time = current_t #Cdouble
                nstep = num.current_i
            
                # phi_array=phL.phi_ele #do not transpose since python row major
                
                # Compute electrical current, interpolate velocity on scalar grid
                compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
        
                #store in us, vs instead of Eus, Evs
                interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,us,vs)

                @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                "i_current_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                "i_current_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,  
                "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                C_NULL::Ptr{Cvoid})::Cint
            
                interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,us,vs)
                    
                iLSpdi = 1 # TODO all LS

                # Exposing data to PDI for IO    
                # if writing "D" array (bulk, interface, border), add "_1D" to the name

                PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
                C_NULL::Ptr{Cvoid})::Cint
        
                
                #TODO debug with volume fraction

                # A = zeros(gv.ny, gv.nx+1)
                # for jplot in 1:gv.ny
                #     for iplot in 1:gv.nx+1
                #     II = CartesianIndex(jplot, iplot) #(id_y, id_x)
                #     pII = lexicographic(II, gp.ny + 1)
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

        ####################################################################################################

        print("\n test",electrolysis_phase_change_case, occursin("levelset",electrolysis_phase_change_case))

        if electrolysis && electrolysis_phase_change_case != "None"
            printstyled(color=:magenta, @sprintf "\n compute_mass_flux!\n")

            compute_mass_flux!(num,grid,phL,opC_pL,1,mass_flux_vec1,mass_flux_vecb,mass_flux_veci,mass_flux)

        end
        
        #    grid.LS[i].α  which is the angle of the outward point normal with respect to the horizontal axis

        for iLS in 1:nLS
            if is_stefan(BC_int[iLS])
                update_stefan_velocity(num, grid, iLS, LS[iLS].u, phS.T, phL.T, periodic_x, periodic_y, λ, Vmean)
            elseif is_fs(BC_int[iLS]) || occursin("levelset",electrolysis_phase_change_case)
                grid.V .= 0.0
                printstyled(color=:green, @sprintf "\n V %.2e max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(V)) maximum(abs.(phL.u)) maximum(abs.(phL.v)))

                if electrolysis_phase_change_case!="none"    
                    if occursin("levelset",electrolysis_phase_change_case)

                        if electrolysis_phase_change_case == "levelset"

                            printstyled(color=:magenta, @sprintf "\n phase-change velocity")

                            #return


                            # plot_electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)

                            # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 

                            varnH2 = -sum(mass_flux) * diffusion_coeff[1] 

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

                            nH2 = new_nH2

                            update_free_surface_velocity_electrolysis(num, grid, grid_u, grid_v, iLS, phL.uD, phL.vD, periodic_x, periodic_y, Vmean, phL.trans_scalD[:,1],diffusion_coeff[1],concentration0[1],electrolysis_phase_change_case,mass_flux)

                        elseif electrolysis_phase_change_case == "levelset_averaged"

                            printstyled(color=:magenta, @sprintf "\n phase-change velocity")

                            #return


                            # plot_electrolysis_velocity!(num, grid, LS, V, TL, MIXED, periodic_x, periodic_y, concentration_scal_intfc)

                            # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 

                            varnH2 = -sum(mass_flux) * diffusion_coeff[1] 

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

                            nH2 = new_nH2

                            iLS = 1
                            #TODO check velocity
                            @inbounds @threads for II in grid.LS[iLS].MIXED
                                grid.V[II] = -sum(mass_flux) * diffusion_coeff[1] *(1.0/num.rho2-1.0/num.rho1).*diffusion_coeff[1].*num.MWH2
                            end

                        end

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

                # print("\n test left H2", vecb_L(phL.trans_scalD[:,1], grid))

               
                # print("\n testvalflux1")
                # for jplot in 1:ny
                #     for iplot in 1:nx
                #         if phL.saved_scal[jplot,iplot,1] !=0.0
                #             printstyled(color=:green, @sprintf "\n j: %5i %5i %.2e\n" iplot jplot phL.saved_scal[jplot,iplot,1])
                #         end
                #     end
                # end





                # printstyled(color=:cyan, @sprintf "\n div(0,grad): %.2e %.2e %.2e %.2e \n" sum(phL.saved_scal[:,:,1]) sum(phL.saved_scal[:,:,2]) sum(phL.saved_scal[:,:,3]) sum(phL.saved_scal[:,:,4]))

                
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


                # mass_fluxO=compute_mass_flux!(num,grid, grid_u, grid_v, phL, phS,  opC_pL, opC_pS,diffusion_coeff,3)

                # Print values on line
                # for iplot in 1:ny
                #     II = CartesianIndex(iplot, 1) #(id_y, id_x)
                #     pII = lexicographic(II, grid.ny)
                #     printstyled(color=:green, @sprintf "\n j: %5i %.2e %.2e %.2e %.2e \n" iplot grid.y[iplot]/num.plot_xscale mass_flux_vec1[pII] mass_flux_vecb[pII] mass_flux_veci[pII])
                # end
       
                #TODO mass_fluxL
                # mass_fluxL=-... * diffusion_coeff[1] 



                # Minus sign because normal points toward bubble and varnH2 for gaz, not liquid phase 
                varnH2 = -sum(mass_flux) * diffusion_coeff[1] 

                #TODO mode_2d==0 flux corresponds to cylinder of length 1
                #2D cylinder reference length
                if mode_2d == 1
                    varnH2 .*= num.ref_thickness_2d
                end

                #Pliquid is the average value of p over the bubble interface plus the ambient operating pressure (P).
                p_liq= num.pres0 + mean(veci(phL.pD,grid,2)) #TODO here one bubble
                # p_g=p_liq + 2 * num.σ / current_radius #3D
                p_g=p_liq + num.σ / current_radius #2D

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
                    return
                end

               
                printstyled(color=:cyan, @sprintf "\n div(0,grad): %.5i %.2e %.2e %.2e %.2e\n" nx num.τ num.L0/nx (current_radius-previous_radius)/(num.L0/nx) sum(phL.saved_scal[:,:,1]))
                
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

            elseif (electrolysis && electrolysis_phase_change_case == "imposed_radius4")

                #CFL 0.5
                current_radius = current_radius + grid.dx[1,1]/4

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
                elseif is_fs(bc) || occursin("levelset",electrolysis_phase_change_case)
                    rhs_LS .= 0.0
                    LS[iLS].A.nzval .= 0.0
                    LS[iLS].B.nzval .= 0.0
                    IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, θ_out, num.τ, periodic_x, periodic_y)
                    BC_LS_interior!(num, grid, grid_u, grid_v, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
                    BC_LS!(grid, LS[iLS].u, LS[iLS].A, LS[iLS].B, rhs_LS, BC_u)
                    utmp .= reshape(gmres(LS[iLS].A, LS[iLS].B * vec(LS[iLS].u) .+ rhs_LS), grid)

                    rhs_LS .= 0.0
                    S2IIOE!(grid, grid_u, grid_v, LS[iLS].A, LS[iLS].B, utmp, LS[iLS].u, θ_out, num.τ, periodic_x, periodic_y)
                    BC_LS_interior!(num, grid, grid_u, grid_v, iLS, LS[iLS].A, LS[iLS].B, rhs_LS, BC_int, periodic_x, periodic_y)
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
                if auto_reinit == 1 && (num.current_i-1)%num.reinit_every == 0
                    for iLS in 1:nLS
                        if !is_wall(BC_int[iLS])
                            ls_rg, rl_rg_v = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)
                            println("$(ls_rg)")
                            printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
                            if ls_rg >= δreinit || num.current_i == 1
                                print("(ls_rg >= δreinit || num.current_i == 1): yes")
                                # println("yes")
                                RK2_reinit!(ls_scheme, grid, ind, iLS, LS[iLS].u, nb_reinit, periodic_x, periodic_y, BC_u, BC_int)
                                
                                ls_rg, rl_rg_v = rg(num, grid, LS[iLS].u, periodic_x, periodic_y, BC_int)
                                println("$(ls_rg) ")
                                printstyled(color=:green, @sprintf "\n ls_rg : %.2e \n" ls_rg)
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
            if free_surface && breakup ==1
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
                return
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
                    # printstyled(color=:red, @sprintf "\n PDI test \n" )
            
                    time = current_t #Cdouble
                    nstep = num.current_i
               
                    # phi_array=phL.phi_ele #do not transpose since python row major
                    
                    # Compute electrical current, interpolate velocity on scalar grid
                    compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current
            

                    #store in us, vs instead of Eus, Evs
                    interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,us,vs)

                    #TODO i_current_mag need cond

                    @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
                    "i_current_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "i_current_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,  
                    "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    C_NULL::Ptr{Cvoid})::Cint

                    # interpolate_grid_liquid!(grid,grid_u,grid_v,phL.Eu, phL.Ev,Eus,Evs)
            
                    interpolate_grid_liquid!(grid,grid_u,grid_v,phL.u,phL.v,us,vs)
                        
                    iLSpdi = 1 # TODO all LS

                    # Exposing data to PDI for IO    
                    # if writing "D" array (bulk, interface, border), add "_1D" to the name

                    PDI_status = @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
                    "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
                    "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
                    "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_p"::Cstring, LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_u"::Cstring, grid_u.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    "levelset_v"::Cstring, grid_v.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
                    "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2"::Cstring, phL.trans_scalD[:,1]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_KOH"::Cstring, phL.trans_scalD[:,2]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "trans_scal_1D_H2O"::Cstring, phL.trans_scalD[:,3]::Ptr{Cdouble}, PDI_OUT::Cint,
                    # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
                    # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
                    "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
                    "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
                    "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint, 
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
    

                print_electrolysis_statistics(nb_transported_scalars,grid,phL)

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
        #TODO save radius
        return # radius
    # elseif flapping
    #     return xc, yc
    else
        return
    end
end

