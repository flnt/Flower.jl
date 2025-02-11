using Revise
using Flower


# PDI

localARGS = ARGS
@show localARGS

# print("\n Arguments ", localARGS)
# print("\n length(localARGS) ",length(localARGS))

if length(localARGS)>0
    yamlfile = localARGS[1]
    printstyled(color=:magenta, @sprintf "\n YAML file ")
    print(yamlfile)
    print("\n")
    yamlpath = yamlfile
    yamlnamelist = split(yamlfile, "/")
    yamlname = yamlnamelist[end]
end

data = YAML.load_file(yamlpath)

# Dictionaries 
prop_dict = PropertyDict(data)
io = PropertyDict(prop_dict.plot) 
flower = PropertyDict(prop_dict.flower)
mesh = PropertyDict(flower.mesh)
sim = PropertyDict(flower.simulation)
phys = PropertyDict(flower.physics)
macros = PropertyDict(flower.macros) #to parse code from .yml

study = PropertyDict(prop_dict.study)

# print parameters by evaluating Julia code stored in .yml   
eval(Meta.parseall(macros.print_parameters))



npts = study.meshes
n_cases = length(npts)
print("\n number of points ", npts, "\n")



if io.pdi>0

    @debug "Before PDI init"
    yml_file = yamlfile
    conf = @ccall "libparaconf".PC_parse_path(yml_file::Cstring)::PC_tree_t
    @debug "after conf"
    getsubyml = @ccall "libparaconf".PC_get(conf::PC_tree_t,".pdi"::Cstring)::PC_tree_t  
    @debug "after getsubyml"
    local pdi_status = @ccall "libpdi".PDI_init(getsubyml::PC_tree_t)::Cint
    @debug "after PDI_init"

    # Send meta-data to PDI
    mpi_coords_x = 1
    mpi_coords_y = 1
    mpi_max_coords_x = 1
    mpi_max_coords_y = 1
    local nx = 32
    local ny = 32
    nstep = 0
    # nx=gp.nx
    # ny=gp.ny

    #TODO check Clonglong ...

    phys_time = 0.0 #Cdouble
    # nstep = num.current_i
    

    local PDI_status = @ccall "libpdi".PDI_multi_expose("init_PDI"::Cstring, 
            "mpi_coords_x"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
            "mpi_coords_y"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
            "mpi_max_coords_x"::Cstring, mpi_max_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
            "mpi_max_coords_y"::Cstring, mpi_max_coords_y::Ref{Clonglong}, PDI_OUT::Cint,
            "nx"::Cstring, nx::Ref{Clonglong}, PDI_OUT::Cint,
            "ny"::Cstring, ny::Ref{Clonglong}, PDI_OUT::Cint,
            "nb_transported_scalars"::Cstring, phys.nb_transported_scalars::Ref{Clonglong}, PDI_OUT::Cint,
            "nb_levelsets"::Cstring, phys.nb_levelsets::Ref{Clonglong}, PDI_OUT::Cint,
            "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

    @debug "after PDI_multi_expose"

    @debug "After full PDI init"

end #if io.pdi>0

# arrays to store errors
l1 = zeros(n_cases)
l2 = zeros(n_cases)
loo = zeros(n_cases)
l1_mixed = zeros(n_cases)
l2_mixed = zeros(n_cases)
loo_mixed = zeros(n_cases)
l1_full = zeros(n_cases)
l2_full = zeros(n_cases)
loo_full = zeros(n_cases)

cell_volume_list = zeros(n_cases)


# Convergence study loop
for (i,n) in enumerate(npts)


    # init regular grid
    scalar_mesh_x = collect(LinRange(mesh.xmin, mesh.xmax, n + 1))    
    scalar_mesh_y = collect(LinRange(mesh.ymin, mesh.ymax, n + 1))

    @debug "Before Numerical"
    global num = Numerical(
        CFL = sim.CFL,
        Re = Re,
        TEND=phys.end_time,
        x = scalar_mesh_x,
        y = scalar_mesh_y,
        xcoord = phys.intfc_x,
        ycoord = phys.intfc_y,
        case = sim.case,
        R = phys.radius,
        max_iterations = sim.max_iter,
        save_every = sim.max_iter,
        ϵ = sim.epsilon, 
        ϵwall = sim.epsilon_wall,
        epsilon_mode = sim.epsilon_mode,
        nLS = phys.nb_levelsets,
        nb_transported_scalars=phys.nb_transported_scalars,
        concentration0=phys.concentration0, 
        epsilon_concentration=phys.epsilon_concentration,
        diffusion_coeff=phys.diffusion_coeff,
        temperature0=phys.temperature0,
        i0=phys.i0,
        phi_ele0=phys.phi_ele0,
        phi_ele1=phys.phi_ele1,
        alpha_c=phys.alpha_c,
        alpha_a=phys.alpha_a,
        Ru=phys.Ru,
        Faraday=phys.Faraday,
        MWH2=phys.MWH2,
        θd=phys.temperature0,
        eps=sim.eps,
        mu1=phys.mu1,
        mu2=phys.mu2,
        rho1=phys.rho1,
        rho2=phys.rho2,
        u_inf = 0.0,
        v_inf = 0.0,
        pres0=phys.pres0,
        g = phys.g,
        β = phys.beta,
        σ = phys.sigma,   
        reinit_every = sim.reinit_every,
        nb_reinit = sim.nb_reinit,
        δreinit = sim.delta_reinit,
        n_ext_cl = sim.n_ext,
        NB = sim.NB,
        plot_xscale = io.scale_x,
        dt0 = sim.dt0,
        concentration_check_factor = sim.concentration_check_factor,
        radial_vel_factor = phys.radial_vel_factor,
        debug = sim.debug,
        v_inlet = phys.v_inlet,
        prediction = sim.prediction,
        null_space = sim.null_space,
        io_pdi = io.pdi,
        bulk_conductivity = sim.bulk_conductivity,
        electrical_potential = sim.electrical_potential,
        contact_angle = sim.contact_angle,
        convection_Cdivu = sim.convection_Cdivu,
        convection_mode = sim.convection_mode,
        advection_LS_mode = sim.advection_LS_mode,
        scalar_bc = sim.scalar_bc,
        scalar_scheme = sim.scalar_scheme,
        solver = sim.solver,
        mass_flux = sim.mass_flux,
        average_liquid_solid = sim.average_liquid_solid,
        index_phase_change = sim.index_phase_change,
        index_electrolyte = sim.index_electrolyte,
        extend_field = sim.extend_field,
        average_velocity = sim.average_velocity,
        laplacian = sim.laplacian,
        electrical_potential_max_iter = sim.electrical_potential_max_iter,
        electrical_potential_residual = sim.electrical_potential_residual,
        )
    Broadcast.broadcastable(num::Numerical) = Ref(num) #do not broadcast num 
    @debug "After Numerical"

    global gp, gu, gv = init_meshes(num)
    #gp, gu, gv = init_meshes(num) does not work for eval(Meta.parseall(macros.boundaries))
    global op, phS, phL = init_fields(num, gp, gu, gv)

    gp.LS[1].u .= 1.0 #deactivate interface

    # Define boundary conditions
    eval(Meta.parseall(macros.boundaries))


    if num.io_pdi>0


        # Send meta-data to PDI
        nx=gp.nx
        ny=gp.ny
        
        try
            local PDI_status = @ccall "libpdi".PDI_multi_expose("init_PDI"::Cstring, 
                    "mpi_coords_x"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                    "mpi_coords_y"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                    "mpi_max_coords_x"::Cstring, mpi_max_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                    "mpi_max_coords_y"::Cstring, mpi_max_coords_y::Ref{Clonglong}, PDI_OUT::Cint,
                    "nx"::Cstring, nx::Ref{Clonglong}, PDI_OUT::Cint,
                    "ny"::Cstring, ny::Ref{Clonglong}, PDI_OUT::Cint,
                    "nb_transported_scalars"::Cstring, phys.nb_transported_scalars::Ref{Clonglong}, PDI_OUT::Cint,
                    "nb_levelsets"::Cstring, phys.nb_levelsets::Ref{Clonglong}, PDI_OUT::Cint,
                    "nstep"::Cstring, num.current_i::Ref{Clonglong}, PDI_OUT::Cint,
                    C_NULL::Ptr{Cvoid})::Cint

        catch
            print("\n Bug init_PDI \n")

            print("\n PDI_status ",PDI_status, "\n")

        end

        @debug "after PDI_multi_expose"

        @debug "After full PDI init"

    end #if num.io_pdi>0


    # Define interfaces (for bubbles, drops...)
    if sim.activate_interface == 1

        gp.LS[1].u .= sqrt.((gp.x .- phys.intfc_x).^2 + (gp.y .- phys.intfc_y).^2) - phys.radius * ones(gp)
    
        #modify velocity field near interface
        su = sqrt.((gv.x .- phys.intfc_x).^2 .+ (gv.y .- phys.intfc_y).^2)
        R1 = phys.radius + 3.0*num.Δ

        bl = 4.0
        for II in gv.ind.all_indices
            if su[II] <= R1
                phL.v[II] = 0.0
            # elseif su[II] > R1
            #     uL[II] = tanh(bl*(su[II]-R1))
            end
        end

    elseif sim.activate_interface == -1
        gp.LS[1].u .= sqrt.((gp.x .- phys.intfc_x).^2 + (gp.y .- phys.intfc_y).^2) - phys.radius * ones(gp)
        gp.LS[1].u .*= -1.0

    else
        gp.LS[1].u .= 1.0
    end

    test_LS(gp)

    # Create segments from interface
    # x,y,field,connectivities,num_vtx = convert_interfacial_D_to_segments(num,gp,phL.T,1)
    # print("\n number of interface points ", num_vtx)
    # # print("\n x",x)
    # # print("\n x",y)
    # # print("\n x",field)
    # print("\n x",connectivities)
    # print("\n x",num_vtx)

    printstyled(color=:green, @sprintf "\n sim.CFL : %.2e dt : %.2e\n" sim.CFL sim.CFL*phys.ref_length/mesh.nx/phys.v_inlet)

    # printstyled(color=:green, @sprintf "\n Initialisation0 \n")
    # print_electrolysis_statistics(num,gp,phL)

    #init Bulk
    phL.T .= phys.temperature0
    phS.T .= phys.temperature0

    vPoiseuille = Poiseuille_fmax.(gv.x,phys.v_inlet,phys.ref_length) 
    vPoiseuilleb = Poiseuille_fmax.(gv.x[1,:],phys.v_inlet,phys.ref_length) 

    phL.u .= 0.0
    phL.v .= vPoiseuille 

    vecb_B(phL.vD,gv) .= vPoiseuilleb

    for iscal=1:phys.nb_transported_scalars
        phL.trans_scal[:,:,iscal] .= phys.concentration0[iscal]
    end

    phL.phi_ele .= phys.phi_ele0

    printstyled(color=:green, @sprintf "\n Initialisation \n")


    printstyled(color=:green, @sprintf "\n TODO timestep sim.CFL scal, and print \n")



    if sim.time_scheme == "FE"
        time_scheme = FE
    else
        time_scheme = CN
    end


    if num.io_pdi>0

        try
            # printstyled(color=:red, @sprintf "\n before pdi \n")

            # printstyled(color=:red, @sprintf "\n PDI test \n" )


            # if num.solve_electrical_potential>0

            # phi_array=phL.phi_ele #do not transpose since python row major
            
            # compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS) #TODO current

            # Eus,Evs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.Eu, phL.Ev)

            # us,vs = interpolate_grid_liquid(grid,grid_u,grid_v,phL.u,phL.v)

            # print("\n before write \n ")


            #TODO "levelset_p_tot"::Cstring, gp.LS[2].u::Ptr{Cdouble}, PDI_OUT::Cint,


            iLSpdi = 1 # all LS iLS = 1 # or all LS ?

            # Exposing data to PDI for IO    
            # if writing "D" array (bulk, interface, border), add "_1D" to the name
            
            printstyled(color=:magenta, @sprintf "\n PDI write_data_start_loop %.5i \n" num.current_i)

            #print("\n size LS wall ", size( gp.LS[2].u))
            LStable = zeros(gp)
            if phys.nb_levelsets>1
                LStable = gp.LS[2].u
            end


            # Compute electrical current, interpolate velocity on scalar grid
            # compute_grad_phi_ele!(num, gp, gu, gv, phL, phS, op.opC_pL, op.opC_pS) #TODO current

            # compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS, elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p0,tmp_vec_p1) #TODO current

            
            us=zeros(gp) #TODO allocate only once
            vs=zeros(gp)

            # #store in us, vs instead of Eus, Evs
            # interpolate_grid_liquid!(gp,gu,gv,phL.Eu, phL.Ev,us,vs)

            # #TODO i_current_mag need cond

            # @ccall "libpdi".PDI_multi_expose("write_data_elec"::Cstring,
            # "i_current_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,  
            # "i_current_mag"::Cstring, phL.i_current_mag::Ptr{Cdouble}, PDI_OUT::Cint,
            # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
            # C_NULL::Ptr{Cvoid})::Cint


            interpolate_grid_liquid!(gp,gu,gv,phL.u,phL.v,us,vs)

            current_radius = phys.radius

            PDI_status = @ccall "libpdi".PDI_multi_expose("write_initialization"::Cstring,
            "nstep"::Cstring, num.current_i::Ref{Clonglong}, PDI_OUT::Cint,
            "time"::Cstring, phys_time::Ref{Cdouble}, PDI_OUT::Cint,
            "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
            "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
            "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_p"::Cstring, gp.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_u"::Cstring, gu.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_v"::Cstring, gv.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
            "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
            # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
            "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
            # "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
            "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
            "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,      
            "radius"::Cstring, current_radius::Ref{Cdouble}, PDI_OUT::Cint,  
            # "intfc_vtx_num"::Cstring, intfc_vtx_num::Ref{Clonglong}, PDI_OUT::Cint, 
            # "intfc_seg_num"::Cstring, intfc_seg_num::Ref{Clonglong}, PDI_OUT::Cint, 
            # "intfc_vtx_x"::Cstring, intfc_vtx_x::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_y"::Cstring, intfc_vtx_y::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_field"::Cstring, intfc_vtx_field::Ptr{Cdouble}, PDI_OUT::Cint,
            # "intfc_vtx_connectivities"::Cstring, intfc_vtx_connectivities::Ptr{Clonglong}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

        catch error
            printstyled(color=:red, @sprintf "\n PDI error \n")
            print(error)
            printstyled(color=:red, @sprintf "\n PDI error \n")
        end

    end #if io_pdi

    # if num.io_pdi>0
    #     iLSpdi = 1 # TODO all grid.LS                
    #     PDI_status = @ccall "libpdi".PDI_multi_expose("write_capacities"::Cstring,                    
    #     "dcap"::Cstring, permutedims(gp.LS[iLSpdi].geoL.dcap, (3, 1, 2))::Ptr{Cdouble}, PDI_OUT::Cint,                            
    #     C_NULL::Ptr{Cvoid})::Cint 
    #     # try                            
    #     #     iLSpdi = 1 # TODO all grid.LS                
    #     #     PDI_status = @ccall "libpdi".PDI_multi_expose("write_capacities"::Cstring,                    
    #     #     "dcap"::Cstring, gp.LS[iLSpdi].geoL.dcap'::Ptr{Cdouble}, PDI_OUT::Cint,                            
    #     #     C_NULL::Ptr{Cvoid})::Cint                           
    #     # catch error
    #     #     printstyled(color=:red, @sprintf "\n PDI error \n")
    #     #     print(error)
    #     #     printstyled(color=:red, @sprintf "\n PDI error \n")
    #     # end
    # end #if io_pdi
    # printstyled(color=:red, @sprintf "\n after pdi \n")

    # printstyled(color=:red, @sprintf "\n before run_forward \n")

    run_forward!(
        num, gp, gu, gv, op, phS, phL;
        periodic_x = (sim.periodic_x == 1),
        periodic_y = (sim.periodic_y == 1),
        BC_uL = BC_uL,
        BC_uS=BC_uS,
        BC_vL = BC_vL,
        BC_vS=BC_vS,
        BC_pL = BC_pL,
        BC_pS=BC_pS,
        BC_u = BC_u,
        BC_int = BC_int,
        BC_trans_scal=BC_trans_scal,
        BC_phi_ele = BC_phi_ele,
        auto_reinit = sim.auto_reinit,
        time_scheme = time_scheme,
        electrolysis = true,
        navier_stokes = true,
        ns_advection = (sim.ns_advection ==1),
        ns_liquid_phase = (sim.solve_Navier_Stokes_liquid_phase == 1),
        verbose = true,
        show_every = sim.show_every,
        electrolysis_convection = (sim.electrolysis_convection ==1),  
        electrolysis_liquid_phase = true,
        electrolysis_phase_change_case = sim.electrolysis_phase_change_case,
        electrolysis_reaction = phys.electrolysis_reaction, 
        imposed_velocity = sim.imposed_velocity,
        adapt_timestep_mode = sim.adapt_timestep_mode,#1,
        non_dimensionalize=sim.non_dimensionalize,
        mode_2d = sim.mode_2d,
        breakup = sim.breakup,    
    )

    @debug "After run"


    # #cut small cells for error
    # number_small_cells_for_error = 0
    # cutoff_for_error_volume = 1e-12 #TODO

    # for II in gp.ind.all_indices
    #     if gp.LS[1].geoL.cap[II,5] < cutoff_for_error_volume
    #         number_small_cells_for_error += 1
    #         Tana[II] = 0.0
    #         T[II] = 0.0
    #     end
    # end

    # printstyled(color=:green, @sprintf "\n number_small_cells_for_error %.3i \n" number_small_cells_for_error)

    LIQUID = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .> (1-1e-16)]
    MIXED = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .<= (1-1e-16) .&& gp.LS[1].geoL.cap[:,:,5] .> 1e-16]


    norm_all = relative_errors(phL.v, vPoiseuille, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
    norm_mixed = relative_errors(phL.v, vPoiseuille, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
    norm_full = relative_errors(phL.v, vPoiseuille, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)

    l1[i] = norm_all[1]
    l2[i] = norm_all[2]
    loo[i] = norm_all[3]

    l1_mixed[i] = norm_mixed[1]
    l2_mixed[i] = norm_mixed[2]
    loo_mixed[i] = norm_mixed[3]

    l1_full[i] = norm_full[1]
    l2_full[i] = norm_full[2]
    loo_full[i] = norm_full[3]

    cell_volume_list[i] = minimum(gp.LS[1].geoL.dcap[:,:,5])

    print("\n analytical ",-0.011655612832847977)


end #convergence

min_cell_volume = minimum(gp.LS[1].geoL.cap[:,:,5])

print("\n min_cell_volume ",min_cell_volume," type ",typeof(min_cell_volume))

local PDI_status = @ccall "libpdi".PDI_multi_expose("convergence_study"::Cstring, 
"n_tests"::Cstring, n_cases::Ref{Clonglong}, PDI_OUT::Cint,
"nx_list"::Cstring, npts::Ptr{Clonglong}, PDI_OUT::Cint,
"cell_volume_list"::Cstring, cell_volume_list::Ptr{Cdouble}, PDI_OUT::Cint,
"l1_rel_error"::Cstring, l1::Ptr{Cdouble}, PDI_OUT::Cint,
"l2_rel_error"::Cstring, l2::Ptr{Cdouble}, PDI_OUT::Cint,
"linfty_rel_error"::Cstring, loo::Ptr{Cdouble}, PDI_OUT::Cint,
"l1_rel_error_full_cells"::Cstring, l1_full::Ptr{Cdouble}, PDI_OUT::Cint,
"l2_rel_error_full_cells"::Cstring, l2_full::Ptr{Cdouble}, PDI_OUT::Cint,
"linfty_rel_error_full_cells"::Cstring, loo_full::Ptr{Cdouble}, PDI_OUT::Cint,
"l1_rel_error_partial_cells"::Cstring, l1_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
"l2_rel_error_partial_cells"::Cstring, l2_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
"linfty_rel_error_partial_cells"::Cstring, loo_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
"domain_length"::Cstring, L0::Ref{Cdouble}, PDI_OUT::Cint,
"min_cell_volume"::Cstring, min_cell_volume::Ref{Cdouble}, PDI_OUT::Cint,
C_NULL::Ptr{Cvoid})::Cint
 
if io.pdi>0
    try
        
        # local PDI_status = @ccall "libpdi".PDI_event("close_pycall"::Cstring)::Cint

        # PDI_event("finalization");
        # local PDI_status = @ccall "libpdi".PDI_event("finalization"::Cstring)::Cint


        local PDI_status = @ccall "libpdi".PDI_finalize()::Cint
        # printstyled(color=:red, @sprintf "\n PDI end\n" )

    catch error
        printstyled(color=:red, @sprintf "\n PDI error \n")
        print(error)
    end
end #if io.pdi>0

printstyled(color=:red, @sprintf "\n After PDI \n")

#Tests 
eval(Meta.parseall(macros.test_end))
