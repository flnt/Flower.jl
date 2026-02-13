# **Author**: Jean-Michel Batto (CEA, Maison de la Simulation)

# using Revise 
using Flower
using YAML
using Printf

# ---------------------------------------------------------------------------------------------
# FUNCTION BARRIER: Wrapper pour exécuter la simulation
# Cela permet d'isoler le code compilé des variables dynamiques (BC) générées par eval()
# ---------------------------------------------------------------------------------------------

function execute_simulation_step(num, gp, gu, gv, op, phS, phL, sim, phys, time_scheme, bc_dict)
    
    # Préparation des arguments optionnels (Keyword Arguments)
    # On construit un NamedTuple uniquement avec les BC qui existent vraiment
    kwargs_bc = Dict{Symbol, Any}()
    
    # Liste des clés possibles
    keys_to_check = [:BC_uL, :BC_uS, :BC_vL, :BC_vS, :BC_pL, :BC_pS, :BC_u, :BC_int, :BC_trans_scal, :BC_phi_ele]
    
    for k in keys_to_check
        if haskey(bc_dict, k)
            kwargs_bc[k] = bc_dict[k]
        end
    end

    # Appel avec le splatting operator (...) qui déballe le dictionnaire en arguments
    run_forward!(
        num, gp, gu, gv, op, phS, phL;
        periodic_x = (sim.periodic_x == 1),
        periodic_y = (sim.periodic_y == 1),
        auto_reinit = sim.auto_reinit,
        time_scheme = time_scheme,
        electrolysis = true,
        navier_stokes = true,
        ns_advection = (sim.ns_advection == 1),
        ns_liquid_phase = (sim.solve_Navier_Stokes_liquid_phase == 1),
        verbose = true,
        show_every = sim.show_every,
        electrolysis_convection = (sim.electrolysis_convection == 1),
        electrolysis_liquid_phase = true,
        electrolysis_phase_change_case = sim.electrolysis_phase_change_case,
        imposed_velocity = sim.imposed_velocity,
        adapt_timestep_mode = sim.adapt_timestep_mode,
        non_dimensionalize = sim.non_dimensionalize,
        mode_2d = sim.mode_2d,
        breakup = sim.breakup,
        # On injecte ici les BC dynamiques présentes
        kwargs_bc...
    )
end

# Lecture des arguments et du fichier YAML
localARGS = ARGS
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

#region study parameters
study = PropertyDict(prop_dict.study)

#mesh connvergence
nb_grid_points = study.meshes
n_cases = length(nb_grid_points)
print("\n number of points ", nb_grid_points, "\n")

#timestep convergence
timesteps = study.timesteps

#endregion study parameters

io = PropertyDict(prop_dict.plot) 
flower = PropertyDict(prop_dict.flower)
mesh = PropertyDict(flower.mesh)
sim = PropertyDict(flower.simulation)
phys = PropertyDict(flower.physics)
macros = PropertyDict(flower.macros) #to parse code from .yml

# boundaries_dict = PropertyDict(macros.boundaries_list)
study_name = ""

# print parameters by evaluating Julia code stored in .yml   
eval(Meta.parseall(macros.print_parameters))


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
    mpi_coords_x = 1; mpi_coords_y = 1
    mpi_max_coords_x = 1; mpi_max_coords_y = 1
    local nx = 32
    local ny = 32
    nstep = 0
    phys_time = 0.0 #Cdouble
    

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
            # "timestep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
            # "nb_Navier_slip_BC"::Cstring, num.nNavier::Ref{Clonglong}, PDI_OUT::Cint,
            C_NULL::Ptr{Cvoid})::Cint

    print("\n PDI INIT status ",PDI_status)

    @debug "after PDI_multi_expose"

    @debug "After full PDI init"

end #if io.pdi>0

# arrays to store errors
error_list_l1 = zeros(n_cases)
error_list_l2 = zeros(n_cases)
error_list_linfty = zeros(n_cases)
error_list_l1_mixed = zeros(n_cases)
error_list_l2_mixed = zeros(n_cases)
error_list_linfty_mixed = zeros(n_cases)
error_list_l1_full = zeros(n_cases)
error_list_l2_full = zeros(n_cases)
error_list_linfty_full = zeros(n_cases)

cell_volume_list = zeros(n_cases)

base_directory = pwd()

#region timestep convergence
for timestep in timesteps
    
    print("\n Timestep ",timestep)
    timestep_to_string = @sprintf "timestep_%.4e" timestep
    # print("\n timestep_to_string ",timestep_to_string)
    timestep_to_string = replace(timestep_to_string, "." => "_")
    # print("\n timestep_to_string ",timestep_to_string)

    mkpath(timestep_to_string)
    cd(timestep_to_string)

    # print("\n pwd() ",pwd())

    #region mesh convergence
    for (i,study_nb_grid_points) in enumerate(nb_grid_points)
        

        mesh_to_string = @sprintf "mesh_%.5i" study_nb_grid_points
       
        mesh_ratio = Int((mesh.ymax - mesh.ymin)/ (mesh.xmax - mesh.xmin) )
        
        print("\n mesh ratio ", mesh_ratio ) 

        study_nb_grid_points_x = study_nb_grid_points
        study_nb_grid_points_y = study_nb_grid_points * mesh_ratio

        print("\n nx ", study_nb_grid_points_x, " ny ", study_nb_grid_points_y)

        mkpath(mesh_to_string)
        cd(mesh_to_string)

        
        # init regular grid
        scalar_mesh_x = collect(LinRange(mesh.xmin, mesh.xmax, study_nb_grid_points_x + 1))    
        scalar_mesh_y = collect(LinRange(mesh.ymin, mesh.ymax, study_nb_grid_points_y + 1))

        @debug "Before Numerical"

        #region test fill struct
        
        #aliases from yml to Flower.jl
        aliases = Dict(
            :x => :scalar_mesh_x,
            :y => :scalar_mesh_y,
            :xcoord => :intfc_x,
            :ycoord => :intfc_y,
            :R => :radius,
            :max_iterations => :max_iter,
            :save_every => :max_iter,
            :ϵ => :epsilon,
            :ϵwall => :epsilon_wall,
            :nLS => :nb_levelsets,
            :θd => :temperature0,
            :β => :beta,
            :σ => :sigma,
            :δreinit => :delta_reinit,
            :n_ext_cl => :n_ext,
            :timestep_0 => :timestep,
            :timestep_n => :timestep,
            :io_pdi => :pdi,
        )

        # fields defined locally, not in yml
        extra = Dict(
            :scalar_mesh_x => scalar_mesh_x,
            :scalar_mesh_y => scalar_mesh_y,
        )


        phase_change_symb = isnothing(sim.phase_change) ? "default" : Symbol(sim.phase_change)

        if isnothing(sim.mass_transfer_redistribute)
            if sim.phase_change_method == 5
                mass_transfer_redistribute_symb = :no
            elseif num.phase_change_method in [6,7]
                mass_transfer_redistribute_symb = :yes
            end
        else
            mass_transfer_redistribute_symb = Symbol(sim.mass_transfer_redistribute)
        end
        # mass_transfer_redistribute_symb = isnothing(sim.mass_transfer_redistribute) ? "no" : Symbol(sim.mass_transfer_redistribute)

       

        default = Numerical{Float64,Int}(
            x = scalar_mesh_x,
            y = scalar_mesh_y,
            timestep_n = timestep,
            timestep_0 = timestep,
            electrolysis_reaction_symb = Symbol(phys.electrolysis_reaction),
            bulk_velocity_symb = Symbol(phys.bulk_velocity),
            phase_change_symb = phase_change_symb,
            mass_transfer_redistribute_symb = mass_transfer_redistribute_symb,
            )  
        # construct default parametric instance with x otherwise L0 and ... not defined in the same way


        # global num_new = safefill_with_aliases(Numerical{Float64, Int}, sim, phys, io,aliases)

        # global num_new = safefill_with_aliases_and_extra(Numerical{Float64,Int},
        #                          sim, phys, io,
        #                          aliases,
        #                          extra)

        print("\n ns_advection ", (sim.ns_advection == 1))

        global num_new = safefill_with_aliases_and_extra_already_init(Numerical{Float64,Int}, default,
                                 sim, phys, io, aliases, extra)

        #endregion test fill struct

    

        global num = num_new

        if num.verbosity >0
            print("\n Parameters num ",num)
        end

        Broadcast.broadcastable(num::Numerical) = Ref(num) #do not broadcast num 
        @debug "After Numerical"

        # num.electrolysis_reaction = Symbol(num.electrolysis_reaction)

        global gp, gu, gv = init_meshes(num)
        #gp, gu, gv = init_meshes(num) does not work for eval(Meta.parseall(macros.boundaries))
        global op, phS, phL = init_fields(num, gp, gu, gv)

        gp.LS[1].u .= 1.0 # deactivate interface

		
		# ------------------------------------------------------------
        # EVALUATION DYNAMIQUE (Macros & BC)
        # ------------------------------------------------------------
        # 1. On exécute TOUTES les macros pour créer les variables dans Main
        eval(Meta.parseall(macros.init_fields))

        # Define boundary conditions
        eval(Meta.parseall(macros.boundaries))
        eval(Meta.parseall(macros.interface)) # <--- DÉPLACÉ ICI (IMPORTANT)

        # 2. Capture des BC dans un dictionnaire pour éviter le "World Age Problem"
        bc_dict = Dict{Symbol, Any}()
        bc_symbols = [:BC_uL, :BC_uS, :BC_vL, :BC_vS, :BC_pL, :BC_pS, :BC_u, :BC_int, :BC_trans_scal, :BC_phi_ele]
        
        for sym in bc_symbols
            if isdefined(Main, sym)
                bc_dict[sym] = getfield(Main, sym)
            else
                # Si la variable n'est pas définie, on ne l'ajoute pas au dict.
                # Cela laissera run_forward! utiliser sa valeur par défaut (si elle existe)
                # ou plantera avec une erreur plus claire si elle est obligatoire.
            end
        end
        # ------------------------------------------------------------

        if num.io_pdi > 0
            # Send meta-data to PDI
            nx = gp.nx
            ny = gp.ny
            
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
                        "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
                        "nb_Navier_slip_BC"::Cstring, num.nNavier::Ref{Clonglong}, PDI_OUT::Cint,
                        "timestep"::Cstring, timestep::Ref{Cdouble}, PDI_OUT::Cint,
                        C_NULL::Ptr{Cvoid})::Cint

            catch error
                print("\n Bug init_PDI \n")

                print("\n PDI_status ",PDI_status, "\n")
                print("\n error",error)

            end

            @debug "after PDI_multi_expose"

            @debug "After full PDI init"

        end #if num.io_pdi>0


        # Define interfaces (for bubbles, drops...)
        eval(Meta.parseall(macros.interface))


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
                
                # printstyled(color=:magenta, @sprintf "\n PDI write_data_start_linftyp %.5i \n" num.current_iter)

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
                "nstep"::Cstring, num.current_iter::Ref{Clonglong}, PDI_OUT::Cint,
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


        velocity_y_bubble = 0.0 .* gp.LS[end].geoS.dcap[:,:,5]
        volume_fraction = gp.LS[end].geoS.dcap[:,:,5]
        center_of_mass_x = 0.0
        center_of_mass_y = 0.0
        circularity = 1.0
        rise_velocity_y = 0.0
        iLSpdi = 1
        velocity_y = velocity_y_bubble
        area_liq = 0.0
        area_gaz = 0.0
        # PDI_status = @ccall "libpdi".PDI_multi_expose("write_postprocessing_rising_bubble"::Cstring,
        # "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        # "center_of_mass_x"::Cstring, center_of_mass_x::Ref{Cdouble}, PDI_OUT::Cint,      
        # "center_of_mass_y"::Cstring, center_of_mass_y::Ref{Cdouble}, PDI_OUT::Cint,      
        # "rise_velocity_y"::Cstring, rise_velocity_y::Ref{Cdouble}, PDI_OUT::Cint,      
        # "circularity"::Cstring, circularity::Ref{Cdouble}, PDI_OUT::Cint,      
        # "velocity_y_bubble"::Cstring, velocity_y_bubble::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "area_gaz"::Cstring, area_gaz::Ref{Cdouble}, PDI_OUT::Cint,      
        # "area_liq"::Cstring, area_liq::Ref{Cdouble}, PDI_OUT::Cint,      
        # C_NULL::Ptr{Cvoid})::Cint

        # PDI_status = @ccall "libpdi".PDI_multi_expose("post_processing_rising_bubble_new"::Cstring,
        # "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        # "velocity_y"::Cstring, velocity_y_bubble::Ptr{Cdouble}, PDI_OUT::Cint,      
        # "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
        # # "volume_liq_cell"::Cstring, grid_p.LS[end].geoL.dcap[:,:,5]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # # "volume_cell"::Cstring, grid_p.LS[end].geoS.dcap[:,:,5]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # "mesh_p_x"::Cstring, gp.x::Ptr{Cdouble}, PDI_OUT::Cint,
        # "mesh_p_y"::Cstring, gp.y::Ptr{Cdouble}, PDI_OUT::Cint,
        # "dcap_1"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # "dcap_2"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # "dcap_3"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # "dcap_4"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        # C_NULL::Ptr{Cvoid})::Cint


        PDI_status = @ccall "libpdi".PDI_multi_expose("post_processing_rising_bubble"::Cstring,
        "nstep"::Cstring, num.current_iter ::Ref{Clonglong}, PDI_OUT::Cint,
        "velocity_y"::Cstring, velocity_y::Ptr{Cdouble}, PDI_OUT::Cint,      
        "volume_fraction"::Cstring, volume_fraction::Ptr{Cdouble}, PDI_OUT::Cint,
        "volume_liq_cell"::Cstring, gp.LS[iLSpdi].geoL.dcap[:,:,5]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        "volume_cell"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,5]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        "mesh_p_x"::Cstring, gp.x::Ptr{Cdouble}, PDI_OUT::Cint,
        "mesh_p_y"::Cstring, gp.y::Ptr{Cdouble}, PDI_OUT::Cint,
        "dcap_1"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,1]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        "dcap_2"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,2]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        "dcap_3"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,3]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        "dcap_4"::Cstring, gp.LS[iLSpdi].geoS.dcap[:,:,4]::Ptr{Cdouble}, PDI_OUT::Cint, #geoS for bubble phase
        C_NULL::Ptr{Cvoid})::Cint

        # -----------------------------------------------------------------------------------------
        # APPEL DU SOLVEUR via FUNCTION BARRIER + INVOKELATEST
        # C'est ici que la correction opère : on passe le dictionnaire 'bc_dict' via invokelatest
        # -----------------------------------------------------------------------------------------
        @debug "Before run_forward wrapper"
        Base.invokelatest(
            execute_simulation_step, 
            num, gp, gu, gv, op, phS, phL, sim, phys, time_scheme, bc_dict
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


        if study.compute_errors != "None" #"Poiseuille"

            LIQUID = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .> (1-1e-16)]
            MIXED = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .<= (1-1e-16) .&& gp.LS[1].geoL.cap[:,:,5] .> 1e-16]

            if study.compute_errors == "Poiseuille"
                norm_all = relative_errors(phL.v, vPoiseuille, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
                norm_mixed = relative_errors(phL.v, vPoiseuille, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
                norm_full = relative_errors(phL.v, vPoiseuille, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)
            elseif study.compute_errors == "diffusion_full_cell"
                
                # Parameters

                L = mesh.xmax-mesh.xmin

                D = num.diffusion_coeff[2] #3.2e-9  #  diffusivity

                analytical_nx = 1000 #accumulating error if n too small so if N big, accumulating too much if nx=100 and N=1000

                diffusion_time_scale = L^2 / D

                println("diffusion time scale", diffusion_time_scale)


                # ele_current_i = minimum(i_current_mag from gradient)
                # print("i mag min ", ele_current_i)

                ele_current_i = -1.59e4

                F1 = ele_current_i / (2 * num.Faraday * D) #BC

                t_max = diffusion_time_scale 

                println("diffusion_time_scale", diffusion_time_scale)

                analytical_dx = L / analytical_nx  # Spatial step size

                analytical_x = collect(0:analytical_dx:L)  # Spatial domain
               
             
                c0 = num.concentration0[2] 

                println("max c without convection", c0 - F1 * L / 2, c0 + F1 * L / 2)

                max_nb_Fourier_series = 1000


                # Special function u1
                function full_cell_u1(x, t)
                    return F1 * x
                    # return (F2-F1)/(2*L) + F1 * x + D*(F2-F1)/L*t #if different left right
                end

                # Initial condition f(x)
                function full_cell_f(x)
                    return c0 - full_cell_u1(x, 0)  # Example initial condition c0-u1?
                end

                # Coefficients A_n
                function full_cell_A_n(n, L, F1, dx,c0)
                    # integral = sum(full_cell_f(xi) * cos(n * π * xi / L) for xi in x) * dx
                    # return (2 / L) * integral #* cos(n * π * x / L)
                    if n == 0
                        return -F1*L +2*c0 #-c0*L*2/L #analytical if just F1 * x 
                    else
                        return (-2*F1*L)/(n*π)^2*((-1)^n-1) #analytical if just F1 * x 
                    end
                end



                # Solution u2
                function full_cell_u2(x, t, L, D, max_nb_Fourier_series,c0,analytical_dx)
                    result = sum(full_cell_A_n(n, L, F1, analytical_dx,c0) * cos(n * π * x / L) * exp(-D * (n * π / L)^2 * t) for n in 1:max_nb_Fourier_series)
                    result += full_cell_A_n(0, L, F1, analytical_dx,c0) * cos(0 * π * x / L) * exp(-D * (0 * π / L)^2 * t) / 2
                    #case 0 special int cos(0) cos(0)dx = L
                    return result
                end

                # Total solution u
                function full_cell_u(x, t, L, D, max_nb_Fourier_series,c0,analytical_dx)
                    return full_cell_u1(x, t) + full_cell_u2(x, t, L, D, max_nb_Fourier_series,c0,analytical_dx)
                end
                

                # print("\n param ",gp.x[2,:]," | ",num.time," | ",L," | ",num.diffusion_coeff[2]," | ",max_nb_Fourier_series," | ",c0," | ",analytical_dx)
                #TODO NX vs n ... need to interpolate
                # concentration_profile = full_cell_u.(analytical_x,num.time,L,num.diffusion_coeff[2],max_nb_Fourier_series,c0,analytical_dx)
                concentration_profile = full_cell_u.(gp.x[2,:],num.time,L,num.diffusion_coeff[2],max_nb_Fourier_series,c0,analytical_dx)
                
                print("\n concentration profile ",phL.trans_scal[2,:,2])

                print("\n concentration profile ",concentration_profile)

                # norm_all = relative_errors(phL.trans_scal[:,:,2], concentration_profile, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
                # norm_mixed = relative_errors(phL.trans_scal[:,:,2], concentration_profile, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
                # norm_full = relative_errors(phL.trans_scal[:,:,2], concentration_profile, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)
               
                
                l1,l2,linfty = relative_errors(phL.trans_scal[:,:,2], concentration_profile, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
                l1_mixed,l2_mixed,linfty_mixed = relative_errors(phL.trans_scal[:,:,2], concentration_profile, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
                l1_full,l2_full,linfty_full = relative_errors(phL.trans_scal[:,:,2], concentration_profile, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)
            end


            error_list_l1[i] = l1
            error_list_l2[i] = l2
            error_list_linfty[i] = linfty

            error_list_l1_mixed[i] = l1_mixed
            error_list_l2_mixed[i] = l2_mixed
            error_list_linfty_mixed[i] = linfty_mixed

            error_list_l1_full[i] = l1_full
            error_list_l2_full[i] = l2_full
            error_list_linfty_full[i] = linfty_full

            cell_volume_list[i] = minimum(gp.LS[1].geoL.dcap[:,:,5])
           
            min_cell_volume = minimum(gp.LS[1].geoL.dcap[:,:,5])

            if study.compute_errors != "None" #"Poiseuille"

                local PDI_status = @ccall "libpdi".PDI_multi_expose("convergence_study_iter"::Cstring, 
                "study_nb_grid_points"::Cstring, study_nb_grid_points::Ref{Clong}, PDI_OUT::Cint,
                "study_timestep"::Cstring, timestep::Ref{Cdouble}, PDI_OUT::Cint,
                # "cell_volume"::Cstring, cell_volume_list::Ptr{Cdouble}, PDI_OUT::Cint,
                "study_l1_rel_error"::Cstring, l1::Ref{Cdouble}, PDI_OUT::Cint,
                "study_l2_rel_error"::Cstring, l2::Ref{Cdouble}, PDI_OUT::Cint,
                "study_linfty_rel_error"::Cstring, linfty::Ref{Cdouble}, PDI_OUT::Cint,
                "study_l1_rel_error_full_cells"::Cstring, l1_full::Ref{Cdouble}, PDI_OUT::Cint,
                "study_l2_rel_error_full_cells"::Cstring, l2_full::Ref{Cdouble}, PDI_OUT::Cint,
                "study_linfty_rel_error_full_cells"::Cstring, linfty_full::Ref{Cdouble}, PDI_OUT::Cint,
                "study_l1_rel_error_partial_cells"::Cstring, l1_mixed::Ref{Cdouble}, PDI_OUT::Cint,
                "study_l2_rel_error_partial_cells"::Cstring, l2_mixed::Ref{Cdouble}, PDI_OUT::Cint,
                "study_linfty_rel_error_partial_cells"::Cstring, linfty_mixed::Ref{Cdouble}, PDI_OUT::Cint,
                "domain_length"::Cstring, L0::Ref{Cdouble}, PDI_OUT::Cint,
                "min_cell_volume"::Cstring, min_cell_volume::Ref{Cdouble}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cint

            end

        end

        current_directory = pwd()
        if current_directory == base_directory*"/"*timestep_to_string*"/"*mesh_to_string
            cd("..")
        else 
            print("\n Error in mesh subdirectory",current_directory," not ",base_directory*"/"*timestep_to_string*"/"*mesh_to_string)
            exit()
        end

    end 
    #endregion mesh convergence

    current_directory = pwd()
    if current_directory == base_directory*"/"*timestep_to_string
        cd("..")
    else
        print("\n Error in timestep subdirectory",current_directory," not ",base_directory*"/"*timestep_to_string)
        exit()

    end

end
#endregion timestep convergence


min_cell_volume = minimum(gp.LS[1].geoL.cap[:,:,5])

print("\n min_cell_volume ",min_cell_volume," type ",typeof(min_cell_volume))

if study.compute_errors != "None" #"Poiseuille"

    local PDI_status = @ccall "libpdi".PDI_multi_expose("convergence_study"::Cstring, 
    "n_tests"::Cstring, n_cases::Ref{Clonglong}, PDI_OUT::Cint,
    "nx_list"::Cstring, nb_grid_points::Ptr{Clonglong}, PDI_OUT::Cint,
    "cell_volume_list"::Cstring, cell_volume_list::Ptr{Cdouble}, PDI_OUT::Cint,
    "l1_rel_error"::Cstring, error_list_l1::Ptr{Cdouble}, PDI_OUT::Cint,
    "l2_rel_error"::Cstring, error_list_l2::Ptr{Cdouble}, PDI_OUT::Cint,
    "linfty_rel_error"::Cstring, error_list_linfty::Ptr{Cdouble}, PDI_OUT::Cint,
    "l1_rel_error_full_cells"::Cstring, error_list_l1_full::Ptr{Cdouble}, PDI_OUT::Cint,
    "l2_rel_error_full_cells"::Cstring, error_list_l2_full::Ptr{Cdouble}, PDI_OUT::Cint,
    "linfty_rel_error_full_cells"::Cstring, error_list_linfty_full::Ptr{Cdouble}, PDI_OUT::Cint,
    "l1_rel_error_partial_cells"::Cstring, error_list_l1_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
    "l2_rel_error_partial_cells"::Cstring, error_list_l2_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
    "linfty_rel_error_partial_cells"::Cstring, error_list_linfty_mixed::Ptr{Cdouble}, PDI_OUT::Cint,
    "domain_length"::Cstring, L0::Ref{Cdouble}, PDI_OUT::Cint,
    "min_cell_volume"::Cstring, min_cell_volume::Ref{Cdouble}, PDI_OUT::Cint,
    C_NULL::Ptr{Cvoid})::Cint

end
 
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
if haskey(macros, "test_end")
    eval(Meta.parseall(macros.test_end))
end