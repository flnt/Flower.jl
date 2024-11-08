# using Revise


# using PrecompileTools_startup

# module testyamlfile3 #enables to perform a test with ARGS to give an input file
# ARGS = String["../examples/levelset_Butler.yml"]
# include("../examples/main_current_folder.jl")
# end

using Flower


# using PrettyTables

# using SnoopCompileCore


makedir = false 

# localARGS = isdefined(:newARGS) ? newARGS : ARGS 
localARGS = ARGS
@show localARGS

# print("\n Arguments ", localARGS)
# print("\n length(localARGS) ",length(localARGS))

if length(localARGS)>0
    yamlfile = localARGS[1]
    printstyled(color=:magenta, @sprintf "\n YAML file ")
    print(yamlfile)
    print("\n")

    prefix = "."

    yamlpath = yamlfile

    yamlnamelist = split(yamlfile, "/")

    yamlname = yamlnamelist[end]
    
    # cp(yamlpath,"./"*yamlname,force=true)
end

data = YAML.load_file(yamlpath)

# Dictionaries 
prop_dict = PropertyDict(data)
io = PropertyDict(prop_dict.plot) 
flower = PropertyDict(prop_dict.flower)
mesh = PropertyDict(flower.mesh)
sim = PropertyDict(flower.simulation)
phys = PropertyDict(flower.physics)


# Simulation parameters
save_every = sim.max_iter

radial_vel_factor = 1e-7

pressure_channel = false

current_radius = 0.0

ns_advection=(sim.ns_advection ==1)



# print("\n ns_advection ",ns_advection," ",sim.ns_advection)

# #TODO restart with PDI
# if sim.restart == 1:
#     cf hello_access.jl
# end

if sim.name == "falling_drop"
    L0x = 4.0
    L0y = 6.0
        
    # x = collect(LinRange(-L0x / 2, L0x / 2, mesh.nx+1))
    x = LinRange(mesh.xmin, mesh.xmax, mesh.nx+1)

    dx = diff(x)[1]
    y = collect(-L0y/2:dx:L0y/2+dx)

    A = zeros(8)
    A[1] = 0.1
      
    Re = 1.0
    # dt in branch free-surface 
    deltax = min(diff(x)..., diff(y)...)
    sim.dt0 = min(sim.CFL*deltax^2*Re, sim.CFL*deltax)

    print("\n dt0 ", sim.dt0)

    # advection deactivated

elseif occursin("sessile",sim.name) #sim.name == "sessile_2LS"
 
    # Physical parameters 
    x = LinRange(mesh.xmin, mesh.xmax, mesh.nx+1)    
    y = collect(LinRange(mesh.ymin, mesh.ymax, mesh.ny + 1))

    h0 = 0.5
    Re = 1.0
 
    # not imposing the angle exactly at the boundary but displaced a cell because ghost cells are not used. 
    # So the exact contact angle cannot be imposed
    # TODO should at some point modify the levelset so it works as the other fields that we have in the code, 
    # with the boundary values in addition to the bulk field. That way we could impose it exactly at the boundary
    # needs a lot of work, not priority
    _θe = acos((0.5 * diff(y)[1] + cos(phys.theta_e * π / 180) * h0) / h0) * 180 / π
    
    # _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
    # println("θe = $(_θe)")

else
    

    # Physical parameters 
    x = LinRange(mesh.xmin, mesh.xmax, mesh.nx+1)
    y = LinRange(mesh.ymin, mesh.ymax, mesh.ny+1)


    # Physics
    mu = phys.mu_cin1 *phys.rho1 #in Pa s = M L^{-1} T^{-1}}

    phys.mu1 = mu
    phys.mu2 = mu

    mu1=mu
    mu2=mu 

    h0 = phys.radius
    ref_thickness_2d = 4.0 / 3.0 *phys.radius 

    c0_H2,c0_KOH,c0_H2O = phys.concentration0

    DH2,DKOH,DH2O= phys.diffusion_coeff

    #not using num.bulk_conductivity like in run.jl as long as the initial concentration is independent from small cell isues
    elec_cond=2*phys.Faraday^2*c0_KOH*DKOH/(phys.Ru*phys.temperature0)

    Re=phys.rho1*phys.v_inlet*phys.ref_length/mu #Reynolds number
    printstyled(color=:green, @sprintf "\n Re : %.2e %.2e %.2e %.2e\n" Re phys.rho1/mu1 phys.rho1 mu1)

    Re=phys.rho1/mu1 #not Reynolds number, but rho1/mu1

    printstyled(color=:green, @sprintf "\n 'Re' i.e. rho/mu : %.2e %.2e %.2e %.2e\n" Re phys.rho1/mu1 phys.rho1 mu1)

    if length(phys.concentration0)!=phys.nb_transported_scalars
        print(@sprintf "nb_transported_scalars: %5i\n" phys.nb_transported_scalars)
        @error ("nb_transported_scalars")
    end

    if length(phys.diffusion_coeff)!=phys.nb_transported_scalars
        print(@sprintf "nb_transported_scalars: %5i\n" phys.nb_transported_scalars)
        @error ("nb_transported_scalars")
    end

    print(@sprintf "nb_transported_scalars: %5i\n" phys.nb_transported_scalars)

    diffusion_t = (phys.radius^2)./phys.diffusion_coeff
    print("\n diffusion time ", diffusion_t,"\n")
    
    # print_table = false
    ## print_table = true

    # if print_table
    #     hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold cyan")

    #     diffusion_t = (phys.radius^2)./phys.diffusion_coeff

    #     pretty_table(vcat(
    #         hcat("Diffusion time",diffusion_t'),
    #         hcat("Diffusion coef",phys.diffusion_coeff'),
    #         hcat("Concentration",phys.concentration0')); 
    #     formatters    = ft_printf("%0.2e", 2:4), #not ecessary , 2:4
    #     header = ["","H2", "KOH", "H2O"], 
    #     highlighters=hl)

    #     @debug "After Table"
    # end


    # printstyled(color=:green, @sprintf "\n Species diffusion timescales: %.2e %.2e %.2e \n" (phys.radius^2)/DH2 (phys.radius^2)/DKOH (phys.radius^2)/DH2O )

    # printstyled(color=:green, @sprintf "\n nmol : \n" phys.concentration0[1]*4.0/3.0*pi*phys.radius^3  )

    current_radius = phys.radius

    p_liq= phys.pres0 #+ mean(veci(phL.pD,grid,2)) #TODO here one bubble
    # p_g=p_liq + 2 * phys.sigma / current_radius
    p_g=p_liq + phys.sigma / current_radius

    c0test = p_g / (phys.temperature0 * phys.Ru) 

    # printstyled(color=:green, @sprintf "\n c0test: %.2e \n" c0test)

    # printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" phys.concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(phys.temperature0*phys.Ru))

    #dummy gp
    gp = Mesh(GridCC, x, y, 1)

    # BC 
    i_butler = gp.x[:,1] .*0.0
    phi_ele =  gp.x[:,1] .*0.0

    # print("\n i_butler ", i_butler)
    # print("\n phi_ele ", phi_ele)


    i_butler=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele,phys.phi_ele1,phys.Ru,phys.temperature0)

    # print(@sprintf "Butler-Volmer %.2e %.2e %.2e %.2e\n" i_butler[1] -i_butler[1]/(2*phys.Faraday*DH2) c0_H2-i_butler[1]/(2*phys.Faraday*DH2)*gp.dx[1,1] c0_H2+i_butler[1]/(2*phys.Faraday*DH2)*gp.dx[1,1])


    # Pressure
    p_top = 0
    p_bottom = p_top + 8*mu1/phys.ref_length*phys.v_inlet

    # H2 boundary condition
    BC_trans_scal_H2 = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[1]),
    top    = Neumann(),
    left   = Neumann(val=-i_butler/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
    right  = Dirichlet(val = phys.concentration0[1]),
    int    = Dirichlet(val = phys.concentration0[1]))

    BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=-i_butler/(2*phys.Faraday*DKOH)),
        right  = Dirichlet(val = phys.concentration0[2]),
        int    = Neumann(val=0.0)) #KOH
        
    BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=i_butler/(phys.Faraday*DH2O)),
        right  = Dirichlet(val = phys.concentration0[3]),
        int    = Neumann(val=0.0)) #Dirichlet(val = phys.concentration0[3])),#Neumann(val=0.0)) 
        #H2O #Dirichlet(val = phys.concentration0[3])

    BC_pL = Boundaries(
            left   = Neumann(val=0.0),
            right  = Neumann(val=0.0),
            bottom = Neumann(val=0.0),
            top    = Dirichlet(),
        )

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Neumann(),
        top    = Neumann(),
    )


    if sim.name == "small_cell"
        #Test case 1: small cells (high concentration at vecb_L)
        sim.max_iter = 1
        save_every = 1

        folder="electrolysis_circle_wall_CFL_small_cell"

    elseif sim.name == "100it"
        #Test case 2: scalar without velocity
        sim.max_iter = 100
        save_every = 10
        save_p = false

        folder="electrolysis_circle_wall_CFL"*sim.name

    elseif sim.name == "radial"
        sim.imposed_velocity = "radial"
        sim.max_iter = 1
        save_every = 1
        mesh.nx = 64
        radial_vel_factor = 1e-7
        folder="electrolysis_circle_wall_CFL"*sim.name

    elseif sim.name == "channel_no_bubble"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25

        # sim.max_iter = 2
        # save_every = 1

    elseif sim.name == "channel_Dirichlet"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25

        # sim.max_iter = 2
        # save_every = 1
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[1]),
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1])) #H2

        phys.electrolysis_reaction = "none"

    elseif sim.name == "channel_Dirichlet_pressure"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25

        sim.max_iter = 2
        save_every = 1
        sim.max_iter = 1

        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[1]),
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1])) #H2


        phys.electrolysis_reaction = "none"
        # sim.imposed_velocity = "constant"

        

        test_v=-(p_top - p_bottom)/phys.ref_length * (phys.ref_length^2)/4/2/mu

        print("\n phys.ref_length ",phys.ref_length)
        
        printstyled(color=:green, @sprintf "\n mu1 : %.2e v %.2e vtest %.2e sim.CFL %.2e \n" mu1 phys.v_inlet test_v phys.v_inlet*sim.dt0/(phys.ref_length/mesh.nx))

        printstyled(color=:green, @sprintf "\n p_bottom : %.2e p_top %.2e grad %.2e grad %.2e \n" p_bottom p_top -(p_top-p_bottom)/phys.ref_length 8*mu1/phys.ref_length^2*phys.v_inlet)


        # print("\n mu1 ",mu1," phys.ref_length ", phys.ref_length)

        pressure_channel = true

        BC_vL= Boundaries(
            left   = Dirichlet(),
            right  = Dirichlet(),
            bottom = Neumann(),
            top    = Neumann(),
        )

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Dirichlet(val = p_bottom),
            top    = Dirichlet(val = p_top),
        )

    elseif sim.name == "channel_Dirichlet_constant_vel"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25

        save_every = sim.max_iter

        sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet/2 
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[1]),
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[2]),
        right  = Dirichlet(val = phys.concentration0[2]),
        int    = Dirichlet(val = phys.concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[3]),
        right  = Dirichlet(val = phys.concentration0[3]),
        int    = Dirichlet(val = phys.concentration0[3])) #H2O


        phys.electrolysis_reaction = "none"
        sim.imposed_velocity = "constant"

        BC_vL= Boundaries(
            left   = Dirichlet(val = phys.v_inlet),
            right  = Dirichlet(val = phys.v_inlet),
            bottom = Dirichlet(val = phys.v_inlet),
            top    = Neumann(val=0.0),
        )

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Neumann(),
            top    = Neumann(),
        )

    elseif sim.name == "channel_Dirichlet_zero_vel"
        sim.activate_interface = 0
    
        sim.max_iter = 1

        save_every = sim.max_iter

        sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet/2 
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[1]),
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[2]),
        right  = Dirichlet(val = phys.concentration0[2]),
        int    = Dirichlet(val = phys.concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[3]),
        right  = Dirichlet(val = phys.concentration0[3]),
        int    = Dirichlet(val = phys.concentration0[3])) #H2O


        phys.electrolysis_reaction = "none"
        sim.imposed_velocity = "none"

        BC_vL= Boundaries(
            left   = Dirichlet(val = phys.v_inlet),
            right  = Dirichlet(val = phys.v_inlet),
            bottom = Dirichlet(val = phys.v_inlet),
            top    = Neumann(val=0.0),
        )

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Neumann(),
            top    = Neumann(),
        )

    elseif sim.name == "channel_Dirichlet_imposed_Poiseuille"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25

        save_every = sim.max_iter

        # sim.max_iter = 1
        # save_every = 1
        # # sim.CFL 1
        # sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet 
        # sim.CFL 0.5
        sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet/2 

        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[1]),
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[2]),
        right  = Dirichlet(val = phys.concentration0[2]),
        int    = Dirichlet(val = phys.concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = phys.concentration0[3]),
        right  = Dirichlet(val = phys.concentration0[3]),
        int    = Dirichlet(val = phys.concentration0[3])) #H2O


        phys.electrolysis_reaction = "none"
        sim.imposed_velocity = "Poiseuille_bottom_top"

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Neumann(),
            top    = Neumann(),
        )

    elseif sim.name == "channel_no_bubble_Cdivu"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25
        sim.convection_Cdivu = true

    elseif sim.name == "channel_no_bubble_no_vel"
        sim.activate_interface = 0
        sim.max_iter = 100
        save_every = 25
    elseif sim.name == "channel"
        sim.max_iter = 100
        save_every = 25
    elseif sim.name == "imposed_radius"
        electrolysis_phase_change = true
        sim.max_iter = 100
        save_every = 25
        sim.electrolysis_phase_change_case = "imposed_radius"
    end


    folder = sim.name

    # Save path
    if makedir
        prefix *= "/"*folder*"/"
        isdir(prefix) || mkdir(prefix)
        # yamlfile
        yamlfile2="flower.yml"
        cp(yamlpath,prefix*yamlfile2,force=true) #copy yaml file to simulation directory
    end

end #if sim.name ==


@debug "Before Numerical"





num = Numerical(
    CFL = sim.CFL,
    Re = Re,
    TEND=phys.end_time,
    x = x,
    y = y,
    xcoord = phys.intfc_x,
    ycoord = phys.intfc_y,
    case = sim.case,
    R = phys.radius,
    max_iterations = sim.max_iter,
    save_every = save_every,
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
    plot_prefix = prefix,
    dt0 = sim.dt0,
    concentration_check_factor = sim.concentration_check_factor,
    radial_vel_factor = radial_vel_factor,
    debug = sim.debug,
    v_inlet = phys.v_inlet,
    prediction = sim.prediction,
    null_space = sim.null_space,
    io_pdi = io.pdi,
    bulk_conductivity = sim.bulk_conductivity,
    electric_potential = sim.electric_potential,
    contact_angle = sim.contact_angle,
    convection_Cdivu = sim.convection_Cdivu,
    convection_mode = sim.convection_mode,
    advection_LS_mode = sim.advection_LS_mode,
    scalar_bc = sim.scalar_bc,
    scalar_scheme = sim.scalar_scheme,
    solver = sim.solver,
    mass_flux = sim.mass_flux,
    average_liquid_solid = sim.average_liquid_solid,
    )

Broadcast.broadcastable(num::Numerical) = Ref(num) #do not broadcast num 

@debug "After Numerical"


#Initialization
gp, gu, gv = init_meshes(num)
op, phS, phL = init_fields(num, gp, gu, gv)

print("\n TODO BC_uL")
BC_uL = Boundaries(
        left   = Dirichlet(),#Navier_cl(λ = 1e-2), #Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Neumann(val=0.0),
    )

BC_phi_ele = BoundariesInt(
    left   = Neumann(val=i_butler./elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
    right  = Dirichlet(),
    bottom = Neumann(val=0.0),
    top    = Neumann(val=0.0),
    int    = Neumann(val=0.0),
    LS = [Neumann(val=0.0)]
)

if sim.name == "falling_drop"

    # r = 0.5
    # gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .- y0).^2) - r * ones(gp)
    # gp.LS[1].u .*= -1.0

    wall_shape = (
        A[1] .* cos.(2π .* gp.x) .+ A[2] .* sin.(2π .* gp.x) .+
        A[3] .* cos.(2π .* gp.x) .^ 2 .+ A[4] .* sin.(2π .* gp.x) .^ 2 .+
        A[5] .* cos.(2π .* gp.x) .^ 3 .+ A[6] .* sin.(2π .* gp.x) .^ 3 .+
        A[7] .* cos.(2π .* gp.x) .^ 4 .+ A[8] .* sin.(2π .* gp.x) .^ 4
    )

    # u1 = sqrt.((gp.x .+ 0.49).^2 + (gp.y .+ 0.2).^2) - r * ones(gp)
    # u2 = sqrt.((gp.x .- 0.49).^2 + (gp.y .+ 0.2).^2) - r * ones(gp)
    # gp.LS[2].u .= combine_levelsets(gp, u1, u2)

    gp.LS[2].u .= gp.y .- (L0y/2 .- 0.75) .+ wall_shape
    gp.LS[2].u .*= -1.0

    phL.u .= 0.0
    phL.v .= 0.0
      
    BC_uL = Boundaries(bottom = Navier_cl(λ = 1e-2),)
    BC_vL = Boundaries(bottom = Dirichlet(),)
    BC_pL = Boundaries()
    BC_u = Boundaries(
        bottom = Neumann_cl(θe = phys.theta_e2 * π / 180),
        top = Neumann(),
        left = Neumann_inh(),
        right = Neumann_inh()
    )
    BC_int = [FreeSurface(), WallNoSlip(θe = phys.theta_e * π / 180)]

elseif sim.name == "levelset_Butler_two_LS"

    gp.LS[2].u .= gp.x .- phys.ls_wall_xmin
    # gp.LS[2].u .*= -1.0

    # BC_vL = Boundaries(left=Dirichlet(val=gv.y[:,1]))

    BC_vL = Boundaries()

    BC_pL = Boundaries()

    BC_int = [WallNoSlip(),WallNoSlip()]

    phi_ele_scal = 0.0 # gp.x[:,1] .*0.0


    i_butler_scal=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele_scal,phys.phi_ele1,phys.Ru,phys.temperature0)


    BC_trans_scal_H2 = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[1]),
    top    = Neumann(),
    left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
    right  = Dirichlet(val = phys.concentration0[1]),
    int    = Dirichlet(val = phys.concentration0[1]),
    # LS     = [Dirichlet(val = 10),Dirichlet(val = -10)]
    LS     = [Dirichlet(val = phys.concentration0[1]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2))]
    ) #H2

    BC_trans_scal_KOH = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[2]),
    top    = Neumann(),
    left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
    right  = Dirichlet(val = phys.concentration0[2]),
    int    = Neumann(val=0.0),
    LS     = [Neumann(val=0.0),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
    # int    = Dirichlet(val = phys.concentration0[2]),
    # LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
    # LS     = [Neumann(val=0.0),Dirichlet(val = phys.concentration0[2])]
    # LS     = [Dirichlet(val=phys.concentration0[2]),Dirichlet(val = phys.concentration0[2])]
    ) #KOH

    # BC_trans_scal_KOH = BoundariesInt(
    #     bottom = Dirichlet(val = phys.concentration0[2]),
    #     top    = Neumann(),
    #     left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
    #     right  = Dirichlet(val = phys.concentration0[2]),
    #     int    = Dirichlet(val = phys.concentration0[2]),
    #     LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))]
    #     ) #KOH
    


    BC_trans_scal_H2O = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[3]),
    top    = Neumann(),
    left   = Neumann(val=i_butler_scal/(phys.Faraday*DH2O)),
    right  = Dirichlet(val = phys.concentration0[3]),
    int    = Neumann(val=0.0),
    LS     = [Neumann(val=0.0),Neumann(val=i_butler_scal/(phys.Faraday*DH2O))]
    ) #H2O


    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=i_butler_scal./elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
        LS = [Neumann(val=0.0),Neumann(val=i_butler_scal./elec_cond)]) #first is interface


    elseif sim.name == "test_levelset_Butler_two_LS"

        gp.LS[2].u .= gp.x .- phys.ls_wall_xmin
        # gp.LS[2].u .*= -1.0
    
        # BC_vL = Boundaries(left=Dirichlet(val=gv.y[:,1]))
    
        BC_vL = Boundaries()
    
        BC_pL = Boundaries()
    
        BC_int = [WallNoSlip(),WallNoSlip()]
    
        phi_ele_scal = 0.0 # gp.x[:,1] .*0.0
    
    
        i_butler_scal=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele_scal,phys.phi_ele1,phys.Ru,phys.temperature0)
    
    
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1]),
        # LS     = [Dirichlet(val = 10),Dirichlet(val = -10)]
        LS     = [Dirichlet(val = phys.concentration0[1]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2))]
        ) #H2
    
        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
        right  = Dirichlet(val = phys.concentration0[2]),
        int    = Neumann(val=0.0),
        LS     = [Neumann(val=0.0),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
        # int    = Dirichlet(val = phys.concentration0[2]),
        # LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
        # LS     = [Neumann(val=0.0),Dirichlet(val = phys.concentration0[2])]
        # LS     = [Dirichlet(val=phys.concentration0[2]),Dirichlet(val = phys.concentration0[2])]
        ) #KOH
    
        # BC_trans_scal_KOH = BoundariesInt(
        #     bottom = Dirichlet(val = phys.concentration0[2]),
        #     top    = Neumann(),
        #     left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
        #     right  = Dirichlet(val = phys.concentration0[2]),
        #     int    = Dirichlet(val = phys.concentration0[2]),
        #     LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))]
        #     ) #KOH
        
    
    
        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=i_butler_scal/(phys.Faraday*DH2O)),
        right  = Dirichlet(val = phys.concentration0[3]),
        int    = Neumann(val=0.0),
        LS     = [Neumann(val=0.0),Neumann(val=i_butler_scal/(phys.Faraday*DH2O))]
        ) #H2O
    
    
        BC_phi_ele = BoundariesInt(
            left   = Neumann(val=i_butler_scal./elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
            right  = Neumann(val=i_butler_scal./elec_cond),
            bottom = Neumann(val=i_butler_scal./elec_cond),
            top    = Neumann(val=i_butler_scal./elec_cond),
            int    = Neumann(val=i_butler_scal./elec_cond),
            LS = [Neumann(val=i_butler_scal./elec_cond),Neumann(val=i_butler_scal./elec_cond)]) #first is interface

  
    elseif sim.name == "test_levelset_Butler_two_LS_3"

            gp.LS[2].u .= gp.x .- phys.ls_wall_xmin
            # gp.LS[2].u .*= -1.0
        
            # BC_vL = Boundaries(left=Dirichlet(val=gv.y[:,1]))
        
            BC_vL = Boundaries()
        
            BC_pL = Boundaries()
        
            BC_int = [WallNoSlip(),WallNoSlip()]
        
            phi_ele_scal = 0.0 # gp.x[:,1] .*0.0
        
        
            i_butler_scal=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele_scal,phys.phi_ele1,phys.Ru,phys.temperature0)
        
        
            BC_trans_scal_H2 = BoundariesInt(
            bottom = Dirichlet(val = phys.concentration0[1]),
            top    = Neumann(),
            left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
            right  = Dirichlet(val = phys.concentration0[1]),
            int    = Dirichlet(val = phys.concentration0[1]),
            # LS     = [Dirichlet(val = 10),Dirichlet(val = -10)]
            LS     = [Dirichlet(val = phys.concentration0[1]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2))]
            ) #H2
        
            BC_trans_scal_KOH = BoundariesInt(
            bottom = Dirichlet(val = phys.concentration0[2]),
            top    = Neumann(),
            left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
            right  = Dirichlet(val = phys.concentration0[2]),
            int    = Neumann(val=0.0),
            LS     = [Neumann(val=0.0),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
            # int    = Dirichlet(val = phys.concentration0[2]),
            # LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))],
            # LS     = [Neumann(val=0.0),Dirichlet(val = phys.concentration0[2])]
            # LS     = [Dirichlet(val=phys.concentration0[2]),Dirichlet(val = phys.concentration0[2])]
            ) #KOH
        
            # BC_trans_scal_KOH = BoundariesInt(
            #     bottom = Dirichlet(val = phys.concentration0[2]),
            #     top    = Neumann(),
            #     left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
            #     right  = Dirichlet(val = phys.concentration0[2]),
            #     int    = Dirichlet(val = phys.concentration0[2]),
            #     LS     = [Dirichlet(val = phys.concentration0[2]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))]
            #     ) #KOH
            
        
        
            BC_trans_scal_H2O = BoundariesInt(
            bottom = Dirichlet(val = phys.concentration0[3]),
            top    = Neumann(),
            left   = Neumann(val=i_butler_scal/(phys.Faraday*DH2O)),
            right  = Dirichlet(val = phys.concentration0[3]),
            int    = Neumann(val=0.0),
            LS     = [Neumann(val=0.0),Neumann(val=i_butler_scal/(phys.Faraday*DH2O))]
            ) #H2O
        
        
            BC_phi_ele = BoundariesInt(
                left   = Neumann(val=i_butler_scal./elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
                right  = Dirichlet(),
                bottom = Neumann(val=0.0),
                top    = Neumann(val=0.0),
                int    = Neumann(val=0.0),
                LS = [Neumann(val=i_butler_scal./elec_cond),Neumann(val=i_butler_scal./elec_cond)]) #first is interface
        

elseif sim.name == "levelset_Butler_two_LS_test_isca1_1_2"

    gp.LS[2].u .= gp.x .- phys.ls_wall_xmin
    # gp.LS[2].u .*= -1.0

    # BC_vL = Boundaries(left=Dirichlet(val=gv.y[:,1]))

    BC_vL = Boundaries()

    BC_pL = Boundaries()

    BC_int = [WallNoSlip(),WallNoSlip()]

    phi_ele_scal = 0.0 # gp.x[:,1] .*0.0


    i_butler_scal=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele_scal,phys.phi_ele1,phys.Ru,phys.temperature0)


    BC_trans_scal_H2 = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[1]),
    top    = Neumann(),
    left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
    right  = Dirichlet(val = phys.concentration0[1]),
    int    = Dirichlet(val = phys.concentration0[1]),
    # LS     = [Dirichlet(val = 10),Dirichlet(val = -10)]
    LS     = [Dirichlet(val = phys.concentration0[1]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2))]
    ) #H2

    # BC_trans_scal_KOH = BoundariesInt(
    # bottom = Dirichlet(val = phys.concentration0[2]),
    # top    = Neumann(),
    # left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH)),
    # right  = Dirichlet(val = phys.concentration0[2]),
    # int    = Neumann(val=0.0),
    # LS     = [Neumann(val=0.0),Neumann(val=-i_butler_scal/(2*phys.Faraday*DKOH))]
    # ) #KOH


    BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = phys.concentration0[1]),
        top    = Neumann(),
        left   = Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2)), #Dirichlet(val = phys.concentration0[1]), #
        right  = Dirichlet(val = phys.concentration0[1]),
        int    = Dirichlet(val = phys.concentration0[1]),
        # LS     = [Dirichlet(val = 10),Dirichlet(val = -10)]
        LS     = [Dirichlet(val = phys.concentration0[1]),Neumann(val=-i_butler_scal/(2*phys.Faraday*DH2))]
        ) #H2

    BC_trans_scal_H2O = BoundariesInt(
    bottom = Dirichlet(val = phys.concentration0[3]),
    top    = Neumann(),
    left   = Neumann(val=i_butler_scal/(phys.Faraday*DH2O)),
    right  = Dirichlet(val = phys.concentration0[3]),
    int    = Neumann(val=0.0),
    LS     = [Neumann(val=0.0),Neumann(val=i_butler_scal/(phys.Faraday*DH2O))]
    ) #H2O


    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=i_butler_scal./elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
        LS = [Neumann(val=0.0),Neumann(val=i_butler_scal./elec_cond)]) #first is interface


elseif sim.name == "sessile"

    BC_uL = Boundaries(
        bottom = Navier_cl(λ = 1e-2),
        top = Dirichlet(),
    )

    BC_vL = Boundaries(
        bottom = Dirichlet(),
        top = Dirichlet(),
    )

    BC_pL = Boundaries(
        left = Dirichlet(),
        right = Dirichlet(),
    )

    BC_u = Boundaries(
        bottom = Neumann_cl(θe = _θe * π / 180),
        top = Neumann_inh(),
        left = Neumann_inh(),
        right = Neumann_inh()
    )

    BC_int = [FreeSurface()]

elseif sim.name == "sessile_2LS"

    gp.LS[2].u .= gp.y .+ 0.85  .+ 0.1 .* cos.(2.0 .*π .* gp.x)

    phL.u .= 0.0
    phL.v .= 0.0

    BC_uL = Boundaries(
    bottom = Dirichlet(),
    top = Dirichlet(),
    )

    BC_vL = Boundaries(
    bottom = Dirichlet(),
    top = Dirichlet(),
    )

    BC_pL = Boundaries()

    BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_inh(),
    right = Neumann_inh()
    )

    BC_int = [FreeSurface(), WallNoSlip(θe = phys.theta_e * π / 180, θadv = phys.theta_adv * π / 180, θrec = phys.theta_rec * π / 180)]

elseif sim.name == "sessile_2LS_adv"

    gp.LS[2].u .= gp.y .+ 0.85  #.+ 0.1 .* cos.(2.0 .*π .* gp.x)

    dx = diff(x)[1]

    v_adv = dx/sim.dt0*0.5

    phL.u .= v_adv
    phL.v .= 0.0

    BC_uL = Boundaries(
    bottom = Dirichlet(val=v_adv),
    top = Dirichlet(val=v_adv),
    left= Dirichlet(val=v_adv),
    right = Dirichlet(val=v_adv),
    )

    BC_vL = Boundaries(
    bottom = Dirichlet(),
    top = Dirichlet(),
    )

    BC_pL = Boundaries()

    BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_inh(),
    right = Neumann_inh()
    )

    BC_int = [FreeSurface(), WallNoSlip(θe = phys.theta_e * π / 180, θadv = phys.theta_adv * π / 180, θrec = phys.theta_rec * π / 180)]

elseif sim.name == "sessile_2LS_inclined"

    gp.LS[2].u .= gp.y .+ 0.85  - 0.1 * gp.x  #.+ 0.1 .* cos.(2.0 .*π .* gp.x)



    dx = diff(x)[1]


    phL.u .= 0.0
    phL.v .= 0.0

    BC_uL = Boundaries(
    bottom = Dirichlet(),
    top = Dirichlet(),
    )

    BC_vL = Boundaries(
    bottom = Dirichlet(),
    top = Dirichlet(),
    )

    BC_pL = Boundaries()

    BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_inh(),
    right = Neumann_inh()
    )

    BC_int = [FreeSurface(), WallNoSlip(θe = phys.theta_e * π / 180, θadv = phys.theta_adv * π / 180, θrec = phys.theta_rec * π / 180)]

elseif sim.name == "sessile_2LS_inclined"

    gp.LS[2].u .= 3.7 .- gp.y

    phL.u .= 0.0
    phL.v .= 0.0

  
    BC_uL = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Navier_cl(λ = 1e-2),
        top = Dirichlet(),
    )

    BC_vL = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Dirichlet(),
        top = Dirichlet()
    )

    BC_pL = Boundaries(
        left = Periodic(),
        right = Periodic(),
    )

    BC_u = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Neumann_cl(θe = π / 2.0),
        top = Neumann_inh(),
    )

    BC_int = [FreeSurface(), WallNoSlip()]
  

    elseif sim.name == "levelset_Butler"

        BC_u = Boundaries(
        bottom = Neumann_inh(),
        top = Neumann_inh(),
        left = Neumann_inh(),
        right = Neumann_inh())



        BC_pL = Boundaries()

    elseif sim.name == "levelset_Butler_3"
        
        BC_vL = Boundaries()
        BC_pL = Boundaries()


    elseif sim.name == "levelset_Butler_Navier_slip"


        BC_vL = Boundaries(left=Navier_cl(λ = phys.Navier_slip_length))


          # BC_u = Boundaries(
        # bottom = Neumann_cl(θe = _θe * π / 180),
        # top = Neumann_inh(),
        # left = Neumann_inh(),
        # right = Neumann_inh())


        # BC_u = Boundaries(
        # bottom = Neumann_inh(),
        # top = Neumann_inh(),
        # left = Neumann_inh(),
        # right = Neumann_inh())

        #gp.LS[2].u .= gp.x


        # BC_uL = Boundaries(
        #     bottom = Navier_cl(λ = 1e-2),
        #     top = Dirichlet(),
        # ),
        # BC_vL = Boundaries(
        #     bottom = Dirichlet(),
        #     top = Dirichlet(),
        # ),
        # BC_pL = Boundaries(
        #     left = Dirichlet(),
        #     right = Dirichlet(),
        # ),
        # BC_u = Boundaries(
        #     bottom = Neumann_cl(θe = _θe * π / 180),
        #     top = Neumann_inh(),
        #     left = Neumann_inh(),
        #     right = Neumann_inh()
        #),
        #BC_int = [FreeSurface()], keep wall with u=0 and Neumann pressure

        BC_pL = Boundaries()

        # BC_uL = Boundaries(bottom = Navier_cl(λ = 1e-2),),
        # BC_vL = Boundaries(bottom = Dirichlet(),),
        # BC_pL = Boundaries(),
        # BC_u = Boundaries(
        #     bottom = Neumann_cl(θe = θe2 * π / 180),
        #     top = Neumann(),
        #     left = Neumann_inh(),
        #     right = Neumann_inh()
        # ),
        # BC_int = [FreeSurface(), WallNoSlip(θe = θe * π / 180)],

    elseif sim.name == "levelset_Butler_Dirichlet"

        BC_vL = Boundaries(left=Dirichlet(val=gv.y[:,1]))
        BC_pL = Boundaries()

end

#Easier for sim.CFL: one velocity, and not phys.v_inlet, 3phys.v_inlet/2
vPoiseuille = Poiseuille_fmax.(gv.x,phys.v_inlet,phys.ref_length) 
vPoiseuilleb = Poiseuille_fmax.(gv.x[1,:],phys.v_inlet,phys.ref_length) 


if sim.imposed_velocity == "radial"
    phL.v .= 0.0
    phL.u .= 0.0
    phS.v .= 0.0
    phS.u .= 0.0

    radial_vel!(phL,gu,gv,radial_vel_factor,phys.intfc_x,phys.intfc_y,phys.radius)
    radial_vel!(phS,gu,gv,radial_vel_factor,phys.intfc_x,phys.intfc_y,phys.radius)

    printstyled(color=:red, @sprintf "\n radial vel: %.2e %.2e sim.CFL %.2e sim.dt0 %.2e dx %.2e \n" maximum(phL.u) maximum(phL.v) max(maximum(phL.u),maximum(phL.v))*sim.dt0/gp.dx[1,1] sim.dt0 gp.dx[1,1])
    vPoiseuilleb = 0.0

elseif sim.imposed_velocity == "constant"
    
    phL.v .= phys.v_inlet
    phL.u .= 0.0

elseif sim.imposed_velocity == "Poiseuille_bottom_top"

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(val = vPoiseuilleb),
        top    = Dirichlet(val = vPoiseuilleb),
    )
    phL.v .=vPoiseuille 
    phL.u .= 0.0
    printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))

elseif sim.imposed_velocity == "Poiseuille_bottom"

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(val = vPoiseuilleb),
        top    = Neumann(),
    )


    # BC_pL= Boundaries(
    #     left   = Dirichlet(),
    #     right  = Dirichlet(),
    #     bottom = Neumann(),
    #     top    = Dirichlet(),
    # )

    BC_pL= Boundaries(
        left   = Neumann(),
        right  = Neumann(),
        bottom = Neumann(),
        top    = Dirichlet(),
    )

    phL.v .=vPoiseuille 
    phL.u .= 0.0
    printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))

    
elseif sim.imposed_velocity == "Poiseuille_pressure_only"

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Neumann(),
        top    = Neumann(),
    )


    # BC_pL= Boundaries(
    #     left   = Dirichlet(),
    #     right  = Dirichlet(),
    #     bottom = Neumann(),
    #     top    = Dirichlet(),
    # )

    BC_pL= Boundaries(
        left   = Neumann(),
        right  = Neumann(),
        bottom = Dirichlet(val = p_bottom),
        top    = Dirichlet(),
    )

    phL.v .= 0.0
    phL.u .= 0.0
    printstyled(color=:red, @sprintf "\n max u %.2e v %.2e p %.2e \n" maximum(phL.u) maximum(phL.v) maximum(phL.p) )
    print("p_bottom ", p_bottom, "\n")
    print(BC_pL)

elseif sim.imposed_velocity == "zero"

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Dirichlet(),
    )

    BC_uL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Dirichlet(),
    )

    BC_pL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Dirichlet(),
    )

    phL.v .= 0.0
    phL.u .= 0.0
    printstyled(color=:red, @sprintf "\n max u %.2e v %.2e p %.2e \n" maximum(phL.u) maximum(phL.v) maximum(phL.p) )
    print("p_bottom ", p_bottom)
    print(BC_pL)
# else

#     phL.v .=vPoiseuille 
#     phL.u .= 0.0

    
#     printstyled(color=:red, @sprintf "\n Poiseuille BC Dir + Neu \n")
        
#     if sim.name != "channel_Dirichlet_pressure"
#         BC_vL= Boundaries(
#         left   = Dirichlet(),
#         right  = Dirichlet(),
#         bottom = Dirichlet(val = copy(vPoiseuilleb)),
#         top    = Neumann(val=0.0),
#         )
#     end
    
#     printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))
end

try
    if phys.boundary_conditions == "dir"
        global BC_trans_scal_H2O = BoundariesInt(
            bottom = Dirichlet(val = phys.concentration0[3]),
            top    = Neumann(),
            left   = Neumann(val=i_butler/(phys.Faraday*DH2O)),
            right  = Dirichlet(val = phys.concentration0[3]),
            int    = Dirichlet(val = phys.concentration0[3]))
            print(BC_trans_scal_H2O)
    end
catch error
    printstyled(color=:red, @sprintf "\n not modifying BC \n")
    # print(error)
    print(BC_trans_scal_H2O)
end

printstyled(color=:green, @sprintf "\n Initialisation0 \n")

print_electrolysis_statistics(num,gp,phL)



if pressure_channel

    # test_Poiseuille(num,phL,gv)

    # vecb_B(phL.vD,gv) .= vPoiseuilleb
    # vecb_T(phL.vD,gv) .= vPoiseuilleb

    # test_Poiseuille(num,phL,gv)

    
    phL.p .= gp.y * (p_top - p_bottom)/phys.ref_length .+ p_bottom

    printstyled(color=:red, @sprintf "\n test x %.2e y %.2e p %.2e p_bottom %.2e \n" gp.x[1,1] gp.y[1,1] phL.p[1,1] p_bottom)


    # p_bottom .+ (p_top .- p_bottom)./phys.ref_length.*x

    printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

    # phL.p[1,:] .= p_bottom 

    # printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))


    printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

    printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

else
    phL.p .= phys.pres0 #0.0

end

phL.T .= phys.temperature0
phS.T .= phys.temperature0

#Initialize for CRASH detection (cf isnan...)
phS.v .= 0.0
phS.u .= 0.0
phS.vD .= 0.0
phS.uD .= 0.0


#Test array order Julia ->PDI
# for j in 1:gp.ny
#     for i in 1:gp.nx
#         pdi_array[j,i]=1000*i+j
#     end
# end
    


if num.io_pdi>0

    @debug "Before PDI init"
        
    # # using MPI
    # MPI.Init()

    # @debug "after MPI.Init"

    # comm = MPI.COMM_WORLD
    # @debug "MPI.COMM_WORLD"

    # print(comm)

    yml_file = yamlfile

    # print("\n yml_file ",yml_file)

    # Version: julia +1.10.4

    conf = @ccall "libparaconf".PC_parse_path(yml_file::Cstring)::PC_tree_t

    @debug "after conf"

    getsubyml = @ccall "libparaconf".PC_get(conf::PC_tree_t,".pdi"::Cstring)::PC_tree_t  

    @debug "after getsubyml"


    # print("\n getsubyml ",getsubyml)

    # @ccall "libpdidummy".PDI_init(getsubyml::PC_tree_t)::Cvoid
    local pdi_status = @ccall "libpdi".PDI_init(getsubyml::PC_tree_t)::Cint

    # print("\n pdi_status ",pdi_status)


    # print(getsubyml)
    # pdi_status = @ccall "libpdi".PDI_init(getsubyml::PC_tree_t)::Cint

    @debug "after PDI_init"


    # print("\n PDI_init ")

    # @ccall "libpdi".PDI_init(conf::PC_tree_t)::Cvoid

    # #python event to plot
    # @ccall "libpdi".PDI_event("testing"::Cstring)::Cvoid

    # Send meta-data to PDI
    mpi_coords_x = 1
    mpi_coords_y = 1
    mpi_max_coords_x = 1
    mpi_max_coords_y = 1

    nx=gp.nx
    ny=gp.ny

    #TODO check Clonglong ...

    phys_time = 0.0 #Cdouble
    nstep = num.current_i

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

end #if num.io_pdi>0



# current_t = 0
# num.current_i = 0

# #PDI (IO)


# if num.io_pdi>0
#     try
#         printstyled(color=:red, @sprintf "\n PDI test \n" )

#         phys_time = current_t #Cdouble
#         nstep = num.current_i
#         # print("\n nstep ",typeof(nstep))
#         # pdi_array =zeros(nx,ny)

#         # for j in 1:gp.ny
#         #     for i in 1:gp.nx
#         #         pdi_array[j,i]=1000*i+j
#         #     end
#         # end

#         vec1(phL.vD,gv) .= vec(phL.v)
#         # vec2(phL.vD,gv) .= 333 #test

#         # phi_array=phL.phi_ele #do not transpose since python row major
        
#         compute_grad_phi_ele!(num, gp, gu, gv, phL, phS, op.opC_pL, op.opC_pS) #TODO current

#         Eus,Evs = interpolate_grid_liquid(gp,gu,gv,phL.Eu, phL.Ev)

#         us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)


#         print("\n before write \n ")
#         #if writing "D" array (bulk, interface, border), add "_1D" to the name
#         @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
#         "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
#         "time"::Cstring, phys_time::Ref{Cdouble}, PDI_OUT::Cint,
#         "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
#         "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
#         "trans_scal_1D"::Cstring, phL.trans_scalD::Ptr{Cdouble}, PDI_OUT::Cint,
#         "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
#         "i_current_x"::Cstring, Eus::Ptr{Cdouble}, PDI_OUT::Cint,   
#         "i_current_y"::Cstring, Evs::Ptr{Cdouble}, PDI_OUT::Cint,   
#         "velocity_x"::Cstring, us::Ptr{Cdouble}, PDI_OUT::Cint,   
#         "velocity_y"::Cstring, vs::Ptr{Cdouble}, PDI_OUT::Cint,   

#         C_NULL::Ptr{Cvoid})::Cvoid

#         print("\n after write \n ")

#         @ccall "libpdi".PDI_finalize()::Cvoid

#         printstyled(color=:red, @sprintf "\n PDI test end\n" )

#     catch error
#         printstyled(color=:red, @sprintf "\n PDI error \n")
#         print(error)
#         printstyled(color=:red, @sprintf "\n PDI error \n")
#     end
# end #if io.pdi>0

# 




# r = 0.5
# gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ phys.ref_lengthy / 2).^2) - r * ones(gp)
# gp.LS[1].u .*= -1.0


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



# x,y,field,connectivities,num_vtx = convert_interfacial_D_to_segments(num,gp,phL.T,1)
# print("\n number of interface points ", num_vtx)
# # print("\n x",x)
# # print("\n x",y)
# # print("\n x",field)
# print("\n x",connectivities)
# print("\n x",num_vtx)


vecb_L(phL.uD, gu) .= 0.0
vecb_B(phL.uD, gu) .= 0.0
vecb_R(phL.uD, gu) .= 0.0
vecb_T(phL.uD, gu) .= 0.0


printstyled(color=:green, @sprintf "\n sim.CFL : %.2e dt : %.2e\n" sim.CFL sim.CFL*phys.ref_length/mesh.nx/phys.v_inlet)


for iscal=1:phys.nb_transported_scalars
    phL.trans_scal[:,:,iscal] .= phys.concentration0[iscal]
end

phL.phi_ele .= phys.phi_ele0

printstyled(color=:green, @sprintf "\n Initialisation \n")

# print_electrolysis_statistics(phys.nb_transported_scalars,gp,phL)
print_electrolysis_statistics(num,gp,phL)

printstyled(color=:green, @sprintf "\n TODO timestep sim.CFL scal, and print \n")


@unpack τ,CFL,Δ,Re,θd=num
# print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ sim.CFL sim.CFL*Δ sim.CFL*Δ^2*Re Re θd)
# τ=sim.CFL*Δ/phys.v_inlet
# num.τ=τ
# print(@sprintf "dt %.2e \n" τ)

#TODO pressure left and right BC not mentioned in the article Khalighi 2023

#TODO need to iterate more for potential since phiele=0 initially?

# print("\n linrange ", x[1,:] .*0.0)

# print("\n phi_ele ", phi_ele)


# eta = phys.phi_ele1 .-phi_ele
#TODO precision: number of digits
# i_butler=phys.i0*(exp(phys.alpha_a*phys.Faraday*eta/(phys.Ru*phys.temperature0))-exp(-phys.alpha_c*phys.Faraday*eta/(phys.Ru*phys.temperature0)))
# i_butler=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele,phys.phi_ele1,phys.Ru,phys.temperature0)




BC_uS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

BC_vS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

BC_pS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

if sim.time_scheme == "FE"
    time_scheme = FE
else
    time_scheme = CN
end

# @debug "Before run"
# print("\n Before run")
# print(i_butler./elec_cond," ",size(i_butler./elec_cond),"\n ")

if phys.nb_levelsets ==1
    BC_int = [WallNoSlip()] #[FreeSurface()]
end

if sim.name == "falling_drop" || occursin("sessile",sim.name)
    BC_trans_scal = Vector{BoundariesInt}() #[Boundaries()]
    BC_phi_ele = Boundaries()
else
    BC_trans_scal = [
        BC_trans_scal_H2, #H2
        BC_trans_scal_KOH, #KOH
        BC_trans_scal_H2O, #H2O       
    ] #)

   
    # left = Dirichlet(val=phys.pres0),
    # right = Dirichlet(val=phys.pres0),

    # BC_u = Boundaries(
    #     bottom = Neumann_inh(),
    #     top = Neumann_inh(),
    #     left = Neumann_cl(θe = _θe * π / 180),
    #     right = Neumann_inh()
    # )
    BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_inh(),
    right = Neumann_inh())
end




if sim.bc_int == "FreeSurface"
    BC_int = [FreeSurface()]
end

if num.io_pdi>0

    try
        printstyled(color=:red, @sprintf "\n before pdi \n")

        # printstyled(color=:red, @sprintf "\n PDI test \n" )


   
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

        # compute_grad_phi_ele!(num, grid, grid_u, grid_v, phL, phS, op.opC_pL, op.opC_pS, elec_cond,tmp_vec_u,tmp_vec_v,tmp_vec_p,tmp_vec_p1,tmp_vec_p1) #TODO current

        
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


        PDI_status = @ccall "libpdi".PDI_multi_expose("write_initialization"::Cstring,
        "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
        "time"::Cstring, phys_time::Ref{Cdouble}, PDI_OUT::Cint,
        "u_1D"::Cstring, phL.uD::Ptr{Cdouble}, PDI_OUT::Cint,
        "v_1D"::Cstring, phL.vD::Ptr{Cdouble}, PDI_OUT::Cint,
        "p_1D"::Cstring, phL.pD::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p"::Cstring, gp.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_u"::Cstring, gu.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_v"::Cstring, gv.LS[iLSpdi].u::Ptr{Cdouble}, PDI_OUT::Cint,
        "levelset_p_wall"::Cstring, LStable::Ptr{Cdouble}, PDI_OUT::Cint,
        # "trans_scal_1DT"::Cstring, phL.trans_scalD'::Ptr{Cdouble}, PDI_OUT::Cint,
        # "phi_ele_1D"::Cstring, phL.phi_eleD::Ptr{Cdouble}, PDI_OUT::Cint,   
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
printstyled(color=:red, @sprintf "\n after pdi \n")


print("\n BC_u ",BC_u)
print("\n BC_int ",BC_int)
print("\n BC_uL ",BC_uL)
print("\n BC_phi_ele ",BC_phi_ele)
print("\n BC_phi_ele ",BC_phi_ele.LS[1])
print("\n BC_phi_ele ",BC_phi_ele.LS[1].val)

print("\n BC_vL ",BC_vL)

print("\n BC_vL ",BC_vL.bottom.val)

BC_vL


print("\n BC_trans_scal ",BC_trans_scal)



printstyled(color=:red, @sprintf "\n before run_forward \n")

run_forward!(
# tinf = @snoop_inference run_forward!(
# tinf = @snoopi_deep run_forward!(
# @profile @time current_i=run_forward(
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
    # BC_TL  = Boundaries(
    # left   = Dirichlet(val = phys.temperature0),
    # right  = Dirichlet(val = phys.temperature0),
    # bottom = Dirichlet(val = phys.temperature0),
    # top    = Dirichlet(val = phys.temperature0),
    # ),
    BC_trans_scal=BC_trans_scal,
    # BC_trans_scal = (
    #     BC_trans_scal_H2, #H2
    #     BC_trans_scal_KOH, #KOH
    #     BC_trans_scal_H2O, #H2O       
    # ),
    BC_phi_ele = BC_phi_ele,
    auto_reinit = sim.auto_reinit,
    # save_length = true,
    time_scheme = time_scheme,
    electrolysis = true,
    navier_stokes = true,
    ns_advection = ns_advection,
    ns_liquid_phase = true,
    Vmean = (sim.mean_velocity == 1),
    verbose = true,
    show_every = sim.show_every,
    electrolysis_convection = true,  
    electrolysis_liquid_phase = true,
    electrolysis_phase_change_case = sim.electrolysis_phase_change_case,
    electrolysis_reaction = phys.electrolysis_reaction, 
    imposed_velocity = sim.imposed_velocity,
    adapt_timestep_mode = sim.adapt_timestep_mode,#1,
    non_dimensionalize=sim.non_dimensionalize,
    mode_2d = sim.mode_2d,
    breakup = sim.breakup,    
)

# @show tinf


@debug "After run"

printstyled(color=:green, @sprintf "\n max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(phL.u)) maximum(abs.(phL.v)))


if occursin("sessile",sim.name)

    V0 = 0.5 * π * 0.5^2
    Vf = volume(gp.LS[1].geoL)
    Vratio = Vf / V0

    mean_rad = 1 / abs(mean(gp.LS[1].κ[gp.LS[1].MIXED[5:end-5]]))
    RR0_sim = mean_rad / 0.5
    RR0_teo = RR0(phys.theta_e * π / 180)

    println("Vratio = $(Vratio)")
    println("mean rad = $(mean_rad)")
    println("RR0_sim = $(RR0_sim)")
    println("RR0_teo = $(RR0_teo)")

end

# vec = Poiseuille_fmax.(gv.x,phys.v_inlet,phys.ref_length)
# print("\n Poiseuille_fmax ",vec[1,:])
# print("\n ")
# print("\n phL.v ",phL.v[1,:])
# print("\n ")


# print("\n rel err",abs.(phL.v[1,:].-vec[1,:])./phys.v_inlet)


# print("\n ", gv.x[1,:])
# print("\n ", gv.x[1,1])

# print("\n P ",Poiseuille_fmax(gv.x[1,1],phys.v_inlet,phys.ref_length)," v ",phL.v[1,1]," test ",Poiseuille_fmax(0.0,phys.v_inlet,phys.ref_length))


#Attention = vs copy

# printstyled(color=:red, @sprintf "\n iMx %.10e Mx %.10e eps %.10e \n" 1/(1e-4/128)^2 (1e-4/128)^2 eps(0.01))

print("\n mesh.nx ",mesh.nx," nx ",gp.nx," ny ",gp.ny)



if sim.name == "channel_Dirichlet_pressure"
    printstyled(color=:green, @sprintf "\n p_bottom : %.10e p_top %.10e grad %.10e grad %.10e \n" p_bottom p_top -(p_top-p_bottom)/phys.ref_length 8*mu1/phys.ref_length^2*phys.v_inlet)
end


if num.io_pdi>0
    try
        local PDI_status = @ccall "libpdi".PDI_finalize()::Cint
        # printstyled(color=:red, @sprintf "\n PDI end\n" )

    catch error
        printstyled(color=:red, @sprintf "\n PDI error \n")
        print(error)
    end
end #if num.io_pdi>0

#Tests 


bubble_volume_0 = 0.5 * π * phys.radius^2

bubble_volume = volume(gp.LS[1].geoS)
equivalent_radius = sqrt(bubble_volume*2/π)

force_buoyancy = - phys.g * (phys.rho1- phys.rho2) * bubble_volume

force_buoyancy_0 = - phys.g * (phys.rho1- phys.rho2) * bubble_volume_0


adim_nb_Bond_final = num.current_radius * phys.g * abs(phys.rho1-phys.rho2)/phys.sigma
adim_nb_Bond_initial = phys.radius * phys.g * abs(phys.rho1-phys.rho2)/phys.sigma


printstyled(color=:cyan, @sprintf "\n Bond number: initial %.2e final %.2e \n" adim_nb_Bond_initial adim_nb_Bond_final)

printstyled(color=:cyan, @sprintf "\n volume: initial %.2e final %.2e R %.2e \n" bubble_volume_0 bubble_volume equivalent_radius)

printstyled(color=:red, @sprintf "\n Test radius R0 %.2e R %.2e \n" phys.radius num.current_radius)

printstyled(color=:cyan, @sprintf "\n Bond number: initial %.2e final %.2e \n" adim_nb_Bond_initial adim_nb_Bond_final)

printstyled(color=:cyan, @sprintf "\n Buoyancy force: initial %.2e final %.2e \n" force_buoyancy_0 force_buoyancy)



# Expression: minimum(phL.phi_eleD) < -1.2e9
# Evaluated: -0.004326033691589959 < -1.2e9


test_tolerance = 1e-14

@testset "Phase change: mass flux" begin

@testset "Phase change: mass flux" begin
    phL.trans_scalD[:,1] .= 1.0 
    mass_flux_vec1 = fzeros(gp)
    mass_flux_vecb = fzeros(gp)
    mass_flux_veci = fzeros(gp)
    mass_flux = zeros(gp)
integrate_mass_flux_over_interface_2_no_writing(num,gp,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,mass_flux)
# @test sum(mass_flux) == 0 
@test sum(mass_flux) ≈ 0 atol=test_tolerance
end #"Phase change: mass flux" begin

@testset "Phase change: mass flux old" begin
    phL.trans_scalD[:,1] .= 1.0 
    mass_flux_vec1 = fzeros(gp)
    mass_flux_vecb = fzeros(gp)
    mass_flux_veci = fzeros(gp)
    mass_flux = zeros(gp)
integrate_mass_flux_over_interface_no_writing(num,gp,op.opC_pL,phL.trans_scalD[:,1],mass_flux_vec1,mass_flux_vecb,mass_flux_veci,mass_flux)
# @test sum(mass_flux) == 0 
@test sum(mass_flux) ≈ 0 atol=test_tolerance
end #"Phase change: mass flux" begin

end #phase change

@testset "Interpolation" begin
gp.V .= 1.0
interpolate_scalar!(gp, gu, gv, gp.V, gu.V, gv.V)

@test minimum(gu.V) ≈ 1.0 atol=test_tolerance
@test maximum(gu.V) ≈ 1.0 atol=test_tolerance
@test minimum(gv.V) ≈ 1.0 atol=test_tolerance
@test maximum(gv.V) ≈ 1.0 atol=test_tolerance

tmp_vec_p = zeros(gp) 
tmp_vec_p0 = zeros(gp) 

tmp_vec_u = zeros(gu) 
tmp_vec_v = zeros(gv) 

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)
# interpolate_grid_solid!(gp,gu,gv,tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# tmp_vec_u .= 1.0 
# tmp_vec_v .= 1.0
# interpolate_grid_liquid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# interpolate_grid_solid_2!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

# @test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
# @test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
# @test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

tmp_vec_u .= 1.0 
tmp_vec_v .= 1.0
interpolate_grid_liquid_solid!(num,gp,gu.LS[end],gv.LS[end],tmp_vec_u, tmp_vec_v,tmp_vec_p,tmp_vec_p0)

@test minimum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p) ≈ 1.0 atol=test_tolerance
@test minimum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance
@test maximum(tmp_vec_p0) ≈ 1.0 atol=test_tolerance

# LS_u =grid_u.LS[1]
# LS_v = grid_v.LS[1]
# us .= (
#     (u[:,2:end] .* LS_u.geoL.dcap[:,2:end,6] .+ 
#     u[:,1:end-1] .* LS_u.geoL.dcap[:,1:end-1,6]) ./ 
#     (LS_u.geoL.dcap[:,1:end-1,6] .+ LS_u.geoL.dcap[:,2:end,6])
# )
# vs .= (
#     (v[2:end,:] .* LS_v.geoL.dcap[2:end,:,7] .+ 
#     v[1:end-1,:] .* LS_v.geoL.dcap[1:end-1,:,7]) ./
#     (LS_v.geoL.dcap[1:end-1,:,7] .+ LS_v.geoL.dcap[2:end,:,7])
# )

# u = phL.Eu
# v = phL.Ev


for j in 1:gp.ny
    for i in 1:gp.nx
        if (tmp_vec_p[j,i] != 1.0) 
            print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
            print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
        end
    end
end

for j in 1:gp.ny
    for i in 1:gp.nx
        if (tmp_vec_p0[j,i] != 1.0) 
            print("\n j",j," i ",i, " tmp_vec_p ",tmp_vec_p[j,i]," tmp_vec_p0 ",tmp_vec_p0[j,i])
            print("\n LS_u.geoL.dcap[:,2:end,6] ",gu.LS[1].geoL.dcap[j,i,6]," ",gu.LS[1].geoL.dcap[j,i+1,6])
        end
    end
end

end

if sim.name == "test_levelset_Butler_two_LS"
    print("\n phi ", minimum(phL.phi_eleD) ," ref ", -1.2e9)
    @testset "Poisson equation with multiple LS" begin
        @test minimum(phL.phi_eleD) < -1.2e9
    end 
else
    @testset "Phase change: radius" begin
        @test num.current_radius > phys.radius
    end #@testset "laplacian"
end


# @test (visc_term - (p_top-p_bottom)/L0)/((p_top-p_bottom)/L0) ≈0 atol=test_tolerance skip=true  
