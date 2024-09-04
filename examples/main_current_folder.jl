# using Revise


# using PrecompileTools_startup

using Flower


using PrettyTables


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


# Sessile from sessile.jl
Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)

θe= 90
θe= 145

if θe < 40
    max_its = 35000
    n_ext = 10
    sim.CFL = 0.5
elseif θe < 100
    max_its = 15000
    n_ext = 10
    sim.CFL = 0.5
else
    max_its = 5000
    n_ext = 10
    sim.CFL = 0.5
end


# Physical parameters 
x = LinRange(mesh.xmin, mesh.xmax, mesh.nx+1)
y = LinRange(mesh.ymin, mesh.ymax, mesh.ny+1)



# Physics
mu = phys.mu_cin1 *phys.rho1 #in Pa s = M L^{-1} T^{-1}}

mu1=mu
mu2=mu 

h0 = phys.radius
ref_thickness_2d = 4.0 / 3.0 *phys.radius 

c0_H2,c0_KOH,c0_H2O = phys.concentration0

DH2,DKOH,DH2O= phys.diffusion_coeff

nb_saved_scalars=1

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

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold cyan")

diffusion_t = (phys.radius^2)./phys.diffusion_coeff

pretty_table(vcat(
    hcat("Diffusion time",diffusion_t'),
    hcat("Diffusion coef",phys.diffusion_coeff'),
    hcat("Concentration",phys.concentration0')); 
formatters    = ft_printf("%0.2e", 2:4), #not ecessary , 2:4
header = ["","H2", "KOH", "H2O"], 
highlighters=hl)

@debug "After Table"


# printstyled(color=:green, @sprintf "\n Species diffusion timescales: %.2e %.2e %.2e \n" (phys.radius^2)/DH2 (phys.radius^2)/DKOH (phys.radius^2)/DH2O )

# printstyled(color=:green, @sprintf "\n nmol : \n" phys.concentration0[1]*4.0/3.0*pi*phys.radius^3  )

current_radius = phys.radius

p_liq= phys.pres0 #+ mean(veci(phL.pD,grid,2)) #TODO here one bubble
# p_g=p_liq + 2 * phys.sigma / current_radius
p_g=p_liq + phys.sigma / current_radius

c0test = p_g / (phys.temperature0 * phys.Ru) 

# printstyled(color=:green, @sprintf "\n c0test: %.2e \n" c0test)

# printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" phys.concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(phys.temperature0*phys.Ru))


# Sessile
# h0 = 0.5
# phys.ref_lengthx = 3.0
# phys.ref_lengthy = 2.0
# mesh.nx = 96

# x = collect(LinRange(-phys.ref_lengthx / 2, phys.ref_lengthx / 2, mesh.nx + 1))
# y = collect(LinRange(-phys.ref_lengthy / 2, 0, mesh.nx ÷ 3 + 1))

# not imposing the angle exactly at the boundary but displaced a cell because ghost cells are not used. 
# So the exact contact angle cannot be imposed
# TODO should at some point modify the levelset so it works as the other fields that we have in the code, 
# with the boundary values in addition to the bulk field. That way we could impose it exactly at the boundary
# needs a lot of work, not priority
_θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# println("θe = $(_θe)")

radial_vel_factor = 1e-7

# BC 
i_butler = x[1,:] .*0.0

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

pressure_channel = false

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

    p_top = 0
    # phys.v_inlet: avg
    # p_bottom = p_top - 8*mu1/phys.ref_length*phys.v_inlet
    # p_bottom = p_top + 8*mu1/phys.ref_length*phys.v_inlet*3/2 #vmoy
    p_bottom = p_top + 8*mu1/phys.ref_length*phys.v_inlet

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

    # sim.max_iter = 2
    # save_every = 1

    # sim.CFL 1
    sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet

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
    # save_every = 25

    save_every = sim.max_iter


    # sim.CFL 1
    sim.dt0 = phys.ref_length/mesh.nx/phys.v_inlet

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
    sim.imposed_velocity = "Poiseuille"

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

    # sim.max_iter = 2
    # save_every = 1
    sim.convection_Cdivu = true

elseif sim.name == "channel_no_bubble_no_vel"
    sim.activate_interface = 0
    sim.max_iter = 100
    save_every = 25

    # sim.max_iter = 2
    # save_every = 1


elseif sim.name == "channel"
    sim.max_iter = 100
    save_every = 25




elseif sim.name == "imposed_radius"
    electrolysis_phase_change = true
    sim.max_iter = 100
    save_every = 25



    sim.electrolysis_phase_change_case = "imposed_radius"

end

# folder="electrolysis_circle_wall_CFL"*"_"*sim.name

folder = sim.name



# Save path
if makedir
    prefix *= "/"*folder*"/"
    isdir(prefix) || mkdir(prefix)
    # yamlfile
    yamlfile2="flower.yml"
    cp(yamlpath,prefix*yamlfile2,force=true) #copy yaml file to simulation directory
end




@debug "Before Numerical"



num = Numerical(
    CFL = sim.CFL,
    Re = Re,
    TEND=phys.end_time,
    x = x,
    y = y,
    xcoord = phys.intfc_x,
    ycoord = phys.intfc_y,
    case = "Cylinder",#"Planar",
    R = phys.radius,
    max_iterations = sim.max_iter,
    save_every = save_every,
    ϵ = sim.epsilon, 
    epsilon_mode = sim.epsilon_mode,
    nLS = phys.nb_levelsets,
    nb_transported_scalars=phys.nb_transported_scalars,
    nb_saved_scalars=nb_saved_scalars,
    concentration0=phys.concentration0, 
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
    mu1=mu,
    mu2=mu,
    rho1=phys.rho1,
    rho2=phys.rho2,
    u_inf = 0.0,
    v_inf = 0.0,
    pres0=phys.pres0,
    # σ=phys.sigma,   
    # reinit_every = 10,
    # nb_reinit = 2,
    # δreinit = 10.0,
    # n_ext_cl = n_ext,
    # NB = 24,
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
    bulk_conductivity = sim.bulk_conductivity
    )

@debug "After Numerical"


#Initialization
gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)


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

elseif sim.imposed_velocity == "Poiseuille"

    BC_vL= Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(val = vPoiseuilleb),
        top    = Dirichlet(val = vPoiseuilleb),
    )
    phL.v .=vPoiseuille 
    phL.u .= 0.0
    printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))

    
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



#PDI attempt (IO)


if io.pdi>0
    try
        @debug "Before PDI init"

        # printstyled(color=:red, @sprintf "\n PDI test \n" )

        # using MPI
        MPI.Init()

        comm = MPI.COMM_WORLD

        @debug "after MPI.Init"


        yml_file = yamlfile

        print("\n yml_file ",yml_file)

        # Version: julia +1.10.4

        conf = @ccall "libparaconf".PC_parse_path(yml_file::Cstring)::PC_tree_t

        @debug "after conf"


        getsubyml = @ccall "libparaconf".PC_get(conf::PC_tree_t,".pdi"::Cstring)::PC_tree_t  

        @debug "after getsubyml"


        # print("\n getsubyml ",getsubyml)
        # @debug "after printing  getsubyml"


        @debug "test dummy"

        @ccall "libpdidummy".PDI_init(getsubyml::PC_tree_t)::Cvoid

        @debug "test dummy"

        # print(getsubyml)
        @ccall "libpdi".PDI_init(getsubyml::PC_tree_t)::Cvoid

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

        @ccall "libpdi".PDI_multi_expose("init_PDI"::Cstring, 
                "mpi_coords_x"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                "mpi_coords_y"::Cstring, mpi_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                "mpi_max_coords_x"::Cstring, mpi_max_coords_x::Ref{Clonglong}, PDI_OUT::Cint,
                "mpi_max_coords_y"::Cstring, mpi_max_coords_y::Ref{Clonglong}, PDI_OUT::Cint,
                "nx"::Cstring, nx::Ref{Clonglong}, PDI_OUT::Cint,
                "ny"::Cstring, ny::Ref{Clonglong}, PDI_OUT::Cint,
                "nb_transported_scalars"::Cstring, phys.nb_transported_scalars::Ref{Clonglong}, PDI_OUT::Cint,
                "nb_levelsets"::Cstring, phys.nb_levelsets::Ref{Clonglong}, PDI_OUT::Cint,
                C_NULL::Ptr{Cvoid})::Cvoid

        @debug "after PDI_multi_expose"

        # print("\n PDI_multi_expose ")


        # time = 0.0
        # nstep = 0
        # pdi_array =zeros(nx,ny)

        # for j in 1:gp.ny
        #     for i in 1:gp.nx
        #         pdi_array[j,i]=1000*i+j
        #     end
        # end

        # print("\n before write \n ")

        # @ccall "libpdi".PDI_multi_expose("write_data"::Cstring,
        #             "nstep"::Cstring, nstep::Ref{Clonglong}, PDI_OUT::Cint,
        #             "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
        #             "main_field"::Cstring, pdi_array::Ptr{Cdouble}, 
        #             PDI_OUT::Cint,
        #             C_NULL::Ptr{Cvoid})::Cvoid

        # print("\n after write \n ")

        # @ccall "libpdi".PDI_finalize()::Cvoid

        # printstyled(color=:red, @sprintf "\n PDI test end\n" )

    catch error
        printstyled(color=:red, @sprintf "\n PDI error \n")
        print(error)
    end
end #if io_pdi

@debug "After full PDI init"


# current_t = 0
# num.current_i = 0
# 
# #PDI (IO)
# 

# if num.io_pdi>0
#     try
#         printstyled(color=:red, @sprintf "\n PDI test \n" )

#         time = current_t #Cdouble
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
#         "time"::Cstring, time::Ref{Cdouble}, PDI_OUT::Cint,
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

print_electrolysis_statistics(phys.nb_transported_scalars,gp,phL)

printstyled(color=:green, @sprintf "\n TODO timestep sim.CFL scal, and print \n")


@unpack τ,CFL,Δ,Re,θd=num
# print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ sim.CFL sim.CFL*Δ sim.CFL*Δ^2*Re Re θd)
# τ=sim.CFL*Δ/phys.v_inlet
# num.τ=τ
# print(@sprintf "dt %.2e \n" τ)

#TODO pressure left and right BC not mentioned in the article Khalighi 2023

#TODO need to iterate more for potential since phiele=0 initially?

phi_ele=gv.x[1,:] .*0.0

# print("\n linrange ", x[1,:] .*0.0)

# print("\n phi_ele ", phi_ele)


# eta = phys.phi_ele1 .-phi_ele
#TODO precision: number of digits
# i_butler=phys.i0*(exp(phys.alpha_a*phys.Faraday*eta/(phys.Ru*phys.temperature0))-exp(-phys.alpha_c*phys.Faraday*eta/(phys.Ru*phys.temperature0)))
i_butler=butler_volmer_no_concentration.(phys.alpha_a,phys.alpha_c,phys.Faraday,phys.i0,phi_ele,phys.phi_ele1,phys.Ru,phys.temperature0)

print(@sprintf "Butler-Volmer %.2e %.2e %.2e %.2e\n" i_butler[1] -i_butler[1]/(2*phys.Faraday*DH2) c0_H2-i_butler[1]/(2*phys.Faraday*DH2)*gp.dx[1,1] c0_H2+i_butler[1]/(2*phys.Faraday*DH2)*gp.dx[1,1])




BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_cl(θe = _θe * π / 180),
    right = Neumann_inh()
)

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

@debug "Before run"

@time current_i=run_forward(
    num, gp, gu, gv, op, phS, phL;
    BC_uL = Boundaries(
        left   = Dirichlet(),#Navier_cl(λ = 1e-2), #Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Neumann(val=0.0),
    ),
    BC_uS=BC_uS,
    BC_vL = BC_vL,
    BC_vS=BC_vS,
    BC_pL = BC_pL,
    BC_pS=BC_pS,
    # left = Dirichlet(val=phys.pres0),
    # right = Dirichlet(val=phys.pres0),

    # BC_u = BC_u,
    # BC_int = [FreeSurface()],
    BC_int = [WallNoSlip()],


    # BC_TL  = Boundaries(
    # left   = Dirichlet(val = phys.temperature0),
    # right  = Dirichlet(val = phys.temperature0),
    # bottom = Dirichlet(val = phys.temperature0),
    # top    = Dirichlet(val = phys.temperature0),
    # ),

    BC_trans_scal = (
        BC_trans_scal_H2, #H2
        BC_trans_scal_KOH, #KOH
        BC_trans_scal_H2O, #H2O       
    ),

    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=i_butler/elec_cond), #TODO -BC in Flower ? so i_butler not -i_butler
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
    ),
     # auto_reinit = true,
    # save_length = true,
    time_scheme = time_scheme,
    electrolysis = true,
    navier_stokes = true,
    ns_advection=true,#false,
    ns_liquid_phase = true,
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
    convection_Cdivu = sim.convection_Cdivu,
    
    # ns_advection = false,
    # auto_reinit = true,
    # # ns_advection = false, #?
    # save_length = true,
)

@debug "After run"


printstyled(color=:green, @sprintf "\n max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(phL.u)) maximum(abs.(phL.v)))


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



if io.write_h5>0
    #test HDF5
    print("\n current_i ", current_i)
    current_i = 2
    striter = @sprintf "%.5i" current_i


    filename="Mx"
    file = prefix*filename*"_"*striter*".h5"
    A = zeros(gv.ny, gv.nx+1)
    for jplot in 1:gv.ny
        for iplot in 1:gv.nx+1
        II = CartesianIndex(jplot, iplot) #(id_y, id_x)
        pII = lexicographic(II, gp.ny + 1)
        A[jplot,iplot] =  1 ./ op.opC_vL.iMx.diag[pII]
        end
    end


    print("\n size A ",size(A))
    # from_jl_p =  PermutedDimsArray(A, (2,1))
    # print("\n A ", A )
    # hf = h5write(file, "data", A)
    # A=transpose(A)
    hf = h5write2(file, "data", A,"w")


    filename="p"
    file = prefix*filename*"_"*striter*".h5"

    # A = phL.p
    # # A=transpose(A)
    # hf = h5write2(file, "data", A,"w")

    filename="v"
    file = prefix*filename*"_"*striter*".h5"
    A = phL.v
    # from_jl_p =  PermutedDimsArray(A, (2,1))
    # print("\n A ", A )
    # hf = h5write(file, "data", A)
    # A=transpose(A)
    hf = h5write2(file, "data", A,"w")

    # hf = h5write2(file, "data", from_jl_p,"w")


    if prediction == 0
        str_prediction = "prediction_Flower"
    elseif prediction == 3
        str_prediction = "prediction_PIII"
    else 
        str_prediction = @sprintf "_%.1i" prediction
    end

    filename="v_err"*"_"*str_prediction*"_"
    file = prefix*filename*"_"*striter*".h5"
    A .= (A .- Poiseuille_fmax.(gv.x,phys.v_inlet,phys.ref_length)) ./phys.v_inlet
    # A=transpose(A)
    h5write2(file, "data", A,"w")

end #io.write_h5>0


# printstyled(color=:red, @sprintf "\n iMx %.10e Mx %.10e eps %.10e \n" 1/(1e-4/128)^2 (1e-4/128)^2 eps(0.01))

print("\n mesh.nx ",mesh.nx," nx ",gp.nx," ny ",gp.ny)



if sim.name == "channel_Dirichlet_pressure"
    printstyled(color=:green, @sprintf "\n p_bottom : %.10e p_top %.10e grad %.10e grad %.10e \n" p_bottom p_top -(p_top-p_bottom)/phys.ref_length 8*mu1/phys.ref_length^2*phys.v_inlet)
end


if num.io_pdi>0
    try
        @ccall "libpdi".PDI_finalize()::Cvoid
        # printstyled(color=:red, @sprintf "\n PDI end\n" )

    catch error
        printstyled(color=:red, @sprintf "\n PDI error \n")
        print(error)
    end
end #if num.io_pdi>0