using Revise
using Flower
using PrettyTables
# using Printf
# using Interpolations
# using PyCall
# using PyPlot

# prefix="/local/home/pr277828/flower/"
# prefixyaml=prefix*"/Flower.jl/examples/"

prefix = "."
prefixyaml = prefix

####################################################################################################
# Sessile from sessile.jl
####################################################################################################
Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)
####################################################################################################


# data = YAML.load_file(prefixyaml*"flower.yml")
# println(data)

####################################################################################################
θe= 90
θe= 145

if θe < 40
    max_its = 35000
    n_ext = 10
    CFL = 0.5
elseif θe < 100
    max_its = 15000
    n_ext = 10
    CFL = 0.5
else
    max_its = 5000
    n_ext = 10
    CFL = 0.5
end
####################################################################################################



test_case = "small_cell"
test_case = "100it"
test_case = "radial"
test_case = "channel_no_bubble"
# test_case = "channel_no_bubble_no_vel"
test_case = "channel_no_bubble_Cdivu"

test_case = "channel_Dirichlet"

test_case = "channel_Dirichlet_constant_vel"




#################################################
#Validation
#################################################

# 1
test_case = "channel_Dirichlet_constant_vel"

# 2 
test_case = "channel_Dirichlet_imposed_Poiseuille"

# 3
test_case = "channel_Dirichlet_pressure"
#


function run_case(test_case,n,max_iter,prefix,prediction,test_tolerance)

    ####################################################################################################
    # Classical arguments


    # max_its=10
    # save_every = max_its÷100
    save_every = 1

   
    L0 = 1e-4 

    # pres0 = 1e5
    pres0 = 0.0

    adapt_timestep_mode=2
    adapt_timestep_mode=3

    dt0 = 5e-5
    dt0 = 2.5e-5 #CFL >1 for n = 128
    dt0 = 1.25e-5 #CFL >0.5 for n = 128
    dt0 = 6.25e-6 #CFL >0.5 for n = 128
    dt0 = 3.125e-6 #radius CFL: 3.54e-01 -7.17e+00% 
    # dt0 = 1e-6 #error dnH2 -2.40e-15
    # dt0 = 1e-6 #dnH2 -2.40e-15


    ###################################################################################################
    # Physical parameters 
    ###################################################################################################

    x = LinRange(0, L0, n+1)
    y = LinRange(0, L0, n+1)


    ####################################################################################################

    rho1=1258.0 #liquid
    #TODO need for 80°C, check with other study
    rho2=0.069575 #"0.7016E-01" in \citet{cohnTABLETHERMODYNAMICPROPERTIES1963} H2
    #Linear interpolation between 350 and 360
    # 350	360	353		B	0.13841	353	350	360
    # 7.02E-02	6.82E-02		-0.00194999999999999	A	-0.000194999999999999	0.069575	0.07016	0.06821



    radius = 3.0e-6 
  


    epsilon = 0.05

    epsilon_mode = 2



    h0 = radius


    ref_thickness_2d = 4.0 / 3.0 *radius 

    # mode_2d = 1 #use equivalent cylinder
    # mode_2d = 2 #mol/meter
    mode_2d = 3

    xcoord = 0.0
    ycoord = L0/2.0
    # xcoord = -xcoord
    # ycoord = -ycoord

    mu_cin=6.7e-7 #m^2/s
    mu = mu_cin *rho1 #in Pa s = M L^{-1} T^{-1}}
    mu1=mu
    mu2=mu #TODO
    i0=1.0
    temperature0=353.0
    pres0= 0.0 #1e5
    sigma=7.7e-2
    KOHwtpercent=30
    phi_ele1=-0.6
    alpha_c=0.5
    alpha_a=0.5
    DH2=5.8e-9
    DKOH=3.2e-9
    DH2O=3.2e-9
    c0_H2=1.6e-1
    c0_KOH=6.7e3
    c0_H2O=4.9e4
    Ru=8.314
    Faraday = 9.64853321233100184e4 #C⋅mol−1
    MWH2 = 2.01568e-3 #kg/mol
    phi_ele0=0.0
    CFL= 0.5 #0.01 #0.5
    Re=1.0
    TEND=7.3#s
    plot_xscale = 1e-6
    concentration_check_factor = 1e-3 # 1e-1 #1e-2

    # elec_cond=1 #TODO
    elec_cond=2*Faraday^2*c0_KOH*DKOH/(Ru*temperature0)


    v_inlet=6.7e-4
    Re=rho1*v_inlet*L0/mu
    printstyled(color=:green, @sprintf "\n Re : %.2e %.2e %.2e %.2e\n" Re rho1/mu1 rho1 mu1)

    #To give in Flower
    Re=rho1/mu1

    printstyled(color=:green, @sprintf "\n Re : %.2e %.2e %.2e %.2e\n" Re rho1/mu1 rho1 mu1)





    #Concentration
    diffusion_coeff=[DH2, DKOH, DH2O]
    concentration0=[0.16, 6700, 49000]
    nb_transported_scalars=3

    # nb_saved_scalars=4 
    # nb_saved_scalars=5
    # nb_saved_scalars=6

    nb_saved_scalars=1


    if length(concentration0)!=nb_transported_scalars
        print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
        @error ("nb_transported_scalars")
    end

    if length(diffusion_coeff)!=nb_transported_scalars
        print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
        @error ("nb_transported_scalars")
    end

    print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)

    # pretty_table(concentration0'; header = ["cH2", "cKOH", "cH2O"])
    # pretty_table(diffusion_coeff'; header = ["DH2", "DKOH", "DH2O"])
    # pretty_table(vcat(hcat("D",diffusion_coeff'),hcat("c",concentration0')); header = ["","H2", "KOH", "H2O"])
    # hl = Highlighter((d,i,j)->d[i,j][1]*d[i,j][2] < 0, crayon"red")

    hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold cyan")

    diffusion_t = (radius^2)./diffusion_coeff

    pretty_table(vcat(
        hcat("Diffusion time",diffusion_t'),
        hcat("Diffusion coef",diffusion_coeff'),
        hcat("Concentration",concentration0')); 
    formatters    = ft_printf("%0.2e", 2:4), #not ecessary , 2:4
    header = ["","H2", "KOH", "H2O"], 
    highlighters=hl)



    # printstyled(color=:green, @sprintf "\n Species diffusion timescales: %.2e %.2e %.2e \n" (radius^2)/DH2 (radius^2)/DKOH (radius^2)/DH2O )


    # printstyled(color=:green, @sprintf "\n nmol : \n" concentration0[1]*4.0/3.0*pi*radius^3  )

    current_radius = radius

    p_liq= pres0 #+ mean(veci(phL.pD,grid,2)) #TODO here one bubble
    # p_g=p_liq + 2 * sigma / current_radius
    p_g=p_liq + sigma / current_radius


    c0test = p_g / (temperature0 * Ru) 


    printstyled(color=:green, @sprintf "\n c0test: %.2e \n" c0test)

    # printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(temperature0*num.Ru))




    ####################################################################################################

    ####################################################################################################
    # not imposing the angle exactly at the boundary but displaced a cell because ghost cells are not used. 
    # So the exact contact angle cannot be imposed
    # TODO should at some point modify the levelset so it works as the other fields that we have in the code, 
    # with the boundary values in addition to the bulk field. That way we could impose it exactly at the boundary
    # needs a lot of work, not priority
    _θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
    # _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
    println("θe = $(_θe)")

    electrolysis_phase_change_case = "None"

    save_u = false
    save_v = false
    save_p = false
    save_zoom = false
    save_big_picture = false 
    save_interface = false 




    imposed_velocity = "none"

    electrolysis_phase_change_case = "Khalighi"

 


    null_space = 0
    # null_space = 1

    advection =0
    advection =1

    time_scheme = FE
    # time_scheme = CN

    testing = "scalar"
    #################################################

    post_processing_python = false
    # post_processing_python = true

    # test_case = "channel"
    # test_case = "imposed_radius"

    # localARGS = isdefined(:newARGS) ? newARGS : ARGS
    # @show localARGS

    # print("\n Arguments ", localARGS)

    # if length(localARGS)>1
    #    test_case = localARGS[1]
    #    printstyled(color=:magenta, @sprintf "\n Test case ")
    #    print(test_case)
    #    print("\n")

    # end

    radial_vel_factor = 1e-7


    activate_interface = false
    activate_interface = true

    save_KOH = true
    save_H2O = true


    fontsize = 2 #2
    printmode = "val"
    plotbc=true


    scalar_debug = false
    # scalar_debug = true

    convection_Cdivu = false


    i_current = x[1,:] .*0.0


    # H2 boundary condition
    BC_trans_scal_H2 = BoundariesInt(
    bottom = Dirichlet(val = concentration0[1]),
    top    = Neumann(),
    left   = Neumann(val=-i_current/(2*Faraday*DH2)), #Dirichlet(val = concentration0[1]), #
    right  = Dirichlet(val = concentration0[1]),
    int    = Dirichlet(val = concentration0[1]))

    BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=-i_current/(2*Faraday*DKOH)),
        right  = Dirichlet(val = concentration0[2]),
        int    = Neumann(val=0.0)) #KOH
        
    BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=i_current/(Faraday*DH2O)),
        right  = Dirichlet(val = concentration0[3]),
        int    = Neumann(val=0.0)) #Dirichlet(val = concentration0[3])),#Neumann(val=0.0)) 
        #H2O #Dirichlet(val = concentration0[3])

    BC_pL = Boundaries(
            left   = Neumann(val=0.0),
            right  = Neumann(val=0.0),
            bottom = Neumann(val=0.0),
            top    = Dirichlet(),
        )



    electrolysis_reaction = "Butler_no_concentration"

    pressure_channel = false

    test_laplacian = false











    #####################################################################################################


    if test_case == "small_cell"
        #Test case 1: small cells (high concentration at vecb_L)
        electrolysis_phase_change_case = "None"
        max_iter = 1
        save_every = 1

        folder="electrolysis_circle_wall_CFL_small_cell"

    elseif test_case == "100it"
        #Test case 2: scalar without velocity
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 10
        save_p = false

        folder="electrolysis_circle_wall_CFL"*test_case

    elseif test_case == "radial"
        electrolysis_phase_change_case = "None"
        imposed_velocity = "radial"
        max_iter = 1
        save_every = 1
        n = 64
        radial_vel_factor = 1e-7
        folder="electrolysis_circle_wall_CFL"*test_case

    elseif test_case == "channel_no_bubble"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25


        # save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4 #1e-3


    elseif test_case == "channel_Dirichlet"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        # max_iter = 2
        # save_every = 1


        # save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4 #1e-3
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[1]),
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])) #H2

        electrolysis_reaction = "none"

    elseif test_case == "channel_Dirichlet_pressure"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        max_iter = 2
        save_every = 1
        max_iter = 1

        test_laplacian = true


        save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = false

        concentration_check_factor = 1e-4 #1e-3
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[1]),
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])) #H2


        electrolysis_reaction = "none"
        # imposed_velocity = "constant"

        p_top = 0
        # v_inlet: avg
        # p_bottom = p_top - 8*mu1/L0*v_inlet
        # p_bottom = p_top + 8*mu1/L0*v_inlet*3/2 #vmoy
        p_bottom = p_top + 8*mu1/L0*v_inlet

        test_v=-(p_top - p_bottom)/L0 * (L0^2)/4/2/mu

        
        printstyled(color=:green, @sprintf "\n mu1 : %.2e v %.2e vtest %.2e CFL %.2e \n" mu1 v_inlet test_v v_inlet*dt0/(L0/n))

        printstyled(color=:green, @sprintf "\n p_bottom : %.2e p_top %.2e grad %.2e grad %.2e \n" p_bottom p_top -(p_top-p_bottom)/L0 8*mu1/L0^2*v_inlet)


        # print("\n mu1 ",mu1," L0 ", L0)

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

        
    elseif test_case == "channel_Dirichlet_constant_vel"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        # max_iter = 100
        save_every = 25

        save_every = max_iter


        # CFL 1
        dt0 = L0/n/v_inlet

        # CFL 0.5
        dt0 = L0/n/v_inlet/2 


        save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = false

        concentration_check_factor = 1e-4 #1e-3
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[1]),
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[2]),
        right  = Dirichlet(val = concentration0[2]),
        int    = Dirichlet(val = concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[3]),
        right  = Dirichlet(val = concentration0[3]),
        int    = Dirichlet(val = concentration0[3])) #H2O


        electrolysis_reaction = "none"
        imposed_velocity = "constant"

        BC_vL= Boundaries(
            left   = Dirichlet(val = v_inlet),
            right  = Dirichlet(val = v_inlet),
            bottom = Dirichlet(val = v_inlet),
            top    = Neumann(val=0.0),
        )

        # BC_pL = Boundaries(
        #     left   = Neumann(),
        #     right  = Neumann(),
        #     bottom = Neumann(),
        #     top    = Neumann(),
        # )
        BC_pL = Boundaries(
            left   = Dirichlet(),
            right  = Dirichlet(),
            bottom = Dirichlet(),
            top    = Dirichlet(),
        )


    elseif test_case == "channel_Dirichlet_constant_vel_half-circle"
        activate_interface = true
        electrolysis_phase_change_case = "None"
        save_every = 25

        save_every = max_iter

        # CFL 1
        dt0 = L0/n/v_inlet

        # CFL 0.5
        dt0 = L0/n/v_inlet/2 


        save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = false

        concentration_check_factor = 1e-4 #1e-3
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[1]),
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[2]),
        right  = Dirichlet(val = concentration0[2]),
        int    = Dirichlet(val = concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[3]),
        right  = Dirichlet(val = concentration0[3]),
        int    = Dirichlet(val = concentration0[3])) #H2O


        electrolysis_reaction = "none"
        imposed_velocity = "constant"

        BC_vL= Boundaries(
            left   = Dirichlet(val = v_inlet),
            right  = Dirichlet(val = v_inlet),
            bottom = Dirichlet(val = v_inlet),
            top    = Neumann(val=0.0),
        )

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Neumann(),
            top    = Neumann(),
        )


    elseif test_case == "channel_Dirichlet_imposed_Poiseuille"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        save_every = max_iter

        # max_iter = 1
        # save_every = 1

        # # CFL 1
        # dt0 = L0/n/v_inlet 

        # CFL 0.5
        dt0 = L0/n/v_inlet/2 

        save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = false

        concentration_check_factor = 1e-4 #1e-3
        
        BC_trans_scal_H2 = BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[1]),
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])) #H2

        BC_trans_scal_KOH = BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[2]),
        right  = Dirichlet(val = concentration0[2]),
        int    = Dirichlet(val = concentration0[2])) #KOH

        BC_trans_scal_H2O = BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Dirichlet(val = concentration0[3]),
        right  = Dirichlet(val = concentration0[3]),
        int    = Dirichlet(val = concentration0[3])) #H2O


        electrolysis_reaction = "none"
        imposed_velocity = "Poiseuille_bottom_top"

        BC_pL = Boundaries(
            left   = Neumann(),
            right  = Neumann(),
            bottom = Neumann(),
            top    = Neumann(),
        )

    elseif test_case == "channel_no_bubble_Cdivu"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        # save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4 #1e-3

        convection_Cdivu = true

    elseif test_case == "channel_no_bubble_no_vel"
        activate_interface = false
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        # save_v = true
        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4 #1e-3

    elseif test_case == "channel"
        electrolysis_phase_change_case = "None"
        max_iter = 100
        save_every = 25

        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4

    elseif test_case == "imposed_radius"
        max_iter = 100
        save_every = 25

        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4

        electrolysis_phase_change_case = "imposed_radius"
    
    elseif test_case == "imposed_radius4"
        max_iter = 100
        save_every = 25

        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4

        electrolysis_phase_change_case = "imposed_radius4"

    elseif test_case == "imposed_radius_dir"
        max_iter = 100
        save_every = 25

        save_KOH = false
        save_H2O = false
        save_zoom = true

        concentration_check_factor = 1e-4

        electrolysis_phase_change_case = "imposed_radius"

    end

    folder="electrolysis_circle_wall_CFL"*"_"*test_case


    show_every = max_iter

    ####################################################################################################
    # Save path
    ####################################################################################################
    prefix *= "/"*folder*"/"
    # isdir(prefix) || mkdir(prefix)
    ####################################################################################################






    num = Numerical(
        CFL = CFL,
        Re = Re,
        TEND=TEND,
        x = x,
        y = y,
        xcoord = xcoord,
        ycoord = ycoord,
        case = "Cylinder",#"Planar",
        R = radius,
        max_iterations = max_iter,
        save_every = save_every,#10,#save_every,#10,
        ϵ = epsilon, 
        epsilon_mode = epsilon_mode,
        nb_transported_scalars=nb_transported_scalars,
        nb_saved_scalars=nb_saved_scalars,
        concentration0=concentration0, 
        diffusion_coeff=diffusion_coeff,
        temperature0=temperature0,
        i0=i0,
        phi_ele0=phi_ele0,
        phi_ele1=phi_ele1,
        alpha_c=alpha_c,
        alpha_a=alpha_a,
        Ru=Ru,
        Faraday=Faraday,
        MWH2=MWH2,
        θd=temperature0,
        eps=1e-12,
        mu1=mu,
        mu2=mu,
        rho1=rho1,
        rho2=rho2,
        u_inf = 0.0,
        v_inf = 0.0,
        pres0=pres0,
        # σ=sigma,   
        # reinit_every = 10,
        # nb_reinit = 2,
        # δreinit = 10.0,
        # n_ext_cl = n_ext,
        # NB = 24,
        plot_xscale = plot_xscale,
        plot_prefix = prefix,
        dt0 = dt0,
        concentration_check_factor = concentration_check_factor,
        radial_vel_factor = radial_vel_factor,
        # scalar_debug =scalar_debug,
        v_inlet = v_inlet,
        prediction = prediction,
        null_space = null_space,
        )
        # ref_thickness_2d = ref_thickness_2d,

    #Initialization
    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)


    @unpack x, nx, ny, ind = gv


    # vPoiseuille = Poiseuille_favg.(x,v_inlet,L0) 
    # vPoiseuilleb = Poiseuille_favg.(x[1,:],v_inlet,L0) 

    #Easier for CFL: one velocity, and not v_inlet, 3v_inlet/2
    # vPoiseuille = zeros(n +1, n) #zeros(gv)
    vPoiseuille = Poiseuille_fmax.(gv.x,v_inlet,L0) 
    vPoiseuilleb = Poiseuille_fmax.(gv.x[1,:],v_inlet,L0) 




    if imposed_velocity == "radial"
        phL.v .= 0.0
        phL.u .= 0.0
        phS.v .= 0.0
        phS.u .= 0.0

        for ju in 1:gu.ny
            for iu in 1:gu.nx
                xcell = gu.x[ju,iu]
                ycell = gu.y[ju,iu]

                vec0 = [xcoord, ycoord]
                vec1 = [xcell, ycell]

                vecr = vec1-vec0
                normr = norm(vecr)
                if normr>radius
                    vecr .*= 1.0/normr
                    factor = 1.0/normr 
                    factor *= radial_vel_factor
                    phL.u[ju,iu] = factor * vecr[1]
                    phS.u[ju,iu] = factor * vecr[1]

                end
            end
        end

        for jv in 1:gv.ny
            for iv in 1:gv.nx    
                xcell = gv.x[jv,iv]
                ycell = gv.y[jv,iv]

                vec0 = [xcoord, ycoord]
                vec1 = [xcell, ycell]

                vecr = vec1-vec0
                normr = norm(vecr)

                if normr>radius

                    vecr .*= 1.0/normr
                    factor = 1.0/normr 
                    factor *= radial_vel_factor

                    phL.v[jv,iv] = factor * vecr[2]
                    phS.v[jv,iv] = factor * vecr[2]

                    #print("\n i ",iv," j ",jv," vec",vecr[2]," y ",ycell," ycoord ",ycoord," v ",phL.v[jv,iv])
                end
            end
        end

        printstyled(color=:red, @sprintf "\n radial vel: %.2e %.2e CFL %.2e dt0 %.2e dx %.2e \n" maximum(phL.u) maximum(phL.v) max(maximum(phL.u),maximum(phL.v))*dt0/gp.dx[1,1] dt0 gp.dx[1,1])

                

        vPoiseuilleb = 0.0

    elseif imposed_velocity == "constant"
        
        phL.v .= v_inlet
        phL.u .= 0.0

    elseif imposed_velocity == "Poiseuille_bottom_top"

        BC_vL= Boundaries(
            left   = Dirichlet(),
            right  = Dirichlet(),
            bottom = Dirichlet(val = vPoiseuilleb),
            top    = Dirichlet(val = vPoiseuilleb),
        )
        phL.v .=vPoiseuille 
        phL.u .= 0.0
        printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))

    else

        phL.v .=vPoiseuille 
        phL.u .= 0.0

        
        printstyled(color=:red, @sprintf "\n Poiseuille BC Dir + Neu \n")
            
        if test_case != "channel_Dirichlet_pressure"
            BC_vL= Boundaries(
            left   = Dirichlet(),
            right  = Dirichlet(),
            bottom = Dirichlet(val = copy(vPoiseuilleb)),
            top    = Neumann(val=0.0),
            )
        end
        
        printstyled(color=:red, @sprintf "\n initialized bulk velocity field %.2e \n" maximum(phL.v))
    end


    if pressure_channel

        # test_Poiseuille(num,phL,gv)

        # vecb_B(phL.vD,gv) .= vPoiseuilleb
        # vecb_T(phL.vD,gv) .= vPoiseuilleb

        # test_Poiseuille(num,phL,gv)

        
        phL.p .= gp.y * (p_top - p_bottom)/L0 .+ p_bottom

        printstyled(color=:red, @sprintf "\n test x %.2e y %.2e p %.2e p_bottom %.2e \n" gp.x[1,1] gp.y[1,1] phL.p[1,1] p_bottom)


        # p_bottom .+ (p_top .- p_bottom)./L0.*x

        printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))

        # phL.p[1,:] .= p_bottom 

        # printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[1,:]) maximum(phL.p[1,:]))


        printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" minimum(phL.p[end,:]) maximum(phL.p[end,:]))

        printstyled(color=:red, @sprintf "\n pressure min %.2e max %.2e\n" BC_pL.bottom.val BC_pL.top.val )

    else
        phL.p .= pres0 #0.0

    end

    phL.T .= temperature0
    phS.T .= temperature0

    #Initialize for CRASH detection (cf isnan...)
    phS.v .= 0.0
    phS.u .= 0.0
    phS.vD .= 0.0
    phS.uD .= 0.0


    ####################################################################################################
    # r = 0.5
    # gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
    # gp.LS[1].u .*= -1.0
    ####################################################################################################

    if activate_interface 


        gp.LS[1].u .= sqrt.((gp.x .- xcoord).^2 + (gp.y .- ycoord).^2) - radius * ones(gp)
  
        su = sqrt.((gv.x .- xcoord).^2 .+ (gv.y .- ycoord).^2)
        R1 = radius + 3.0*num.Δ

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


    vecb_L(phL.uD, gu) .= 0.0
    vecb_B(phL.uD, gu) .= 0.0
    vecb_R(phL.uD, gu) .= 0.0
    vecb_T(phL.uD, gu) .= 0.0

    printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL CFL*L0/n/v_inlet)



    for iscal=1:nb_transported_scalars
        phL.trans_scal[:,:,iscal] .= concentration0[iscal]
    end

    phL.phi_ele .= phi_ele0

    printstyled(color=:green, @sprintf "\n Initialisation \n")

    print_electrolysis_statistics(num,gp,phL)

    printstyled(color=:green, @sprintf "\n TODO timestep CFL scal, and print \n")


    @unpack τ,CFL,Δ,Re,θd=num
    # print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ CFL CFL*Δ CFL*Δ^2*Re Re θd)
    # τ=CFL*Δ/v_inlet
    # num.τ=τ
    # print(@sprintf "dt %.2e \n" τ)

    #TODO pressure left and right BC not mentioned in the article Khalighi 2023

    #TODO need to iterate more for potential since phiele=0 initially?

    phi_ele=gv.x[1,:] .*0.0

    # print("\n linrange ", x[1,:] .*0.0)

    # print("\n phi_ele ", phi_ele)


    # eta = phi_ele1 .-phi_ele
    #TODO precision: number of digits
    # i_current=i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    i_current=butler_volmer_no_concentration.(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)

    print(@sprintf "Butler-Volmer %.2e %.2e %.2e %.2e\n" i_current[1] -i_current[1]/(2*Faraday*DH2) c0_H2-i_current[1]/(2*Faraday*DH2)*gp.dx[1,1] c0_H2+i_current[1]/(2*Faraday*DH2)*gp.dx[1,1])



    ####################################################################################################
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

    print("\n before run_forward \n")

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
        # left = Dirichlet(val=pres0),
        # right = Dirichlet(val=pres0),

        # BC_u = BC_u,
        # BC_int = [FreeSurface()],
        BC_int = [WallNoSlip()],


        # BC_TL  = Boundaries(
        # left   = Dirichlet(val = temperature0),
        # right  = Dirichlet(val = temperature0),
        # bottom = Dirichlet(val = temperature0),
        # top    = Dirichlet(val = temperature0),
        # ),

        BC_trans_scal = (
            BC_trans_scal_H2, #H2
            BC_trans_scal_KOH, #KOH
            BC_trans_scal_H2O, #H2O       
        ),

        BC_phi_ele = BoundariesInt(
            left   = Neumann(val=-i_current/elec_cond),
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
        show_every = show_every,
        electrolysis_convection = true,  
        electrolysis_liquid_phase = true,
        # electrolysis_phase_change = electrolysis_phase_change,
        # electrolysis_phase_change_case = "levelset",
        electrolysis_phase_change_case = electrolysis_phase_change_case,
        electrolysis_reaction = electrolysis_reaction, 
        imposed_velocity = imposed_velocity,
        # testing = testing,
        adapt_timestep_mode = adapt_timestep_mode,#1,
        non_dimensionalize=0,
        mode_2d = mode_2d,
        convection_Cdivu = convection_Cdivu,
        test_laplacian = test_laplacian,
        
        # ns_advection = false,

        # auto_reinit = true,
        # # ns_advection = false, #?
        # save_length = true,
    )


    ##################################################
    #Tests for operators 
    ##################################################
    if test_case == "channel_Dirichlet_pressure"
        visc_term = current_i
        
        # Test if the viscous term equals the pressure term for Poiseuille
        # skip =true : skip as long as the discretization does not handle variable spacings
        # here: v_{i+1} v_i and v_{i-1} with h and h/2 spacings
        # factor 4/3 missing, we have 3L/4 + O(h)
        @testset "laplacian" begin
        #TODO 
        @test (visc_term - (p_top-p_bottom)/L0)/((p_top-p_bottom)/L0) ≈0 atol=test_tolerance skip=true  
        end #@testset "laplacian"

        @testset "laplacian 4/3" begin
        #For now factor 4/3 missing
        @test (4*visc_term/3 - (p_top-p_bottom)/L0)/((p_top-p_bottom)/L0) ≈0 atol=test_tolerance 
        end #@testset "laplacian 4/3"

        printstyled(color=:green, @sprintf "\n visc_tem : %.2e grad p %.2e\n" visc_term (p_top-p_bottom)/L0 )
        current_i = 2

    end
    ##################################################

    printstyled(color=:green, @sprintf "\n max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(phL.u)) maximum(abs.(phL.v)))


    if test_case == "channel_Dirichlet_constant_vel" || 
       test_case == "imposed_radius"

        for iscal=1:nb_transported_scalars
            @testset "min" begin
                @test (minimum(phL.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]≈0 atol=test_tolerance
            end
   
            @testset "max" begin
                @test (maximum(phL.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]≈0 atol=test_tolerance
            end
        end
    
    # elseif test_case == "imposed_radius"

    #     @testset "min" begin
    #         @test (minimum(phL.trans_scal[:,:,iscal])-concentration0[iscal])/concentration0[iscal]≈0 atol=test_tolerance
    #     end

    elseif test_case == "channel_Dirichlet_pressure"
        @testset "Poiseuille corner velocity" begin
            @test (Poiseuille_fmax(gv.x[1,1],num.v_inlet,num.L0)-phL.v[1,1])/v_inlet≈0 atol=test_tolerance
        end    

    end





    # vec = Poiseuille_fmax.(gv.x,num.v_inlet,num.L0)
    # print("\n Poiseuille_fmax ",vec[1,:])
    # print("\n ")
    # print("\n phL.v ",phL.v[1,:])
    # print("\n ")


    # print("\n rel err",abs.(phL.v[1,:].-vec[1,:])./num.v_inlet)


    # print("\n ", gv.x[1,:])
    # print("\n ", gv.x[1,1])

    print("\n P ",Poiseuille_fmax(gv.x[1,1],num.v_inlet,num.L0)," v ",phL.v[1,1]," test ",Poiseuille_fmax(0.0,num.v_inlet,num.L0))


    #Attention = vs copy




    if test_case == "channel_Dirichlet_pressure"
        
        if write_h5

            if prediction == 0
                str_prediction = "prediction_Flower"
            elseif prediction == 4
                str_prediction = "prediction_PIII"
            else 
                str_prediction = @sprintf "_%.1i" prediction
            end

            filename="v_err"*"_"*str_prediction*"_"
            file = prefix*filename*"_"*striter*".h5"
            A .= (A .- Poiseuille_fmax.(gv.x,num.v_inlet,num.L0)) ./num.v_inlet
            # A=transpose(A)
            h5write2(file, "data", A,"w")
        
        end # write_h5

        printstyled(color=:red, @sprintf "\n iMx %.10e Mx %.10e eps %.10e \n" 1/(1e-4/128)^2 (1e-4/128)^2 eps(0.01))

        print("\n n ",n," nx ",nx," ny ",ny)

        printstyled(color=:green, @sprintf "\n p_bottom : %.10e p_top %.10e grad %.10e grad %.10e \n" p_bottom p_top -(p_top-p_bottom)/L0 8*mu1/L0^2*v_inlet)


    end #test_case == "channel_Dirichlet_pressure"






end #run_case function


#################################################
#Tests
#################################################

n = 64 
n = 128

max_iter=100

prediction = 0 #default in Flower
# prediction = 1 #pressure in prediction
# prediction = 2
# prediction = 3 #PIII


@testset verbose = true "Poiseuille imposed pressure" begin
    test_case = "channel_Dirichlet_pressure"
    # TODO not working with eps() (nearly 2e-16)
    run_case(test_case,n,1,prefix,prediction,1e-14)

    #TODO test pressure-velocity coupling methods 

end

@testset "Constant velocity with half circle" begin
    test_case = "channel_Dirichlet_constant_vel_half-circle"
    # TODO not working with eps() (nearly 2e-16)
    run_case(test_case,n,1,prefix,prediction,1e-14)
    # run_case(test_case,n,100,prefix,prediction,1e-14)
end


@testset "Imposed radius with Dirichlet" begin
    test_case = "imposed_radius_dir"
    run_case(test_case,n,100,prefix,prediction,1e-14)
end

@testset "Imposed radius with Butler-Volmer" begin
    test_case = "imposed_radius"
    run_case(test_case,n,100,prefix,prediction,1e-14)
end

@testset "Imposed radius with Butler-Volmer4" begin
    test_case = "imposed_radius4"
    run_case(test_case,n,100,prefix,prediction,1e-14)
end

@testset "Constant velocity" begin
    test_case = "channel_Dirichlet_constant_vel"
    # TODO not working with eps() (nearly 2e-16)
    run_case(test_case,n,1,prefix,prediction,1e-14)
    run_case(test_case,n,100,prefix,prediction,7e-14) #TODO error accumulates, > 1e-14 
end

#Half circle : which velocity ?
# @testset "Constant velocity with half circle" begin
#     test_case = "channel_Dirichlet_constant_vel_half-circle"
#     # TODO not working with eps() (nearly 2e-16)
#     run_case(test_case,n,1,prefix,prediction,1e-14)
#     # run_case(test_case,n,100,prefix,prediction,1e-14)
# end

# @testset "Poiseuille imposed pressure" begin
#     test_case = "channel_Dirichlet_pressure"
#     # TODO not working with eps() (nearly 2e-16)
#     run_case(test_case,n,1,prefix,prediction,1e-14)
# end


# test_Poiseuille(num,vD,grid_v)
# #TODO
# compute_grad_p!(num,grid, grid_u, grid_v, pD, opC_p, opC_u, opC_v)