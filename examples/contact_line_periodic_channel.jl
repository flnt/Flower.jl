using Revise, Flower
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

Ly = 1
ny = 32

Lx = 1
nx = 32

dx= Lx/nx
dy= Ly/ny

num = Numerical( # defined in types.jl
    # DDM parameters
    subdomains=2,
    overlaps=1,
    tolu=1e-8, # velocity
    tolp=1e-8, # pressure
    tolt=1e-8, # temperature
    x=LinRange(-Lx / 2, Lx / 2, nx + 1),
    y=LinRange(-Ly / 2, Ly / 2, ny + 1),
    # physical parameters
    Re=100.0,
    CFL=0.5, # backwards Euler
    max_iterations=2,
    u_inf=0.0,
    v_inf=1.0,
    save_every=1, #
    σ=1.0, # surface tension
    ϵ=0.05, # 
    g=0.0, # gravity
    β=0, # angle of inclination
    nb_reinit=ny,

    case="Nothing",
    shifted=-1.e-6,
    NB=4, # number of cells that the velocity is extended from the interface
    # Contact angle parameters
    Ca=0.1, # Capillary number
    λCA=3*dy, # slip lenght > dy 
    εCA=3*dx, # width
    θe=70.0, # prescribed contact angle
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdL, fwdS = init_fields(num, gp, gu, gv)

# initial condition of the tracer 
tracer = gp.x .+ num.shifted

# initial condition of the velocity field -> Poiseuille profile u = 1 - y^2
for i in 1:gu.ny
    phL.u[i,:] .= 1.0 .- gu.y[1,:].^2
end

@time run_forward_one_phase(num, gp, gu, gv, opL, opC_pL, opC_uL, opC_vL, phL, fwd, fwdL, tracer;
    periodic_x=true,
    BC_uL=Boundaries(
        left  =Boundary(t=per, f=periodic),
        right =Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=0.0),
        top   =Boundary(t=nav, f=navier, val=0.0),
        ),
    BC_vL=Boundaries(
        left=  Boundary(t=per, f=periodic),
        right =Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=0.0),
        top   =Boundary(t=nav, f=navier, val=0.0),
    ),
    BC_pL=Boundaries( # pressure 
        left  =Boundary(t=per, f=periodic),
        right =Boundary(t=per, f=periodic),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
        top   =Boundary(t=neu, f=neumann, val=0.0),
    ),
    BC_u=Boundaries( # levelset/tracer 
        left  =Boundary(t=per, f=periodic),
        right =Boundary(t=per, f=periodic),
        # bottom=Boundary(t=neu, f=neumann, val=0.0),
        # top   =Boundary(t=neu, f=neumann, val=0.0),
    ),
    advection=true,
    ns_advection=true,
    navier_stokes=true,
    levelset=true,
    verbose=true,
    adaptative_t=false,
    show_every=100
    )

time = num.max_iterations # time to plot figure 

yres = 500 # figure vertical resolutin
xres=yres*float(Lx/Ly)
fu = Figure(resolution = (xres, yres))

ax = Axis(fu[1, 1],title="θe=$(num.θe), t=$(round(fwd.t[time], digits=3))",xlabel="x",ylabel="y")
colsize!(fu.layout, 1, Aspect(1,float(Lx/Ly)))

hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwdL.u[time, :, :]')
Colorbar(fu[1, 2], hm, label="u")

# contour!(gu.x[1, :], gu.y[:, 1], fwdL.u[time, :, :]', levels=0:0, color=:black, linewidth=2, linestyle=:solid)
# lines!(ax, [x, x], [minimum(gp.y), maximum(gp.y)], linewidth=0.5, color=:black)

# mid = nx/2
# lines!(gu.x[mid,:], fwdL.u[1, mid, 2])

# contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[2, :, :]', levels=0:0, color=:red, linewidth=2, linestyle=:dash)
# contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[time, :, :]', levels=0:0, color=:red, linewidth=2, linestyle=:solid)

#resize_to_layout!(fu)

fname = "contact_line_periodic_channel_fu_theta_e_$(num.θe).png"
@show fname
Makie.save(fname, fu)