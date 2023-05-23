using Revise, Flower
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)


dx = 0.07
Ly = 2.
ny = floor(Int, Ly / dx)

Lx = pi
nx = floor(Int, Lx / dx)

dx= Lx/nx
dy= Ly/ny
@show nx,ny,dy

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
    max_iterations=30,
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

# default 
gu.θe[:,:] .= num.θe*pi/180.
gv.θe[:,:] .= num.θe*pi/180.
gp.θe[:,:] .= num.θe*pi/180.

# initial conditions
for i in 1:gu.nx
# velocity field with Couette profile u = y
    phL.u[:,i] .= gu.y[:,1]
end

for i in 1:gu.nx
    gu.θe[end, :] .= (i > gu.nx ÷ 2) ? 76.4 * π / 180.0 : 85.3 * π / 180.0
    gu.θe[1, :] .= (i > gu.nx ÷ 2) ? 85.3 * π / 180.0 : 76.4 * π / 180.0
end
for i in 1:gp.nx
    gp.θe[end, :] .= (i > gp.nx ÷ 2) ? 76.4 * π / 180.0 : 85.3 * π / 180.0
    gp.θe[1, :] .= (i > gp.nx ÷ 2) ? 85.3 * π / 180.0 : 76.4 * π / 180.0
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
ax = Axis(fu[1, 1],title="t=$(round(fwd.t[time], digits=3))",xlabel="x",ylabel="y")
colsize!(fu.layout, 1, Aspect(1,float(Lx/Ly)))

hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwdL.u[time, :, :]')
Colorbar(fu[1, 2], hm, label="u")

for x in gp.x[1, :]
    lines!(ax, [x, x], [minimum(gp.y), maximum(gp.y)], linewidth=0.2, color=:black)
end
for y in gp.y[:, 1]
    lines!(ax, [minimum(gp.x), maximum(gp.x)], [y, y], linewidth=0.2, color=:black)
end

mid = nx ÷ 2
lines!(fwdL.u[1,:,mid], gu.y[:,mid], color=:gray)

contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[2, :, :]', levels=0:0, color=:red, linewidth=1.5, linestyle=:dotted)
contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[time, :, :]', levels=0:0, color=:red, linewidth=2.5, linestyle=:solid)

# plot Young's shear stress at the bottom (blue)
yss = gu.Young[1,:];yss = yss./maximum(yss)/3
lines!(gu.x[1,:], yss.+gu.y[1, 1], color=:blue, linewidth=4)

# plot Young's shear stress at the top (blue)
yss = gu.Young[end,:];yss = yss./maximum(yss)/3
lines!(gu.x[1,:], yss.+gu.y[1, 1], color=:yellow, linewidth=4)


fname = "contact_line_Couette.png"
@show fname
Makie.save(fname, fu)