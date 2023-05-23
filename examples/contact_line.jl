using Revise, Flower
#set_theme!(fontsize_theme)

Lx = 1
Ly = 1
nx = 32
ny = 32
dx= Lx/nx
dy= Ly/ny
@show dx, dy

num = Numerical( # defined in types.jl
    # DDM parameters
    subdomains=2,
    overlaps=1,
    tolu=1e-8, # velocity
    tolp=1e-8, # pressure
    tolt=1e-8, # temperature
    # discrete ranges
    x=LinRange(-Lx / 2, Lx / 2, nx + 1),
    y=LinRange(-Ly / 2, Ly / 2, ny + 1),
    # physical parameters
    Re=20.0,
    CFL=0.5, # backwards Euler
    max_iterations=200,
    save_every=1, #
    σ=1.0, # surface tension
    ϵ=0.05, # 
    g=0.0, # gravity
    β=0, # angle of inclination
    nb_reinit=ny,
    case="Square",
    shifted=-1.e-6,
    R = float(Lx/20),

    # case="Drop", # params: R, A, shifter 
    # shifted=0.0, 
    # R=h0, # radius of drop
    # A=0.02, # amplitude of perturbation

    NB=4, # number of cells that the velocity is extended from the interface

    # Contact angle parameters
    Ca=0.12, # Capillary number
    λCA=3*dy, # slip lenght > dy 
    εCA=3*dx, # width
    θe=70.0, # prescribed contact angle
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdL, fwdS = init_fields(num, gp, gu, gv)

# Set initial condition to the level set 
# Done in init.jl -> cases 

# Set initial condition to the tracer 
tracer = gp.x .+ num.shifted .+ float(Lx/4)

@time run_forward_one_phase(num, gp, gu, gv, opL, opC_pL, opC_uL, opC_vL,
    phL, fwd, fwdL, tracer;
    periodic_x=true,
    BC_uL=Boundaries(
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=1.0),
    ),
    BC_vL=Boundaries(
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=0.0),
    ),
    BC_pL=Boundaries( # pressure 
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        top=Boundary(t=dir, f=dirichlet, val=0.0),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
    ),
    BC_u=Boundaries( # levelset/tracer 
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
    ),
    advection=true, #move the level set
    ns_advection=true,
    navier_stokes=true,
    levelset=true,
    verbose=true,
    adaptative_t=false,
    show_every=100
    )


time = num.max_iterations
fu = Figure(resolution = (700, 500));
ax = Axis(fu[1, 1],title="θe=$(num.θe), t=$(round(fwd.t[time], digits=3))",xlabel="x",ylabel="y")
colsize!(fu.layout, 1, Aspect(1, 1.0))
hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwdL.u[time, :, :]')
max_vel = maximum(fwdL.u[time, :, :])
min_vel = minimum(fwdL.u[time, :, :])
@show min_vel, max_vel
Colorbar(fu[1, 2], hm, label="u", ticks=[round(min_vel,digits=5),round(max_vel,digits=5)])
for x in gp.x[1, :]
    lines!(ax, [x, x], [minimum(gp.y), maximum(gp.y)], linewidth=0.2, color=:black)
end
for y in gp.y[:, 1]
    lines!(ax, [minimum(gp.x), maximum(gp.x)], [y, y], linewidth=0.2, color=:black)
end
# Plot object (i.e. the levelset gp.u)
#contour!(gp.x[1, :], gp.y[:, 1], fwdL.u[end, :, :]', levels = 0:0, color=:gray, linewidth=7)
contour!(gp.x[1, :], gp.y[:, 1], gp.u', levels = 0:0, color=:gray, linewidth=1)
# Plot the initial and current tacer 
contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[2, :, :]', levels=0:0, color=:red, linewidth=2.0, linestyle=:dash)
contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[time, :, :]', levels=0:0, color=:red, linewidth=3)

yss = gu.Young[1,:]; yss = yss./maximum(yss)/3
lines!(gu.x[1,:], yss.+gu.y[1, 1], color=:green, linewidth=3)

#fname = "contact_line_fu_theta_e_$(num.θe).png"
fname = "contact_line.png"
@show fname
Makie.save(fname, fu)

# function gif(fwdL, gp, gu, num)
#     function create_frame(fwdL, gp, gu, θe, time)
#         fu = Figure(resolution = (500, 500))
#         ax = Axis(fu[1, 1], title="θe=$θe, t=$(round(time, digits=3))", xlabel="x", ylabel="y")
#         colsize!(fu.layout, 1, Aspect(1, 1.0))
#         hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwdL.u[time, :, :]')
#         Colorbar(fu[1, 2], hm, label="u")
#         contour!(gp.x[1, :], gp.y[:, 1], gp.u', levels = 0:0.001:0, color=:gray, linewidth=7)
#         contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[2, :, :]', levels=0:0.001:0, color=:red, linewidth=2.0, linestyle=:dash)
#         contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[time, :, :]', levels=0:0.001:0, color=:red, linewidth=3)
#         yss = gu.Young[1,:]
#         yss = yss ./ maximum(yss) / 3
#         lines!(gu.x[1,:], yss .+ gu.y[1, 1], color=:green, linewidth=3)
#         return fu
#     end
    
#     frames = [create_frame(fwdL, gp, gu, num.θe, time) for time in 1:num.max_iterations]
    
#     record("contact_line_fu_theta_e_$(num.θe).mp4", frames, framerate = 15)
# end

# prefix = "/Users/r/Desktop/Codes/Flower.jl/examples/"
# suffix = "_Re$(abs(Re))_λ$(λ)_θe$(num.θe)"
# make_video(gu, fwd, u; title_prefix=prefix,title_suffix=suffix, xlabel="x", ylabel="y", colormap=:viridis,minv=0.0, maxv=0.0, limitsx=(-Lx/2,Lx/2), limitsy=(-Ly/2,Ly/2), framerate=24, step=1, step0=1, stepf=size(field_u,1))