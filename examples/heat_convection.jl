using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

L0 = 6.
n = 64
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

num = Numerical(case = "Cylinder",
    Re = 10.0,    
    CFL = 0.5,
    x = x,
    y = y,
    max_iterations = 60,
    u_inf = 1.0,
    R = 1.5,
    save_every = 1,
    ϵ = 0.05)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)
phL.T .= 0.

@time MIXED, SOLID, LIQUID = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_TL = Boundaries(
        bottom = Dirichlet(val = 1.0),
        top = Dirichlet(val = 1.0)
    ),
    BC_uL = Boundaries(
        left = Dirichlet(val = num.u_inf),
        bottom = Dirichlet(val = num.u_inf),
        top = Dirichlet(val = num.u_inf)
    ),
    BC_vL = Boundaries(
        left = Dirichlet(val = num.v_inf),
        bottom = Dirichlet(val = num.v_inf),
        top = Dirichlet(val = num.v_inf)
    ),
    BC_pL = Boundaries(
        right = Dirichlet(),
    ),
    time_scheme = FE,
    stefan = true,
    advection = true,
    heat = true,
    heat_convection = true,
    heat_solid_phase = false,
    heat_liquid_phase = true,
    navier_stokes = true,
    ns_advection = true,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1
)

f = heatmap!(gp.x[1,:], gp.y[:,1], phL.T', colormap= Reverse(:ice))
f = contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
f = current_figure()

lim = num.L0 / 2

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)
contourf!(gp.x[1,:], gp.y[:,1], phL.T', colormap=:dense, colorrange=(0.2, 1.0))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 5);
contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

# make_video(num, fwd, gp, "T"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", maxv = 1.0, framerate=40)

# make_video(num, fwd, gu, "u"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", maxv = 15.0, framerate=40)

# make_video(num, fwd, gp, "p"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", framerate=40)

make_video(gu, fwd.ux, fwdL.u; title_prefix="u_field",
        title_suffix="", framerate=240)
make_video(gv, fwd.uy, fwdL.v; title_prefix="v_field",
        title_suffix="", framerate=240)
make_video(gp, fwd.u, fwdL.T; title_prefix="T_field",
        title_suffix="", framerate=240)