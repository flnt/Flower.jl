using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

num = Numerical(case = "Cylinder",
    L0 = 4.0,
    n = 64,
    CFL = 0.5,
    max_iterations = 250,
    u_inf = 10.0,
    R = 1.5,
    save_every = 10)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
phL.T .= 0.

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    BC_TL = Boundaries(
        bottom = Boundary(t = dir, f = dirichlet, val = 1.0),
        top = Boundary(t = dir, f = dirichlet, val = 1.0)
    ),
    BC_uL = Boundaries(
        left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
        bottom = Boundary(t = dir, f = dirichlet, val = num.u_inf),
        top = Boundary(t = dir, f = dirichlet, val = num.u_inf)
    ),
    BC_vL = Boundaries(
        left = Boundary(t = dir, f = dirichlet, val = num.v_inf),
        bottom = Boundary(t = dir, f = dirichlet, val = num.v_inf),
        top = Boundary(t = dir, f = dirichlet, val = num.v_inf)
    ),
    stefan = true,
    advection = true,
    heat = true,
    heat_convection = true,
    heat_solid_phase = false,
    heat_liquid_phase = true,
    navier_stokes = true,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1
)

f = heatmap!(num.H, num.H, phL.T', colormap= Reverse(:ice))
f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
f = current_figure()

lim = num.L0 / 2

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
contourf!(gp.x[1,:], gp.y[:,1], phL.T', colormap=:dense, colorrange=(0.2, 1.0))
contour!(num.H, num.H, gp.u', levels = 0:0, color=:red, linewidth = 5);
contour!(gp.x[1,:], gp.y[:,1], fwd.usave[1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

pavg = mean(phL.p[LIQUID].*num.τ)
pstd = std(phL.p[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(num.H, num.H, gp.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

# make_video(num, fwd, gp, "T"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", maxv = 1.0, framerate=40)

# make_video(num, fwd, gu, "u"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", maxv = 15.0, framerate=40)

# make_video(num, fwd, gp, "p"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
#         title_suffix="_2", framerate=40)