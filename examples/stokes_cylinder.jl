using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

num = Numerical(case = "Cylinder",
    L0 = 6.,
    n = 64,
    Re = 1.0,
    CFL = 1.0,
    TEND = 0.5,
    R = 0.5,
    u_inf = 1)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
opS, opL, phS, phL, fwd,
BC_uL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    # bottom = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    # top = Boundary(t = dir, f = dirichlet, val = num.u_inf)
    ),
BC_vL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = 0.0),
    # bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
    # top = Boundary(t = dir, f = dirichlet, val = 0.0)
    ),
stefan = false,
advection = false,
heat = false,
navier_stokes = true,
ns_advection = false,
ns_solid_phase = false,
ns_liquid_phase = true,
verbose = true,
show_every = 1
)

tcks = -num.L0/2:2:num.L0
lim = num.L0 / 2

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
contour!(num.H, num.H, gp.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fCd = Figure(resolution = (1600, 1000))
colsize!(fCd.layout, 1, Aspect(1, 1.0))
ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")  # customized as you see fit
lines!(fwd.Cd)
limits!(ax, 0, size(fwd.Cd, 1), 0, 10)
resize_to_layout!(fCd)
