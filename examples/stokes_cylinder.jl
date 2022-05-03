using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

num = Numerical(case = "Sphere",
    L0 = 6.,
    n = 64,
    Re = 1.0,
    CFL = 1.0,
    TEND = 1.0,
    R = 0.5)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
BC_uL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = 1.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 1.0),
    top = Boundary(t = dir, f = dirichlet, val = 1.0)),
BC_vL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = 0.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
    top = Boundary(t = dir, f = dirichlet, val = 0.0)),
stefan = false,
advection = false,
heat = false,
navier_stokes = true,
solid_phase = false,
liquid_phase = true,
verbose = true,
show_every = 5
)

lim = num.L0 / 2

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
resize_to_layout!(fu)
heatmap!(num.Xu[1,:], num.Yu[:,1], fwd.uL')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
resize_to_layout!(fv)
heatmap!(num.Xv[1,:], num.Yv[:,1], fwd.vL')
contour!(num.Xv[1,:], num.Yv[:,1], fwd.uv', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)

pavg = mean(fwd.pL[LIQUID].*num.τ)
pstd = std(fwd.pL[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
resize_to_layout!(fp)
heatmap!(num.X[1,:], num.Y[:,1], (fwd.pL.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
