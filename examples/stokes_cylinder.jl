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

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

@time MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
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
heatmap!(num.Xu[1,:], num.Yu[:,1], fwd.uL')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(num.Xv[1,:], num.Yv[:,1], fwd.vL')
contour!(num.Xv[1,:], num.Yv[:,1], fwd.uv', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

pavg = mean(fwd.pL[LIQUID].*num.τ)
pstd = std(fwd.pL[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(num.X[1,:], num.Y[:,1], (fwd.pL.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fCd = Figure(resolution = (1600, 1000))
colsize!(fCd.layout, 1, Aspect(1, 1.0))
ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")  # customized as you see fit
lines!(fwd.Cd)
limits!(ax, 0, size(fwd.Cd, 1), 0, 10)
resize_to_layout!(fCd)
