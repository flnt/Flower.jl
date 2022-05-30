using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

num = Numerical(case = "Crystal",
    L0 = 4.0,
    n = 64,
    CFL = 1.0,
    TEND = 10.0,
    u_inf = 30.0,
    R = 0.5,
    A = -0.4,
    N = 6,
    save_every = 10)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)
fwd.TL .= 0.

MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = 1.0)),
BC_uL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    bottom = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    top = Boundary(t = dir, f = dirichlet, val = num.u_inf)),
BC_vL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = num.v_inf),
    bottom = Boundary(t = dir, f = dirichlet, val = num.v_inf),
    top = Boundary(t = dir, f = dirichlet, val = num.v_inf)),
stefan = false,
advection = false,
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

f = heatmap!(num.H, num.H, (fwd.TL)', colormap= Reverse(:ice))
f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
f = current_figure()

lim = num.L0 / 2

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
contourf!(num.X[1,:], num.Y[:,1], fwd.TL', colormap=:dense, colorrange=(0.2, 1.0))
contour!(num.X[1,:], num.Y[:,1], fwd.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
contourf!(num.Xu[1,:], num.Yu[:,1], fwd.uL')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
contourf!(num.Xv[1,:], num.Yv[:,1], fwd.vL')
contour!(num.Xv[1,:], num.Yv[:,1], fwd.uv', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

pavg = mean(fwd.pL[LIQUID].*num.τ)
pstd = std(fwd.pL[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
heatmap!(num.X[1,:], num.Y[:,1], (fwd.pL.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

make_video(num, fwd, "T"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/convection/",
        maxv = 1.0, framerate=40)