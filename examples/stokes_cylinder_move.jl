using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

L0 = 4.
n = 128
Δ = L0/(n-1)
x = [-L0 / 2 - Δ / 2 + i * Δ for i = 0:n]
y = [-L0 / 2 - Δ / 2 + i * Δ for i = 0:n]
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

num = Numerical(case = "Sphere",
    x = x,
    y = y,
    CFL = 0.5,
    u_inf = 0.0,
    R = 0.5,
    max_iterations = 30,
    ϵ = 0.05)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
phL.T .= 0.

@time MIXED, SOLID, LIQUID  = run_forward(num, gp, gu, gv,
opS, opL, phS, phL, fwd,
BC_pL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = 0.0),
    left = Boundary(t = dir, f = dirichlet, val = 0.0),
    right = Boundary(t = dir, f = dirichlet, val = 0.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
),
# BC_uL = Boundaries(
#     left = Boundary(t = dir, f = dirichlet, val = 1.0),
#     bottom = Boundary(t = dir, f = dirichlet, val = 1.0),
#     top = Boundary(t = dir, f = dirichlet, val = 1.0)),
# BC_vL = Boundaries(
#     left = Boundary(t = dir, f = dirichlet, val = 0.0),
#     bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
#     top = Boundary(t = dir, f = dirichlet, val = 0.0)),
stefan = false,
advection = true,
heat = false,
navier_stokes = true,
ns_solid_phase = false,
ns_liquid_phase = true,
speed = -num.Δ / 1. / num.τ,
verbose = true,
show_every = 1
)

lim = 1.5
lim = num.L0 / 2
# lim = 1.0

pref = "/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/moving/rigid_solid/"
suff = ""
make_video(num, fwd, gu, "u"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -0.5, maxv = 0.5)
make_video(num, fwd, gu, "ucorr"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -0.5, maxv = 0.5)
make_video(num, fwd, gv, "v"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -0.3, maxv = 0.3)
make_video(num, fwd, gp, "p"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -0.003, maxv = 0.003)
make_video(num, fwd, gp, "ϕ"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -0.003, maxv = 0.003)
# make_video(num, fwd, gp, "T"; title_prefix=pref,
#         title_suffix=suff, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))

lim = 1.0

fgx = Figure(resolution = (1600, 1000))
colsize!(fgx.layout, 1, Aspect(1, 1.0))
ax = Axis(fgx[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], reshape(opL.CUTGxp, (gu.ny, gu.nx))')
hmap = heatmap!(gu.x[1,:], gu.y[:,1], reshape(phL.Gxm1, (gu.ny, gu.nx))')
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], reshape(opL.CUTu, (gu.ny, gu.nx))')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fgx[1,2] = Colorbar(fgx, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fgx)

ftmp = Figure(resolution = (1600, 1000))
colsize!(ftmp.layout, 1, Aspect(1, 1.0))
ax = Axis(ftmp[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.tmp')
# hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.tmp')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = ftmp[1,2] = Colorbar(fgx, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(ftmp)
