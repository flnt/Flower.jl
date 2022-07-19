using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

L0 = 2.
n = 64
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
    R = 0.7,
    max_iterations = 200)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
phL.T .= 0.

MIXED, SOLID, LIQUID  = run_forward(num, gp, gu, gv,
opS, opL, phS, phL, fwd,
BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = 1.0),
    left = Boundary(t = dir, f = dirichlet, val = 1.0),
    right = Boundary(t = dir, f = dirichlet, val = 1.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 1.0),
),
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

# f = heatmap!(num.H, num.H, (fwd.TL+fwd.TS)', colormap= Reverse(:ice))
# f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

# for j = num.max_iterations÷5:num.max_iterations÷5:num.max_iterations
#     f = contour!(num.H, num.H, fwd.usave[j,:,:]', levels = 0:0, color=:black, linewidth = 3);
# end

# f = current_figure()

lim = 1.0

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.T')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fT[1,2] = Colorbar(fp, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fu[1,2] = Colorbar(fu, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fv[1,2] = Colorbar(fv, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fp[1,2] = Colorbar(fp, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fphi = Figure(resolution = (1600, 1000))
colsize!(fphi.layout, 1, Aspect(1, 1.0))
ax = Axis(fphi[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.ϕ')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
cbar = fphi[1,2] = Colorbar(fphi, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fphi)

# pref = "/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/shrinking/"
# suff = ""
# make_video(num, fwd, gu, "u"; title_prefix=pref,
#         title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0), minv = -0.5, maxv = 0.5)
# make_video(num, fwd, gv, "v"; title_prefix=pref,
#         title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.3, maxv = 0.3)
# make_video(num, fwd, gp, "p"; title_prefix=pref,
#         title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0), minv = -0.003, maxv = 0.003)
# make_video(num, fwd, gp, "T"; title_prefix=pref,
#         title_suffix=suff, framerate=20)
