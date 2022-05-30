using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

num = Numerical(case = "Sphere",
    L0 = 2.,
    n = 64,
    CFL = 1.0,
    u_inf = 0.0,
#     TEND = 0.5,
    # TEND = 0.1,
    R = 0.7,
    max_iterations = 200)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)
fwd.TL .= 0.

MIXED, SOLID, LIQUID  = run_forward(num, idx, idxu, idxv, tmp, fwd,
BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = 1.0),
    left = Boundary(t = dir, f = dirichlet, val = 1.0),
    right = Boundary(t = dir, f = dirichlet, val = 1.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 1.0),
),
stefan = true,
advection = true,
heat = true,
heat_solid_phase = false,
heat_liquid_phase = true,
navier_stokes = true,
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

fV = Figure(resolution = (1600, 1000))
colsize!(fV.layout, 1, Aspect(1, 1.0))
ax = Axis(fV[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.X[1,:], num.Y[:,1], fwd.V')
contour!(num.X[1,:], num.Y[:,1], fwd.u', levels = 0:0, color=:red, linewidth = 3);
cbar = fV[1,2] = Colorbar(fV, hmap)
# limits!(ax, 0.0, lim, -lim/2, lim/2)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fV)

fVu = Figure(resolution = (1600, 1000))
colsize!(fVu.layout, 1, Aspect(1, 1.0))
ax = Axis(fVu[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.Xu[1,:], num.Yu[:,1], fwd.Vu')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
cbar = fVu[1,2] = Colorbar(fVu, hmap)
# limits!(ax, 0.0, lim, -lim/2, lim/2)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fVu)

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.X[1,:], num.Y[:,1], fwd.TL')
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidth = 3);
cbar = fT[1,2] = Colorbar(fp, hmap)
# limits!(ax, 0.0, lim, -lim/2, lim/2)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.Xu[1,:], num.Yu[:,1], fwd.uL')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
# scatter!([num.Xu[CartesianIndex(25,12)]], [num.Yu[CartesianIndex(25,12)]], color=:red)
cbar = fu[1,2] = Colorbar(fu, hmap)
limits!(ax, -lim, lim, -lim, lim)
# limits!(ax, -1.0, -0.5, -0.5, 0.0)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.Xv[1,:], num.Yv[:,1], fwd.vL')
contour!(num.Xv[1,:], num.Yv[:,1], fwd.uv', levels = 0:0, color=:red, linewidth = 3);
cbar = fv[1,2] = Colorbar(fv, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = -4:0.5:4, yticks = -4:0.5:4)  # customized as you see fit
hmap = heatmap!(num.X[1,:], num.Y[:,1], (fwd.pL.*num.τ)')
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidth = 3);
cbar = fp[1,2] = Colorbar(fp, hmap)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

# pref = "/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/shrinking/"
# make_video(num, fwd, "u"; title_prefix=pref,
#         title_suffix="", framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0), minv = -0.5, maxv = 0.5)
# make_video(num, fwd, "v"; title_prefix=pref,
#         title_suffix="", framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.3, maxv = 0.3)
# make_video(num, fwd, "p"; title_prefix=pref,
#         title_suffix="", framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.01, maxv = 0.01)
# make_video(num, fwd, "T"; title_prefix=pref,
#         title_suffix="", framerate=20)
