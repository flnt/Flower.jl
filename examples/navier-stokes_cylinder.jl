using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

num = Numerical(case = "Cylinder",
    L0 = 20.,
    n = 192,
    Re = 100.0,
    CFL = 0.5,
    TEND = 350.0,
    # TEND = 200.0,
    R = 0.5,
    u_inf = 1.0,
    save_every = 10,
    shifted = 5.0)
    # shifted = 3.0)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

# data = load_field("/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/data_Re100.0_L016.0_n128.jld")
# fwd.uL .= data["u"]
# fwd.vL .= data["v"]
# fwd.pL .= data["p"]

# @profview MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
@time MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
BC_uL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    bottom = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    top = Boundary(t = dir, f = dirichlet, val = num.u_inf)),
BC_vL = Boundaries(
    left = Boundary(t = dir, f = dirichlet, val = 0.0),
    bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
    top = Boundary(t = dir, f = dirichlet, val = 0.0)),
stefan = false,
advection = false,
heat = false,
navier_stokes = true,
ns_advection = true,
ns_solid_phase = false,
ns_liquid_phase = true,
verbose = true,
show_every = 1
)

tcks = -num.L0/2:2:num.L0
lim = num.L0 / 2
# lim = 1.0

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

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/"
suffix = "_Re$(num.Re)_L0$(num.L0)_n$(num.n).jld"

save_field(prefix*"data"*suffix, num, fwd)

suffix = "_Re$(num.Re)_L0$(num.L0)_n$(num.n)"
make_video(num, fwd, "u"; title_prefix=prefix,
        title_suffix=suffix, framerate=100, minv=0.0, maxv=1.2, limitsx=(-lim,lim), limitsy=(-lim,lim))

make_video(num, fwd, "v"; title_prefix=prefix,
        title_suffix=suffix, framerate=100, minv=-0.4, maxv=0.4, limitsx=(-lim,lim), limitsy=(-lim,lim))

using Peaks
pks, vals = findmaxima(fwd.Uxsave[:,num.n÷2+10,num.n-50])
f = 1 / ((pks[end-1]-pks[end-2])*num.τ)
