using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

num = Numerical(case = "Sphere",
    L0 = 9.,
    n = 192,
    Re = 1.0,
    CFL = 1.0,
    R = 0.5,
    u_inf = 0.0,
    max_iterations = 90)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

# (pL, uL, vL, ApL, AuL, AvL, LpL, LuL, LvL,
# LCUTp, LCUTu, LCUTv,
# ksppL, kspuL, kspvL, nsL, ns_vecL,
# GxpL, GypL, DxuL, DyvL, LCUTDx, LCUTDy,
# MpL, iMpL, MuL, MvL, iMGxL, iMGyL, iMDxL, iMDyL,
# τ, iRe, n, MIXED, LIQUID) = run_forward(num, idx, idxu, idxv, tmp, fwd,
MIXED, SOLID, LIQUID, iMDx, iMDy, CUTDx, CUTDy = run_forward(num, idx, idxu, idxv, tmp, fwd,
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
# speed = -num.Δ / 32. / num.τ,
verbose = true,
show_every = 1
)

lim = 1.5
lim = num.L0 / 2
# lim = 1.0

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/moving/rigid_solid/"
suffix = "_192_1"
# suffix = ""

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
heatmap!(num.Xu[1,:], num.Yu[:,1], fwd.uL')
contour!(num.Xu[1,:], num.Yu[:,1], fwd.uu', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

make_video(num, fwd, "u"; title_prefix=prefix,
        title_suffix=suffix, framerate=20, minv=0.0, maxv=25.6, limitsx=(-lim,lim), limitsy=(-lim,lim))
 
fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
heatmap!(num.Xv[1,:], num.Yv[:,1], fwd.vL')
contour!(num.Xv[1,:], num.Yv[:,1], fwd.uv', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

# make_video(num, fwd, "v"; title_prefix=prefix,
#         title_suffix=suffix, framerate=20, minv=-0.4, maxv=0.4, limitsx=(-lim,lim), limitsy=(-lim,lim))

pavg = mean(fwd.pL[LIQUID].*num.τ)
pstd = std(fwd.pL[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
heatmap!(num.X[1,:], num.Y[:,1], (fwd.pL.*num.τ)')
contour!(num.H, num.H, fwd.u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

# make_video(num, fwd, "p"; title_prefix=prefix,
#         # minv=pavg-pstd, maxv=pavg+pstd,
#         # minv = -0.25, maxv = 0.25,
#         step0 = 2,
#         title_suffix=suffix, framerate=20, limitsx=(-lim,lim), limitsy=(-lim,lim))

divu = reshape(iMDx * (tmp.DxuL * vec(fwd.uL) .+ CUTDx) .+ iMDy * (tmp.DyvL * vec(fwd.vL) .+ CUTDy), (num.n, num.n))

fdu = Figure(resolution = (1600, 1000))
colsize!(fdu.layout, 1, Aspect(1, 1.0))
ax = Axis(fdu[1,1], aspect = 1, xticks = -4:1:4, yticks = -4:1:4)  # customized as you see fit
hmap = heatmap!(num.X[1,:], num.Y[:,1], divu')
contour!(num.X[1,:], num.Y[:,1], fwd.u', levels = 0:0, color=:red, linewidth = 3);
cbar = fdu[1,2] = Colorbar(fdu, hmap, labelpadding=0)
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fdu)