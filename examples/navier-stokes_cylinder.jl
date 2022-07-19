using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

L0 = 16.
n = 128
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

# s1 = stretching(n÷2, 0.025, 0.2, 0.6)
# s2 = stretching(n+n÷2, 0.025, 0.1, 0.6)
# x = vcat(-reverse(s1),s2[2:end])
# y = vcat(-reverse(s1),s1[2:end])

num = Numerical(case = "Cylinder",
    Re = 100.0,
    CFL = 0.5,
    TEND = 50.0,
    x = x,
    y = y,
    R = 0.5,
    u_inf = 1.0,
    save_every = 10,
    shifted = 2.0)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

# data = load_field("/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/data_Re100.0_L016.0_n128.jld")
# fwd.uL .= data["u"]
# fwd.vL .= data["v"]
# fwd.pL .= data["p"]
# fwd.Gxm1L .= data["gx"]
# fwd.Gym1L .= data["gy"]

# @profview MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
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
    ns_advection = true,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1
)

tcks = -L0/2:2:L0
lim = L0 / 2

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

pavg = mean(phL.p[LIQUID].*num.τ)
pstd = std(phL.p[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/"
# suffix = "_Re$(num.Re)_L0$(num.L0)_n$(num.n)_p.jld"

# save_field(prefix*"data"*suffix, num, fwd)

# suffix = "_Re$(num.Re)_L0$(num.L0)_n$(num.n)_p"
# make_video(num, fwd, "u"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, minv=0.0, maxv=1.2, limitsx=(-lim,lim), limitsy=(-lim,lim))

# make_video(num, fwd, "v"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, minv=-0.4, maxv=0.4, limitsx=(-lim,lim), limitsy=(-lim,lim))

# using Peaks
# pks, vals = findmaxima(fwd.Uxsave[:,num.n÷2+5,num.n-25])
# f = 1 / ((pks[end-1]-pks[end-2])*num.τ*num.save_every)

fCd = Figure(resolution = (1600, 1000))
colsize!(fCd.layout, 1, Aspect(1, 1.0))
ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")
lines!(fwd.Cd)
limits!(ax, 0, size(fwd.Cd, 1), 0, 2)
resize_to_layout!(fCd)

fCl = Figure(resolution = (1600, 1000))
colsize!(fCl.layout, 1, Aspect(1, 1.0))
ax = Axis(fCl[1,1], xlabel="it", ylabel="Cl")
lines!(fwd.Cl)
limits!(ax, 0, size(fwd.Cl, 1), -1, 1)
resize_to_layout!(fCl)
