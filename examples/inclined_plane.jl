using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

h0 = 0.5

Ny = 2
L0x = 2.0
L0y = 2.0 * Ny
n = 64

x = LinRange(-L0x/2, L0x/2, n+1)
y = LinRange(-L0y/2, L0y/2, Ny*n+1)

# d1 = x[2] - x[1]
# d0 = d1 / 2
# s = stretching(n, d0, d1, 0.1, 16, 16, 0.08)
# x = vcat(-reverse(s),s[2:end])

num = Numerical(
    case = "Drop",
    # case = "Cylinder",
    x = x,
    y = y,
    Re = 1.0,
    CFL = 0.7,
    max_iterations = 200,
    u_inf = 0.0,
    v_inf = 0.0,
    shifted = 0.0,
    save_every = 20,
    σ = 1.0,
    ϵ = 0.05,
    g = 10.0,
    # β = 3π/8,
    # β = π/4,
    β = 0,
    R = h0,
    A = 0.3,
    NB = 12,
    subdomains = 2,
    overlaps = 1
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

phL.u .= num.u_inf
phL.v .= num.v_inf

# set inital pDirichlet field so that there's a pressure gradient at the bottom
# @inbounds @threads for II in gp.ind.b_bottom[1]
#     pII = lexicographic(II, gp.ny)
#     phL.pD.data[2][pII] = - num.g*cos(num.β) * gp.dy[1,1] / 2
# end

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
# @profview MIXED, MIXED_u, MIXED_v, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd, fwdS, fwdL,
    periodic_x = true,
    BC_uL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
        top = Boundary(t = dir, f = dirichlet, val = num.u_inf),
    ),
    BC_vL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
        top = Boundary(t = dir, f = dirichlet, val = num.v_inf)
    ),
    BC_pL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = neu, f = neumann, val = 0.0)#num.g*cos(num.β)),
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    stefan = false,
    advection = true,
    heat = false,
    heat_convection = false,
    navier_stokes = true,
    ns_advection = false,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    free_surface = true,
    verbose = true,
    show_every = 1,
    save_length = true
)

tcks = -num.L0/2:2:num.L0
lim = (num.L0 + num.Δ) / 2
lim = 1.0

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.p')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fk = Figure(resolution = (1600, 1000))
colsize!(fk.layout, 1, Aspect(1, 1.0))
ax = Axis(fk[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
hmap = heatmap!(gp.x[1,:], gp.y[:,1], gp.κ')#, colorrange=(-1.5, 1.5))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
cbar = fk[1,2] = Colorbar(fk, hmap, labelpadding=0)
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fk)

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/free_surface/"
suffix = "_dripping_β$(num.β)_g$(num.g)_A$(num.A)_it$(num.max_iterations)_ext_st"

limx = num.x[end]
limy = num.y[end]
make_video(gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field",
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(-limy,limy))
make_video(gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field",
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy-num.Δ/2,limy+num.Δ/2))
make_video(gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field",
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))
make_video(gp, fwd.u, fwd.κ; title_prefix=prefix*"k_field",
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))#, minv=-4.0, maxv=4.0)
make_video(gp, fwd.u; title_prefix=prefix*"none",
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))
