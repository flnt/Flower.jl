using Revise
using Flower

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

h0 = 0.8

Ny = 4
L0x = 2.0
L0y = 2.0 * Ny
n = 64

x = LinRange(-L0x/2, L0x/2, n+1)
y = LinRange(-L0y/2, L0y/2, Ny*n+1)

max_its = 17400

num = Numerical(
    case = "Drop",
    x = x,
    y = y,
    Re = 1.0,
    CFL = 0.25,
    max_iterations = max_its,
    u_inf = 0.0,
    v_inf = 0.0,
    shifted = 0.0,
    save_every = max_its÷100,
    # save_every = 1,
    reinit_every = 3,
    nb_reinit = 2,
    δreinit = 10.0,
    σ = 1.0,
    ϵ = 0.03,
    g = 10.0,
    β = 0.0,
    R = h0,
    A = 0.3,
    NB = 24,
    nLS = 2
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv);

gp.LS[2].u .= 3.7 .- gp.y
# combine_levelsets!(num, gp)

phL.u .= num.u_inf
phL.v .= num.v_inf

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = true,
    BC_uL = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Navier_cl(λ = 1e-2),
        top = Dirichlet(),
    ),
    BC_vL = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Dirichlet(),
        top = Dirichlet()
    ),
    BC_pL = Boundaries(
        left = Periodic(),
        right = Periodic(),
    ),
    BC_u = Boundaries(
        left = Periodic(),
        right = Periodic(),
        bottom = Neumann_cl(θe = π / 2.0),
        top = Neumann_inh(),
    ),
    BC_int = [FreeSurface(), Wall()],
    time_scheme = FE,
    navier_stokes = true,
    ns_advection = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    save_length = true,
)

tcks = -num.L0/2:2:num.L0
lim = (num.L0 + num.Δ) / 2
lim = 1.0

fu = Figure(size = (200, 800))
ax = Axis(fu[1,1], aspect = DataAspect(), xticks = tcks, yticks = tcks)
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.LS[end].u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fu)

fv = Figure(size = (200, 800))
ax = Axis(fv[1,1], aspect = DataAspect(), xticks = tcks, yticks = tcks)
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.LS[end].u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fv)

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/free_surface/2LS/"
suffix = "_drip_β$(num.β)_g$(num.g)_A$(num.A)_it$(num.max_iterations)"

limx = num.x[end]
limy = num.y[end]
# make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field",
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(-limy,limy))
# make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field",
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy-num.Δ/2,limy+num.Δ/2))
# make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field",
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))
# make_video(num, gp, fwd.u, fwd.κ[end,:,:,:]; title_prefix=prefix*"k_field",
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))#, minv=-4.0, maxv=4.0)
# make_video(num, gp, fwd.u title_prefix=prefix*"none",
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))#, stepf=444)

nothing