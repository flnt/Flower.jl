using Revise
using Flower

fontsize_theme = Theme(fontsize=30)
set_theme!(fontsize_theme)

T = Float64 # not used 

Lx = 1.5
Ly = 4.0
nx = 32
ny = 32

num = Numerical( # defined in types.jl
    # DDM parameters
    subdomains=2,
    overlaps=1,
    tolu=1e-8, # velocity
    tolp=1e-8, # pressure
    tolt=1e-8, # temperature
    # discrete ranges
    x=LinRange(0, Lx, nx + 1),
    y=LinRange(-Ly / 2, Ly / 2, ny + 1),
    # physical parameters
    Re=1.0,
    CFL=0.1, # backwards Euler
    max_iterations=1000,
    u_inf=0.0,
    v_inf=1.0,
    save_every=1, #
    σ=1.0, # surface tension
    ϵ=0.05, # 
    g=0.0, # gravity
    β=0, # angle of inclination
    # case="Cylinder", # params: shifted
    # shifted=0,
    # R=1.0,
    # shift_y=1.5,

    case="Planar", # params: shifted
    shifted=-1.e-3,

    # case="Drop", # params: R, A, shifter 
    # shifted=0.0, 
    # R=h0, # radius of drop
    # A=0.02, # amplitude of perturbation

    NB=nx, # number of cells that the velocity is extended from the interface
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd = init_fields(num, gp, gu, gv)

# phL.u .= num.u_inf
# phL.v .= num.v_inf

# set inital pDirichlet field so that there's a pressure gradient at the bottom
# @inbounds @threads for II in gp.ind.b_bottom[1]
#     pII = lexicographic(II, gp.ny)
#     phL.pD.data[2][pII] = - num.g*cos(num.β) * gp.dy[1,1] / 2
# end

@time MIXED, MIXED_u, MIXED_v, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd,
    # periodic_x=false, # default
    # periodic_y=false, # default

    # BCs 
    # neumann is default 

    # u velocity (Liquid)
    BC_uL=Boundaries(
        left=Boundary(t=dir, f=dirichlet, val=0.0),
        right=Boundary(t=neu, f=neumann, val=0.0),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
        top=Boundary(t=neu, f=neumann, val=0.0),
    ),
    # v velocity (Liquid)
    BC_vL=Boundaries(
        left=Boundary(t=dir, f=dirichlet, val=num.v_inf),
        right=Boundary(t=neu, f=neumann, val=0.0),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
        top=Boundary(t=neu, f=neumann, val=0.0),
    ),
    # # u velocity (Solid)
    # BC_uS=Boundaries(
    #     left=Boundary(t=dir, f=dirichlet, val=0.0),
    #     right=Boundary(t=neu, f=neumann, val=0.0),
    #     bottom=Boundary(t=neu, f=neumann, val=0.0),
    #     top=Boundary(t=neu, f=neumann, val=0.0),
    # ),
    # # v velocity (Solid)
    # BC_vS=Boundaries(
    #     left=Boundary(t=dir, f=dirichlet, val=num.v_inf),
    #     right=Boundary(t=neu, f=neumann, val=0.0),
    #     bottom=Boundary(t=neu, f=neumann, val=0.0),
    #     top=Boundary(t=neu, f=neumann, val=0.0),
    # ),
    # pressure (Liquid)
    # dirichlet wherever you put Neumann in the velocity and viceversa
    BC_pL=Boundaries(
        left=Boundary(t=neu, f=neumann, val=0.0),
        right=Boundary(t=dir, f=dirichlet, val=0.0),
        bottom=Boundary(t=dir, f=dirichlet, val=0.0),
        top=Boundary(t=dir, f=dirichlet, val=0.0),
    ),
    BC_pS=Boundaries(
        left=Boundary(t=neu, f=neumann, val=0.0),
        right=Boundary(t=dir, f=dirichlet, val=0.0),
        bottom=Boundary(t=dir, f=dirichlet, val=0.0),
        top=Boundary(t=dir, f=dirichlet, val=0.0),
    ),
    # levelset field -> start with Neumann
    BC_u=Boundaries(
        left=Boundary(t=neu, f=neumann, val=0.0),
        right=Boundary(t=neu, f=neumann, val=0.0),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
        top=Boundary(t=neu, f=neumann, val=0.0),
    ),
    # defined in run.jl 
    stefan=false,
    advection=true,
    heat=false,
    heat_convection=false,
    ns_advection=false,
    navier_stokes=true,
    heat_liquid_phase=false,
    heat_solid_phase=false,
    ns_liquid_phase=true,
    ns_solid_phase=false,
    free_surface=true,
    hill=false,
    Vmean=false,
    levelset=true,
    one_phase=true,
    speed=0,
    analytical=false,
    verbose=true,
    show_every=10,
    save_length=true,
    save_radius=false,
    adaptative_t=true,
    Ra=0,
    λ=1
)

# tcks = -num.L0/2:2:num.L0
# lim = (num.L0 + num.Δ) / 2
# lim = 1.0

# 3 different grids: 
# gp - we use for temperature and pressure, 
# gu - we use for the velocity u
# gv - we use for the velocity v
# Each grid has a levelset (gp.u, gu.u and gv.u) 
#  that is used to compute the corresponding capacities
# The actual levelset is gp.u, which is the one that is advected
# gu.u and gv.u are just interpolated fields from gp.u

fu = Figure();
ax = Axis(fu[1, 1],
    title="t=$(round(fwd.tv[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwd.Uxsave[end, :, :]')
Colorbar(fu[1, 2], hm, label="u")
contour!(gu.x[1, :], gu.y[:, 1], fwd.uusave[end, :, :]', # interpolated
    levels=0:0, color=:blue, linewidth=3)
contour!(gp.x[1, :], gp.y[:, 1], fwd.usave[end, :, :]', # real one
    levels=0:0, color=:red, linewidth=3)
# Makie.save("inclined_plane_fu.png", fu)

fv = Figure();
ax = Axis(fv[1, 1],
    title="t=$(round(fwd.tv[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gv.x[1, :], gv.y[:, 1], fwd.Uysave[end, :, :]')
Colorbar(fv[1, 2], hm, label="v")
contour!(gp.x[1, :], gp.y[:, 1], gp.u',
    levels=0:0, color=:red, linewidth=3);
# Makie.save("inclined_plane_fv.png", fv)

fp = Figure()
ax = Axis(fp[1, 1],
    title="t=$(round(fwd.tv[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gp.x[1, :], gp.y[:, 1], fwd.psave[end, :, :]')
Colorbar(fp[1, 2], hm, label="p")
contour!(gp.x[1, :], gp.y[:, 1], gp.u',
    levels=0:0, color=:red, linewidth=3);
# Makie.save("inclined_plane_fp.png", fp)

# fk = Figure(resolution=(1600, 1000))
# colsize!(fk.layout, 1, Aspect(1, 1.0))
# ax = Axis(fk[1, 1], aspect=1, xticks=tcks, yticks=tcks)  # customized as you see fit
# hmap = heatmap!(gp.x[1, :], gp.y[:, 1], gp.κ')#, colorrange=(-1.5, 1.5))
# contour!(gp.x[1, :], gp.y[:, 1], gp.u', levels=0:0, color=:red, linewidth=3);
# cbar = fk[1, 2] = Colorbar(fk, hmap, labelpadding=0)
# # limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fk)

# prefix = ""
# suffix = "_dripping_β$(num.β)_g$(num.g)_A$(num.A)_it$(num.max_iterations)_ext_2"

# Makie.save(suffix, fk)
# limx = lim
# limy = lim * Ny
# make_video(num, fwd, gu, "u"; title_prefix=prefix,
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(-limy,limy))
# make_video(num, fwd, gv, "v"; title_prefix=prefix,
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy-num.Δ/2,limy+num.Δ/2))
# make_video(num, fwd, gp, "p"; title_prefix=prefix,
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))
# make_video(num, fwd, gp, "κ"; title_prefix=prefix,
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))#, minv=-4.0, maxv=4.0)
# make_video(num, fwd, gp, "none"; title_prefix=prefix,
#         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-limx,limx), limitsy=(-limy,limy))

# fp = Figure()
# ax = Axis(fp[1,1])
# heatmap!(fwd.psave[417,:,:]')
# scatter!(3, 53, color=:red)
# resize_to_layout!(fp)