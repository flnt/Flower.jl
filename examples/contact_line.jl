using Revise
using Flower
include("../src/run_one_phase.jl")

fontsize_theme = Theme(fontsize=30)
set_theme!(fontsize_theme)

T = Float64 # not used 

Lx = 4
Ly = 4
nx = 64
ny = nx

num = Numerical( # defined in types.jl
    # DDM parameters
    subdomains=2,
    overlaps=1,
    tolu=1e-8, # velocity
    tolp=1e-8, # pressure
    tolt=1e-8, # temperature
    # discrete ranges
    x=LinRange(-Lx / 2, Lx / 2, nx + 1),
    y=LinRange(-Ly / 2, Ly / 2, ny + 1),
    # physical parameters
    Re=1.0,
    CFL=1.0, # backwards Euler
    max_iterations=200,
    u_inf=0.0,
    v_inf=1.0,
    save_every=1, #
    σ=1.0, # surface tension
    ϵ=0.05, # 
    g=0.0, # gravity
    β=0, # angle of inclination
    # shifted=0,
    # R=1.0,
    # shift_y=1.5,
    nb_reinit = ny,
    case="Planar", # params: shifted
    shifted=-1.e-3,

    # case="Drop", # params: R, A, shifter 
    # shifted=0.0, 
    # R=h0, # radius of drop
    # A=0.02, # amplitude of perturbation

    NB=0, # number of cells that the velocity is extended from the interface
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@. gp.u = gp.x + num.shifted
tracer = copy(gp.u)
gp.u .= 1.; #level set is always equal to 1 => only LIQUID phase
initial_tracer = copy(tracer)

@time MIXED, MIXED_u, MIXED_v, SOLID, LIQUID = run_forward_one_phase(num, gp, gu, gv,
    opL, opC_pL, opC_uL, opC_vL, 
    phL, fwd, tracer;
    periodic_x = true,
    BC_uL=Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom=Boundary(t=dir, f=dirichlet, val=num.v_inf),
    ),
    BC_vL=Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom=Boundary(t=dir, f=dirichlet, val=0.0),
    ),
    BC_pL=Boundaries(
        bottom=Boundary(t=neu, f=neumann, val=0.0),
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        top=Boundary(t=dir, f=dirichlet, val=0.0),
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    advection = true, #move the level set
    ns_advection = true,
    navier_stokes = true,
    levelset = true,
    verbose = true,
    adaptative_t = false,
    show_every = 100
    )

# @show (initial_tracer == tracer)

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
#contour!(gu.x[1, :], gu.y[:, 1], fwd.uusave[end, :, :]',levels=0:0, color=:blue, linewidth=3)

contour!(gp.x[1, :], gp.y[:, 1], fwd.Tsave[2,:,:]', levels=0:0, color=:red, linewidth=2.,linestyle=:dash)
contour!(gp.x[1, :], gp.y[:, 1], tracer', levels=0:0, color=:red, linewidth=3)
Makie.save("contact_line_fu.png", fu)

fv = Figure();
ax = Axis(fv[1, 1],
    title="t=$(round(fwd.tv[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gv.x[1, :], gv.y[:, 1], fwd.Uysave[end, :, :]')
Colorbar(fv[1, 2], hm, label="v")
contour!(gp.x[1, :], gp.y[:, 1], tracer',
    levels=0:0, color=:red, linewidth=3);
Makie.save("contact_line_fv.png", fv)

fp = Figure()
ax = Axis(fp[1, 1],
    title="t=$(round(fwd.tv[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gp.x[1, :], gp.y[:, 1], fwd.psave[end, :, :]')
Colorbar(fp[1, 2], hm, label="p")
contour!(gp.x[1, :], gp.y[:, 1], tracer',
    levels=0:0, color=:red, linewidth=3);
Makie.save("contact_line_fp.png", fp)