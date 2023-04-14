using Revise
using Flower

fontsize_theme = Theme(fontsize=30)
set_theme!(fontsize_theme)

T = Float64 # not used 

Lx = 1
Ly = 1
nx = 32
ny = 16
dx= Lx/nx
dy= Ly/ny

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
    max_iterations=400,
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
    nb_reinit=ny,
    case="Planar", # params: shifted
    shifted=-1.e-3,

    # case="Drop", # params: R, A, shifter 
    # shifted=0.0, 
    # R=h0, # radius of drop
    # A=0.02, # amplitude of perturbation

    NB=0, # number of cells that the velocity is extended from the interface

    # Contact angle parameters
    Ca=0.1, # Capillary number
    εCA=0.03125, # width
    λCA=0.03125, # slip lenght
    θe=70.0, # prescribed contact angle
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdL, fwdS = init_fields(num, gp, gu, gv)

@. gp.u = gp.x + num.shifted
tracer = copy(gp.u)
gp.u .= 1.0; #level set is always equal to 1 => only LIQUID phase
initial_tracer = copy(tracer)

@time run_forward_one_phase(num, gp, gu, gv, opL, opC_pL, opC_uL, opC_vL,
    phL, fwd, fwdL, tracer;
    periodic_x=true,
    BC_uL=Boundaries(
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=1.0),
    ),
    BC_vL=Boundaries(
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        bottom=Boundary(t=nav, f=navier, val=0.0),
    ),
    BC_pL=Boundaries( # pressure 
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
        top=Boundary(t=dir, f=dirichlet, val=0.0),
        bottom=Boundary(t=neu, f=neumann, val=0.0),
    ),
    BC_u=Boundaries( # levelset/tracer 
        left=Boundary(t=per, f=periodic),
        right=Boundary(t=per, f=periodic),
    ),
    advection=true, #move the level set
    ns_advection=true,
    navier_stokes=true,
    levelset=true,
    verbose=true,
    adaptative_t=false,
    show_every=100
    )

    fu = Figure();
ax = Axis(fu[1, 1],
    title="t=$(round(fwd.t[end], digits=3))",
    xlabel="x",
    ylabel="y")
hm = heatmap!(gu.x[1, :], gu.y[:, 1], fwdL.u[end, :, :]')
Colorbar(fu[1, 2], hm, label="u")
#contour!(gu.x[1, :], gu.y[:, 1], fwd.uusave[end, :, :]',levels=0:0, color=:blue, linewidth=3)

contour!(gp.x[1, :], gp.y[:, 1], fwdL.T[2, :, :]', levels=0:0, color=:red, linewidth=2.0, linestyle=:dash)
contour!(gp.x[1, :], gp.y[:, 1], tracer', levels=0:0, color=:red, linewidth=3)

# for x in gp.x[1, :]
#     lines!(ax, [x, x], [minimum(gp.y), maximum(gp.y)], linewidth=1, color=:red)
# end
# for y in gp.y[:, 1]
#     lines!(ax, [minimum(gp.x), maximum(gp.x)], [y, y], linewidth=1, color=:red)
# end

# youngstress= zeros(size(gp.ind.b_bottom[1]))
# youngstress = compute_young_stress(gp,num,gp.ind.b_bottom[1])
# @show youngstress
# lines!(gp.x[1,:], youngstress)

#bell = zeros(size(gp.ind.b_bottom[1]))
bell = test_bell(gp,num)
lines!(gp.x[1,:], ((bell/10) .+gp.y[1, 1]))
Makie.save("contact_line_fu.png", fu)

# fv = Figure();
# ax = Axis(fv[1, 1],
#     title="t=$(round(fwd.t[end], digits=3))",
#     xlabel="x",
#     ylabel="y")
# hm = heatmap!(gv.x[1, :], gv.y[:, 1], fwdL.v[end, :, :]')
# Colorbar(fv[1, 2], hm, label="v")
# contour!(gp.x[1, :], gp.y[:, 1], tracer',
#     levels=0:0, color=:red, linewidth=3);
# Makie.save("contact_line_fv.png", fv)

# fp = Figure()
# ax = Axis(fp[1, 1],
#     title="t=$(round(fwd.t[end], digits=3))",
#     xlabel="x",
#     ylabel="y")
# hm = heatmap!(gp.x[1, :], gp.y[:, 1], fwdL.p[end, :, :]')
# Colorbar(fp[1, 2], hm, label="p")
# contour!(gp.x[1, :], gp.y[:, 1], tracer',
# levels=0:0, color=:red, linewidth=3);
# Makie.save("contact_line_fp.png", fp)