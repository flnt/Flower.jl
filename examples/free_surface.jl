using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)


L0 = 2.0
n = 128

x = [-1.5, -0.5, 0.5, 1.5]
y = [-1.5, -0.5, 0.5, 1.5]

x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

num = Numerical(
    # case = "Mullins_cos",
    # case = "Jet",
    case = "Ellipse",
    # case = "Cylinder",
    x = x,
    y = y,
    Re = 1.0,
    CFL = 1.0,
    max_iterations = 2000,
    u_inf = 0.0,
    v_inf = 0.0,
    shifted = 1e-8,
    A = -0.05,
    # A = -0.2,
    N = 2,
    save_every = 10,
    σ = 200.0,
    ϵ = 0.05,
    R = 0.2,
)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

gp.u .*= -1.
# gp.u .+= 0.5
phL.u .= num.u_inf
phL.v .= num.v_inf

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
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
        # bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
        # top = Boundary(t = dir, f = dirichlet, val = 0.0),
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
    show_every = 1
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

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/stokes/free_surface/"
suffix = "_$(num.case)_2"

make_video(num, fwd, gu, "u"; title_prefix=prefix,
        title_suffix=suffix, framerate=50, limitsx=(-lim-num.Δ/2,lim+num.Δ/2), limitsy=(-lim,lim))
make_video(num, fwd, gv, "v"; title_prefix=prefix,
        title_suffix=suffix, framerate=50, limitsx=(-lim,lim), limitsy=(-lim-num.Δ/2,lim+num.Δ/2))
make_video(num, fwd, gp, "p"; title_prefix=prefix,
        title_suffix=suffix, framerate=50, limitsx=(-lim,lim), limitsy=(-lim,lim))#, minv = -1e-4, maxv = 1e-4)
