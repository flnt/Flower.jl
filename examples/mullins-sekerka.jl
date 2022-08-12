using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

L0 = 2.0
n = 64
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

num = Numerical(case = "Mullins_cos",
    x = x,
    y = y,
    T_inf = 1.0,
    CFL = 0.5,
    TEND = 0.5,
    A = -0.05,
    N = 2,
    ϵ = 0.05
)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    periodic_x = true,
    BC_TL = Boundaries(
        top = Boundary(t = dir, f = dirichlet, val = -num.T_inf),
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_TS = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_uL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_vL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_pL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = dir, f = dirichlet, val = 0.),
        top = Boundary(t = dir, f = dirichlet, val = 0.),
    ),
    stefan = true,
    advection = true,
    heat = true,
    heat_solid_phase = true,
    heat_liquid_phase = true,

    heat_convection = false,
    navier_stokes = true,
    ns_advection = false,
    ns_liquid_phase = true,
    ns_solid_phase = false,

    verbose = true,
    show_every = 1
)

f = heatmap!(gp.x[1,:], gp.y[:,1], (phL.T+phS.T)', colormap= :ice)
f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

for j = num.max_iterations÷5:num.max_iterations÷5:num.max_iterations
    f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[j,:,:]', levels = 0:0, color=:black, linewidth = 3);
end

f = current_figure()

pref = "/Users/alex/Documents/PhD/Cutcell/New_ops/heat/"
suff = ""
make_video(num, fwd, gu, "u"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.5, maxv = 0.5)
make_video(num, fwd, gv, "v"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.3, maxv = 0.3)
make_video(num, fwd, gp, "p"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))#, minv = -0.03, maxv = 0.03)
make_video(num, fwd, gp, "T"; title_prefix=pref,
        title_suffix=suff, framerate=20)
make_video(num, fwd, gp, "ϕ"; title_prefix=pref,
        title_suffix=suff, framerate=20)
make_video(num, fwd, gu, "ucorr"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))
make_video(num, fwd, gv, "vcorr"; title_prefix=pref,
        title_suffix=suff, framerate=20, limitsx=(-1.0,1.0), limitsy=(-1.0,1.0))