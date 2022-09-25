using Revise
using Flower

ratio = 1
L0 = 2.
nx = 128
ny = ratio * nx

x = LinRange(-L0/2, L0/2, nx+1)
y = LinRange(-ratio*L0/2, ratio*L0/2, ny+1)

num = Numerical(case = "Sphere",
    CFL = 0.5,
    TEND = 0.1,
    max_iterations = 1,
    x = x,
    y = y,
    R = 0.8,
    save_every = 1,
    T_inf = 1.0,
    ϵ = 0.1,
    θd = 1.0,
    A = -0.2,
    N = 6,
    NB = 2
)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    # BC_TS = Boundaries(
    #     bottom = Boundary(t = dir, f = dirichlet, val = 1.0)
    # ),
    stefan = true,
    heat = true,
    heat_solid_phase = true,
    heat_liquid_phase = false,
    verbose = true,
    show_every = 1,
)

tcks = -ratio*L0/2:2:ratio*L0/2
lim = L0 / 2

fsolid = Figure(resolution = (1600, 1000))
colsize!(fsolid.layout, 1, Aspect(1, 1.0))
ax = Axis(fsolid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phS.T', colorrange=(0, num.θd))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fsolid)

#fsolid = current_figure()

# fliquid = Figure(resolution = (1600, 1000))
# colsize!(fliquid.layout, 1, Aspect(1, 1.0))
# ax = Axis(fliquid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gp.x[1,:], gp.y[:,1], phL.T', colorrange=(0, num.θd))
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# # limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fliquid)

# # fliquid = current_figure()