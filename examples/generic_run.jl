using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ratio = 1
L0 = 16.
nx = 128
ny = ratio * nx

x = LinRange(-L0/2, L0/2, nx+1)
y = LinRange(-ratio*L0/2, ratio*L0/2, ny+1)

num = Numerical(case = "Sphere",
    CFL = 0.1,
    TEND = 1.0,
    x = x,
    y = y,
    R = 1.56,
    save_every = 1,
    ϵ = 0.05,
    θd = 0.0,
    T_inf = -0.5
)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    stefan = true,
    advection = true,
    heat = true,
    heat_solid_phase = false,
    heat_liquid_phase = true,
    verbose = true,
    show_every = 1,
)

tcks = -ratio*L0/2:2:ratio*L0/2
lim = L0 / 2

fsolid = Figure(resolution = (1600, 1000))
colsize!(fsolid.layout, 1, Aspect(1, 1.0))
ax = Axis(fsolid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phS.T')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fsolid)

# fsolid = current_figure()

fliquid = Figure(resolution = (1600, 1000))
colsize!(fliquid.layout, 1, Aspect(1, 1.0))
ax = Axis(fliquid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.T')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
contour!(gp.x[1,:], gp.y[:,1], fwd.usave[1,:,:]', levels = 0:0, color=:blue, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fliquid)

# fliquid = current_figure()
