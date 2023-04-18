using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

Lx = 2
Ly = 1
n = 64
nx = Lx * n
ny = n

x = LinRange(-Lx/2, Lx/2, nx+1)
y = LinRange(-Ly/2, Ly/2, ny+1)

nu = 0.0044

num = Numerical(
    case = "Monte",
    x = x,
    y = y,
    Re = 1000.0,
    CFL = 0.8,
    max_iterations = 1000,
    u_inf = 1.0,
    v_inf = 0.0,
    save_every = 10,
    ϵ = 0.05,
    NB = 12,
    subdomains = 2,
    overlaps = 1
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

uy = zeros(ny)
uy[ny÷4+2:end] .= 1.0/(4nu) * ((3*Ly/8) ^ 2 .- (gp.y[ny÷4+2:end,1] .- Ly./8.) .^ 2)
uy ./= findmax(uy)[1]

phL.u .= num.u_inf
phL.v .= num.v_inf

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd, fwdS, fwdL,
    BC_uL = Boundaries(
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet),
        left = Boundary(t = dir, f = dirichlet, val=uy),
        right = Boundary(t = neu, f = neumann)
    ),
    BC_vL = Boundaries(
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet),
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = neu, f = neumann)
    ),
    BC_pL = Boundaries(
        left = Boundary(t = neu, f = neumann),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = neu, f = neumann),
        top = Boundary(t = neu, f = neumann)
    ),
    BC_u = Boundaries(
    ),
    stefan = false,
    advection = false,
    heat = false,
    heat_convection = false,
    navier_stokes = true,
    ns_advection = false,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xlabel="x", ylabel="y")  # customized as you see fit
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1)  # customized as you see fit
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], fwdL.p[end,:,:]')
contour!(gp.x[1,:], gp.y[:,1], fwd.u[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)


prefix = "/Users/alex/Documents/PhD/conferences/2023/Montestigliano_2023/project/"
suffix = "_Re$(num.Re)"

limx = num.x[end]
limy = num.y[end]

make_video(gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(-limy,limy), minv=0.0, maxv=1.0)
make_video(gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx,limx), limitsy=(-limy-num.Δ/2,limy+num.Δ/2))
make_video(gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx,limx), limitsy=(-limy,limy))
