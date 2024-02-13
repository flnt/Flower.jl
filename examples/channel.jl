using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/contact_line/channel/"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

L0x = 1.0
L0y = 1.0
n = 64

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
y = collect(LinRange(-L0y / 2, L0y / 2, n + 1))
num = Numerical(
    case = "Nothing",
    x = x,
    y = y,
    Re = 1.0,
    CFL = 1.0,
    max_iterations = 500,
    save_every = 10,
    ϵ = 0.05,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

gp.LS[1].u .= 1.0

λ = 1e-2
R = L0y / 2.0

uPoiseuille = R^2 .- gu.y[:,1].^2

phL.u .= 1.0
phL.v .= 0.0

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left = Dirichlet(val = copy(uPoiseuille)),
        bottom = Navier(λ = 5λ),
        top = Dirichlet(),
    ),
    BC_vL = Boundaries(
        left = Dirichlet(),
        bottom = Dirichlet(),
        top = Dirichlet(),
    ),
    BC_pL = Boundaries(
        right = Dirichlet(),
    ),
    BC_int = [Wall()],
    time_scheme = CN,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
)

tcks = -0.5:0.1:0.5

fu = Figure(size = (1200, 1000))
ax = Axis(fu[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding = 0)
colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fu)

fv = Figure(size = (1200, 1000))
ax = Axis(fv[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding = 0)
colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fv)

fp = Figure(size = (1200, 1000))
ax = Axis(fp[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.p')
cbar = fp[1,2] = Colorbar(fp, hmap, labelpadding = 0)
colsize!(fp.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fp.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fp)

utop = vecb_T(phL.uD, gu)
ubottom = vecb_B(phL.uD, gu)

fpr = Figure(size = (1600, 1000))
colsize!(fpr.layout, 1, Aspect(1, 1.0))
ax  = Axis(fpr[1,1], aspect = 1, xlabel = L"U _ x", ylabel = L"y",
    xtickalign = 0,  ytickalign = 0, yticks = tcks)
lines!(vcat(ubottom[end-10], phL.u[:,end-10], utop[end-10]), vcat(-0.5, gu.y[:,1], 0.5), linewidth = 3)
limits!(ax, 0.0, 1.1 * maximum(phL.u[:,end-10]), -0.5, 0.5)
resize_to_layout!(fpr)