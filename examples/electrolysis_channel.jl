using Revise
using Flower
using Printf

#From channel.jl

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

prefix="/local/home/pr277828/flower/test/"


#Khalighi 2023: "Hydrogen bubble growth in alkaline water electrolysis: An immersed boundary simulation study"

L0 = 1e-4 
n = 64

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)

function f(x)
    v_inlet=6.7e-4
    return 8/3*v_inlet*x/L0*(1-x/L0)
end

rho=1258
mu=6.7e-7
# i0=1
# theta0=353
# pres0=1e5
# sigma=7.7e-2
# KOHwtpercent=30
# phi1=-0.6
# alpha_c=0.5
# alpha_a=0.5
# DH2=5.8e-9
# DKOH=3.2e-9
# DH2O=3.2e-9
# C0H2=1.6e-1
# C0KOH=6.7e3
# C0H2O=4.9e4
Re=rho*v_inlet*L0/mu

print(@sprintf "Re = %.2e\n" Re)


num = Numerical(
    case = "Nothing",
    x = x,
    y = y,
    Re = Re,
    CFL = 1.0,
    max_iterations = 500,
    save_every = 10,
    ϵ = 0.05,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

#Initialization

#Levelset 1 everywhere
gp.LS[1].u .= 1.0 

vPoiseuille = zeros(gv)
@unpack x, nx, ny, ind = gv
vPoiseuille=f.(x)

phL.v .=vPoiseuille
phL.u .= 0.0
phL.p .= 0.0

vPoiseuilleb=f.(gv.x[1,:])

#Neumann by default?

#TODO pressure left and right BC not mentioned in the article Khalighi 2023

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left = Dirichlet(),
        right = Dirichlet(),
        bottom = Dirichlet(),
        top=Neumann(val=0.0),

    ),
    BC_vL = Boundaries(
        left = Dirichlet(),
        right=Dirichlet(),
        bottom = Dirichlet(val = copy(vPoiseuilleb)),
        top=Neumann(val=0.0),

    ),
    BC_pL = Boundaries(
        left  = Neumann(val=0.0),
        right = Neumann(val=0.0),
        bottom=Neumann(val=0.0),
        top= Dirichlet(),
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


make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field",
        title_suffix="", framerate=240)
make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field",
        title_suffix="", framerate=240)
make_video(num, gp, fwd.u, fwdL.T; title_prefix=prefix*"T_field",
        title_suffix="", framerate=240)

make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field",
title_suffix="", framerate=240)