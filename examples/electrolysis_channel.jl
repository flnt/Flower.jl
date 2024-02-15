using Revise
using Flower
using Printf

#From channel.jl

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

prefix="/local/home/pr277828/flower/test/"




#Khalighi 2023

L0 = 1e-4 
n = 64
# x = LinRange(-L0/2, L0/2, n+1)
# y = LinRange(-L0/2, L0/2, n+1)

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)

# xb=x
# print("test",length(x))



function f(x)
    v_inlet=6.7e-4
    # vPoiseuille=8/3*v_inlet*gv.x[1,:]/L0*(1-gv.x[1,:]/L0)
    return 8/3*v_inlet*x/L0*(1-x/L0)
end


# function dirichlet_bcs!(gp, D)
#     @unpack x, nx, ind = gp
#     # @inbounds @threads for II in ind.
#     #     D[II] = f(x[II])
#     # end
#     for i = 1:nx
#         D[i]=f(x[i])
#     end
# end



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

#Levelset 1 everywhere
gp.LS[1].u .= 1.0 

#init v=vPoiseuille
vPoiseuille = zeros(gv)
# dirichlet_bcs!(gv, vPoiseuille)

@unpack x, nx, ny, ind = gv
   
# for i = 1:nx
#     vPoiseuille[i]=f(x[i])
# end

vPoiseuille=f.(x)

# print(vPoiseuille)

# myprint(x) = print(round.(x; sigdigits=2))
# myprint(vPoiseuille)
# print(string.(round.(vPoiseuille; digits=2)))

# print(Printf.format.(Ref(Printf.Format("%.2e")), vPoiseuille))
# print(length(vPoiseuille),length(x))

# for j = 1:ny
#     phL.v[:,j]=vPoiseuille
# end

phL.v .=vPoiseuille
phL.u .= 0.0
phL.p .= 0.0

vPoiseuilleb=f.(gv.x[1,:])

# print(vPoiseuilleb)
# break 


#Neumann by default?


# λ = 1e-2
# R = L0y / 2.0

# uPoiseuille = R^2 .- gu.y[:,1].^2

# phL.u .= 1.0
# phL.v .= 0.0

# @time run_forward(
#     num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
#     BC_uL = Boundaries(
#         left = Dirichlet(val = copy(uPoiseuille)),
#         bottom = Navier(λ = 5λ),
#         top = Dirichlet(),
#     ),
#     BC_vL = Boundaries(
#         left = Dirichlet(),
#         bottom = Dirichlet(),
#         top = Dirichlet(),
#     ),
#     BC_pL = Boundaries(
#         right = Dirichlet(),
#     ),
#     BC_int = [Wall()],
#     time_scheme = CN,
#     navier_stokes = true,
#     ns_advection = true,
#     ns_liquid_phase = true,
#     verbose = true,
#     show_every = 1,
# )




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
        #not the bottom?
        # top = Dirichlet(val = copy(vPoiseuilleb)), 
        top=Neumann(val=0.0),

    ),
    BC_pL = Boundaries(
        left  = Neumann(val=0.0),
        right = Neumann(val=0.0),
        #not the bottom?
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