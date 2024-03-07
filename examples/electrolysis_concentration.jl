using Revise
using Flower
using Printf
using PrettyTables

#From channel.jl and heat_convection.jl
#Khalighi 2023: "Hydrogen bubble growth in alkaline water electrolysis: An immersed boundary simulation study"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

prefix="/local/home/pr277828/flower/"

folder="test"

prefix *= "/"*folder*"/"

isdir(prefix) || mkdir(prefix)

L0 = 1e-4 
n = 64

h0 = 0.25*L0 #TODO h0

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)

function f(x,v_inlet)
    # v_inlet=6.7e-4
    return 8/3*v_inlet*x/L0*(1-x/L0)
end


####################################################################################################
# Sessile from sessile.jl
####################################################################################################

Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)

θe = 45
max_its = 3000

if max_its <= 100
    save_every = 1
else
    save_every = max_its÷100
end
n_ext = 10
CFL = 0.5

_θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
println("θe = $(_θe)")


# num = Numerical(
#         case = "Planar",
#         x = x,
#         y = y,
#         Re = 1.0,
#         CFL = CFL,
#         max_iterations = max_its,
#         u_inf = 0.0,
#         v_inf = 0.0,
#         save_every = max_its÷100,
#         # save_every = 1,
#         reinit_every = 3,
#         nb_reinit = 2,
#         δreinit = 10.0,
#         σ = 1.0,
#         ϵ = 0.05,
#         n_ext_cl = n_ext,
#         NB = 24
#     )

####################################################################################################


rho=1258
#TODO need for 80°C, check with other study
rhoH2=0.069575 #"0.7016E-01" in \citet{cohnTABLETHERMODYNAMICPROPERTIES1963} 
#Linear interpolation between 350 and 360
# 350	360	353		B	0.13841	353	350	360
# 7.02E-02	6.82E-02		-0.00194999999999999	A	-0.000194999999999999	0.069575	0.07016	0.06821



mu=6.7e-7
i0=1.0
temp0=353.0
pres0=1e5
sigma=7.7e-2
KOHwtpercent=30
phi_ele1=-0.6
alphac=0.5
alphaa=0.5
DH2=5.8e-9
DKOH=3.2e-9
DH2O=3.2e-9
c0_H2=1.6e-1
c0_KOH=6.7e3
c0_H2O=4.9e4
Ru=8.314
Faraday = 9.64853321233100184e4 #C⋅mol−1
MWH2 = 2.01568e-3 #kg/mol

# elec_cond=1 #TODO

elec_cond=2*Faraday^2*c0_KOH*DKOH/(Ru*temp0)


print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")


v_inlet=6.7e-4
Re=rho*v_inlet*L0/mu

diffusion_coeff=[DH2, DKOH, DH2O]
concentration0=[0.16, 6700, 49000]

# concentration0=(0.16, 6700, 49000)
# diffusion_coeff=(DH2, DKOH, DH2O)


nb_transported_scalars=3

if length(concentration0)!=nb_transported_scalars
    print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
    @error ("nb_transported_scalars")
end

if length(diffusion_coeff)!=nb_transported_scalars
    print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
    @error ("nb_transported_scalars")
end

print(@sprintf "Re = %.2e\n" Re)

# print(@sprintf "Transported scalars : %2i\n" nb_transported_scalars)

print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)

# pretty_table(concentration0'; header = ["cH2", "cKOH", "cH2O"])

# pretty_table(diffusion_coeff'; header = ["DH2", "DKOH", "DH2O"])

# pretty_table(vcat(hcat("D",diffusion_coeff'),hcat("c",concentration0')); header = ["","H2", "KOH", "H2O"])

# hl = Highlighter((d,i,j)->d[i,j][1]*d[i,j][2] < 0, crayon"red")

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold blue")


pretty_table(vcat(hcat("Diffusion coef",diffusion_coeff'),hcat("Concentration",concentration0')); header = ["","H2", "KOH", "H2O"], highlighters=hl)




num = Numerical(
    case = "Planar", #"Nothing",
    x = x,
    y = y,
    Re = 1.0,
    CFL = 1.0,
    # max_iterations = 500, #max_its
    max_iterations = 5,
    save_every = 10,
    ϵ = 0.05,
    nb_transported_scalars=3,
    electrolysis = true,
    concentration0=concentration0, 
    diffusion_coeff=diffusion_coeff,
    temp0=temp0,
    i0=i0,
    phi_ele1=phi_ele1,
    alphac=alphac,
    alphaa=alphaa,
    Ru=Ru,
    Faraday=Faraday,
    MWH2=MWH2,

    u_inf = 0.0,
    v_inf = 0.0,
    # save_every = max_its÷100,
    reinit_every = 3,
    nb_reinit = 2,
    δreinit = 10.0,
    # σ = 1.0,
    # ϵ = 0.05,
    n_ext_cl = n_ext,
    NB = 24
    
)

    # TEND=7.3,
    # T_inf=353,
    # case="Electrolysis_concentration"

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

####################################################################################################
# Sessile from sessile.jl
####################################################################################################

# r = 0.5
# gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
# gp.LS[1].u .*= -1.0

# phL.u .= 0.0
# phL.v .= 0.0
####################################################################################################


phL.T .= 0.

#Initialization

#Levelset 1 everywhere
gp.LS[1].u .= 1.0 

vPoiseuille = zeros(gv)
@unpack x, nx, ny, ind = gv
vPoiseuille=f.(x,v_inlet)

phL.v .=vPoiseuille
phL.u .= 0.0
phL.p .= 0.0

vPoiseuilleb=f.(gv.x[1,:],v_inlet)

#Neumann by default?

#TODO pressure left and right BC not mentioned in the article Khalighi 2023


#TODO need to iterate more for potential since phiele=0 initially?

phi_ele=gv.x[1,:] .*0.0
# eta = phi_ele1 .-phi_ele
#TODO precision: number of digits
# i_current=i0*(exp(alphaa*Faraday*eta/(Ru*temp0))-exp(-alphac*Faraday*eta/(Ru*temp0)))
i_current=butler_volmer_no_concentration.(alphaa,alphac,Faraday,i0,phi_ele,phi_ele1,Ru,temp0)

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
    BC_int = [FreeSurface()], #[Wall()],


    BC_trans_scal = (
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top = Neumann(),
        left = Neumann(val=i_current/(2*Faraday*DH2)),
        right = Dirichlet(val = concentration0[1]),
        int=Neumann(val=0.0)), #H2
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top = Neumann(),
        left = Neumann(val=i_current/(2*Faraday*DKOH)),
        right = Neumann(),
        int=Neumann(val=0.0)), #KOH
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top = Neumann(),
        left = Neumann(val=-i_current/(Faraday*DH2O)),
        right = Neumann(),
        int=Neumann(val=0.0)) #H2O
    ),

    BC_phi_eleL = BoundariesInt(
        left = Neumann(val=-i_current/elec_cond),
        right = Dirichlet(),
        bottom = Dirichlet(),
        top=Neumann(val=0.0),
        int=Neumann(val=0.0),
    ),


    time_scheme = CN, #or FE?
    electrolysis = true,

    auto_reinit = true,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    electrolysis_convection = true,  
    electrolysis_liquid_phase = true,
    electrolysis_phase_change = true,

    # ns_advection = false, #?
    save_length = true,
)

####################################################
# BC_int = [FreeSurface()],
# time_scheme = FE,
# auto_reinit = true,
# navier_stokes = true,
# ns_advection = false,
# ns_liquid_phase = true,
# verbose = true,
# show_every = 1,
# save_length = true,
####################################################


# tcks = -0.5:0.1:0.5

# fu = Figure(size = (1200, 1000))
# ax = Axis(fu[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding = 0)
# colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fu)

# fv = Figure(size = (1200, 1000))
# ax = Axis(fv[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
# hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
# cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding = 0)
# colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fv)

# fp = Figure(size = (1200, 1000))
# ax = Axis(fp[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
# hmap = heatmap!(gp.x[1,:], gp.y[:,1], phL.p')
# cbar = fp[1,2] = Colorbar(fp, hmap, labelpadding = 0)
# colsize!(fp.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fp.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fp)

# utop = vecb_T(phL.uD, gu)
# ubottom = vecb_B(phL.uD, gu)

# fpr = Figure(size = (1600, 1000))
# colsize!(fpr.layout, 1, Aspect(1, 1.0))
# ax  = Axis(fpr[1,1], aspect = 1, xlabel = L"U _ x", ylabel = L"y",
#     xtickalign = 0,  ytickalign = 0, yticks = tcks)
# lines!(vcat(ubottom[end-10], phL.u[:,end-10], utop[end-10]), vcat(-0.5, gu.y[:,1], 0.5), linewidth = 3)
# limits!(ax, 0.0, 1.1 * maximum(phL.u[:,end-10]), -0.5, 0.5)
# resize_to_layout!(fpr)

xlabel = L"x \left(\mu m \right)"
ylabel = L"y \left(\mu m \right)"

xscale = 1e-6
yscale = xscale

# xticks = 0:20:num.L0/xscale
# yticks = 0:20:num.L0/yscale

xticks = 0:20:100
yticks = 0:20:100

velscale = 1e-4 


make_video_vec(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u",
        title_suffix="", framerate=240, 
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=velscale, scalelabel=L"\times 10^{-4}", xticks=xticks, yticks=yticks, scalticks=0:2:10)
make_video_vec(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v",
        title_suffix="", framerate=240, 
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=velscale, scalelabel=L"\times 10^{-4}", xticks=xticks, yticks=yticks, scalticks=0:2:10)
make_video_vec(num, gp, fwd.u, fwdL.T; title_prefix=prefix*"T",
        title_suffix="", framerate=240, 
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=1, xticks=xticks, yticks=yticks)

make_video_vec(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p",
title_suffix="", framerate=240, 
xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=1, xticks=xticks, yticks=yticks)

make_video_vec(num, gp, fwd.u, fwdL.trans_scal[:,:,:,1]; title_prefix=prefix*"concentration_H2",
title_suffix="", framerate=240, 
xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=c0_H2, xticks=xticks, yticks=yticks, scalticks=0:500:3000)

make_video_vec(num, gp, fwd.u, fwdL.trans_scal[:,:,:,2]; title_prefix=prefix*"concentration_KOH",
title_suffix="", framerate=240, 
xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=c0_KOH, xticks=xticks, yticks=yticks)

make_video_vec(num, gp, fwd.u, fwdL.trans_scal[:,:,:,3]; title_prefix=prefix*"concentration_H2O",
title_suffix="", framerate=240, 
xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=c0_H2O, xticks=xticks, yticks=yticks)

make_video_vec(num, gp, fwd.u, fwdL.phi_ele; title_prefix=prefix*"phi_ele",
title_suffix="", framerate=240, 
xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=1, xticks=xticks, yticks=yticks)

#Vector plot

fvector = Figure(size = (800, 800))
# Axis(f[1, 1], 
# # backgroundcolor = "black"
# )



ax = Axis(fvector[1,1], aspect = DataAspect(), xlabel = L"x \left(\mu m \right)", ylabel = L"y \left(\mu m \right)",
        xticks = xticks, yticks = xticks)

# ax  = Axis(fpr[1,1], aspect = 1, xlabel = L"U _ x", ylabel = L"y",
#     xtickalign = 0,  ytickalign = 0, yticks = tcks)



us, vs = interpolate_regular_grid(gp,fwdL)

# print(us)

# print(vs)

arrowscale = 0.3
arrowscale = 1e-4


#https://docs.makie.org/stable/reference/plots/arrows/
arrows!(gp.x[1,:] ./xscale, gp.y[:,1] ./yscale, us[1,:,:], vs[1,:,:], 
arrowsize = 10, 
lengthscale = arrowscale,
# arrowcolor = strength, 
# linecolor = strength,
)

fvector





Makie.save(prefix*"vector.pdf", fvector)



# testfig = Figure(size = (800, 800))
# Axis(testfig[1, 1], backgroundcolor = "black")

# xs = LinRange(0, 2pi, 20)
# ys = LinRange(0, 3pi, 20)
# us = [sin(x) * cos(y) for x in xs, y in ys]
# vs = [-cos(x) * sin(y) for x in xs, y in ys]
# strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))

# arrows!(xs, ys, us, vs, arrowsize = 10, lengthscale = 0.3,
#     arrowcolor = strength, linecolor = strength)

# testfig

# Makie.save(prefix*"testfig.pdf", testfig)
