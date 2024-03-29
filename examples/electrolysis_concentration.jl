using Revise
using Flower
using Printf
using PrettyTables
using Interpolations


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
# max_iter=5
# n = 10
# max_iter=1

max_iter=100

# max_iter=1


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

# Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
# RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
# center(r, θ) = r * cos(π - θ)

# θe = 45
# max_its = 3000

# if max_its <= 100
#     save_every = 1
# else
#     save_every = max_its÷100
# end
# n_ext = 10
# CFL = 0.5

# _θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# # _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# println("θe = $(_θe)")

####################################################################################################

rho=1258.0
#TODO need for 80°C, check with other study
rhoH2=0.069575 #"0.7016E-01" in \citet{cohnTABLETHERMODYNAMICPROPERTIES1963} 
#Linear interpolation between 350 and 360
# 350	360	353		B	0.13841	353	350	360
# 7.02E-02	6.82E-02		-0.00194999999999999	A	-0.000194999999999999	0.069575	0.07016	0.06821


radius=2.5e-5 
# xcoord = 0.0
# ycoord = 0.0
# xcoord = radius + 1e-6
# ycoord=5e-5
xcoord = radius + 2e-6 #or 3e-6? not written in the article
ycoord=5e-5

xcoord = 5e-5
ycoord=5e-5

xcoord = -xcoord
ycoord = -ycoord

mu=6.7e-7
i0=1.0
temperature0=353.0
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
phi_ele0=0.0
CFL=0.5
Re=1.0
TEND=7.3#s

# elec_cond=1 #TODO

elec_cond=2*Faraday^2*c0_KOH*DKOH/(Ru*temperature0)


print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")


v_inlet=6.7e-4
Re=rho*v_inlet*L0/mu

#Concentration
diffusion_coeff=[DH2, DKOH, DH2O]
concentration0=[0.16, 6700, 49000]
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


print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)

# pretty_table(concentration0'; header = ["cH2", "cKOH", "cH2O"])
# pretty_table(diffusion_coeff'; header = ["DH2", "DKOH", "DH2O"])
# pretty_table(vcat(hcat("D",diffusion_coeff'),hcat("c",concentration0')); header = ["","H2", "KOH", "H2O"])
# hl = Highlighter((d,i,j)->d[i,j][1]*d[i,j][2] < 0, crayon"red")

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold blue")

pretty_table(vcat(hcat("Diffusion coef",diffusion_coeff'),hcat("Concentration",concentration0')); 
header = ["","H2", "KOH", "H2O"], highlighters=hl)

num = Numerical(
    case = "Cylinder", #"Sphere", #"Nothing",
    shifted = xcoord,
    shifted_y = ycoord,
    R = radius,
    x = x,
    y = y,
    Re = Re,
    CFL = CFL,
    # max_iterations = 500, #max_its
    max_iterations = max_iter,
    save_every = 10,
    ϵ = 0.05,
    nb_transported_scalars=nb_transported_scalars,
    concentration0=concentration0, 
    diffusion_coeff=diffusion_coeff,
    temperature0=temperature0,
    i0=i0,
    phi_ele0=phi_ele0,
    phi_ele1=phi_ele1,
    alphac=alphac,
    alphaa=alphaa,
    Ru=Ru,
    Faraday=Faraday,
    MWH2=MWH2,
    θd=temperature0,
    eps=1e-12,
    TEND=TEND,
    mu1=mu,
    mu2=mu,
    rho1=rho,
    rho2=rhoH2,

    ####################################
    # Sessile
    ####################################
    # u_inf = 0.0,
    # v_inf = 0.0,
    # # save_every = max_its÷100,
    # reinit_every = 3,
    # nb_reinit = 2,
    # δreinit = 10.0,
    # # σ = 1.0,
    # # ϵ = 0.05,
    # n_ext_cl = n_ext,
    # NB = 24,
    ####################################
)

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


phL.T .= temperature0
# phL.TD .= temperature0
phS.T .= temperature0


#Initialization

figname0=""

#Levelset 1 everywhere
gp.LS[1].u .= 1.0
figname0="no_intfc"

vPoiseuille = zeros(gv)
@unpack x, nx, ny, ind = gv
vPoiseuille=f.(x,v_inlet)

phL.v .=vPoiseuille
phL.u .= 0.0
phL.p .= 0.0

vPoiseuilleb=f.(gv.x[1,:],v_inlet)

for iscal=1:nb_transported_scalars
    phL.trans_scal[:,:,iscal] .= concentration0[iscal]
end

phL.phi_ele .= phi_ele0
# phL.phi_eleD .= 0.0 #TODO  phL.phi_eleD .= 0.0 ?

printstyled(color=:green, @sprintf "\n Initialisation \n")

print_electrolysis_statistics(nb_transported_scalars,phL)


printstyled(color=:green, @sprintf "\n TODO timestep CFL scal, and print \n")


@unpack τ,CFL,Δ,Re,θd=num
print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ CFL CFL*Δ CFL*Δ^2*Re Re θd)
# τ=CFL*Δ/v_inlet
# num.τ=τ
# print(@sprintf "dt %.2e \n" τ)


#Neumann by default?

#TODO pressure left and right BC not mentioned in the article Khalighi 2023


#TODO need to iterate more for potential since phiele=0 initially?

phi_ele=gv.x[1,:] .*0.0
# eta = phi_ele1 .-phi_ele
#TODO precision: number of digits
# i_current=i0*(exp(alphaa*Faraday*eta/(Ru*temperature0))-exp(-alphac*Faraday*eta/(Ru*temperature0)))
i_current=butler_volmer_no_concentration.(alphaa,alphac,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)

print(@sprintf "Butler-Volmer %.2e \n" i_current[1])

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Neumann(val=0.0),
    ),
    BC_vL = Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(val = copy(vPoiseuilleb)),
        top    = Neumann(val=0.0),
    ),
    BC_pL = Boundaries(
        left   = Neumann(val=0.0),
        right  = Neumann(val=0.0),
        bottom = Neumann(val=0.0),
        top    = Dirichlet(),
    ),
    # BC_int = [FreeSurface()], #[Wall()],
    BC_int = [Wall()],

    BC_TL  = Boundaries(
    left   = Dirichlet(val = temperature0),
    right  = Dirichlet(val = temperature0),
    bottom = Dirichlet(val = temperature0),
    top    = Dirichlet(val = temperature0),
    ),

    BC_trans_scal = (
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Neumann(val=i_current/(2*Faraday*DH2)),
        right  = Dirichlet(val = concentration0[1]),
        int    = Neumann(val=0.0)), #H2
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=i_current/(2*Faraday*DKOH)),
        right  = Neumann(),
        int    = Neumann(val=0.0)), #KOH
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=-i_current/(Faraday*DH2O)),
        right  = Neumann(),
        int    = Neumann(val=0.0)) #H2O
    ),

    BC_phi_eleL = BoundariesInt(
        left   = Neumann(val=-i_current/elec_cond),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
    ),


    time_scheme = CN, #or FE?
    electrolysis = true,

    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    electrolysis_convection = true,  
    electrolysis_liquid_phase = true,
    electrolysis_phase_change = true,
    electrolysis_reaction = "Butler_no_concentration",
    # electrolysis_reaction = "nothing",
    adapt_timestep_mode = 1,

    # ns_advection = false,

    # auto_reinit = true,
    # # ns_advection = false, #?
    # save_length = true,
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
# ax  = Axis(fpr[1,1], aspect = DataAspect(), xlabel = L"U _ x", ylabel = L"y",
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

make_video_vec(num, gp, fwd.u, (fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2; title_prefix=prefix*"concentration_H2_normalized",
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


#  #Vector plot

#  fvector = Figure(size = (800, 800))
#  # Axis(f[1, 1], 
#  # # backgroundcolor = "black"
#  # )


#  ax = Axis(fvector[1,1], 
#  # aspect = DataAspect(), 
#  aspect=1,
#  xlabel = L"x \left(\mu m \right)", ylabel = L"y \left(\mu m \right)",
#  xticks = xticks, yticks = xticks)

#  # ax  = Axis(fpr[1,1], aspect = DataAspect(), xlabel = L"U _ x", ylabel = L"y",
#  #     xtickalign = 0,  ytickalign = 0, yticks = tcks)

#  us, vs = interpolate_regular_grid(gp,fwdL)

#  # print(us)

#  # print(vs)

#  arrowscale = 1.0 #0.3
#  # arrowscale = 1e-4
#  arrowscale = 1.0/v_inlet #1.e4

#  strength = vec(sqrt.(us[1,:,:] .^ 2 .+ vs[1,:,:] .^ 2))

#  #https://docs.makie.org/stable/reference/plots/arrows/
#  arrows!(gp.x[1,:] ./xscale, gp.y[:,1] ./yscale, us[1,:,:], vs[1,:,:], 
#  arrowsize = 10, 
#  lengthscale = arrowscale,
#  # arrowcolor = strength, 
#  # linecolor = strength,
#  )

#  fvector

#  Makie.save(prefix*"vector.pdf", fvector)

us, vs = interpolate_regular_grid(gp,fwdL)


function plot_velocity_vector(iter,us,vs)
    #Vector plot

    fvector = Figure(size = (800, 800))
    # Axis(f[1, 1], 
    # # backgroundcolor = "black"
    # )


    ax = Axis(fvector[1,1], 
    # aspect = DataAspect(), 
    aspect=1,
    xlabel = L"x \left(\mu m \right)", ylabel = L"y \left(\mu m \right)",
    xticks = xticks, yticks = xticks)

    # ax  = Axis(fpr[1,1], aspect = DataAspect(), xlabel = L"U _ x", ylabel = L"y",
    #     xtickalign = 0,  ytickalign = 0, yticks = tcks)


    # print(us)

    # print(vs)

    arrowscale = 1.0 #0.3
    # arrowscale = 1e-4
    arrowscale = 1.0/v_inlet #1.e4

    strength = vec(sqrt.(us[iter,:,:] .^ 2 .+ vs[iter,:,:] .^ 2))

    #https://docs.makie.org/stable/reference/plots/arrows/
    arrows!(gp.x[1,:] ./xscale, gp.y[:,1] ./yscale, us[iter,:,:], vs[iter,:,:], 
    arrowsize = 10, 
    lengthscale = arrowscale,
    # arrowcolor = strength, 
    # linecolor = strength,
    )

    fvector

    Makie.save(prefix*"vector.pdf", fvector)

end

# plot_velocity_vector(1,us,vs)
# plot_velocity_vector(max_iter,us,vs)







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


# lim = num.L0 / 2
# fH2 = Figure(size = (1600, 1000))
# ax = Axis(fH2[1,1], aspect = DataAspect(), xticks = -4:1:4, yticks = -4:1:4)
# contourf!(gp.x[1,:], gp.y[:,1], phL.T', colormap=:dense, colorrange=(0.2, 1.0))
# contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
# contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
# limits!(ax, -lim, lim, -lim, lim)
# colsize!(fH2.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fH2.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fH2)

# ax, hm = contourf(fig[2, 2][1, 1], xs, ys, zs,
#     colormap = :Spectral, levels = [-1, -0.5, -0.25, 0, 0.25, 0.5, 1])
# Colorbar(fig[2, 2][1, 2], hm, ticks = -1:0.25:1)


phL.trans_scal[:,:,1] .= (1500.0+1)*c0_H2

# lim = num.L0 / 2
fH2 = Figure(size = (1600, 1000))
ax = Axis(fH2[1,1], aspect = DataAspect(), xticks = xticks, yticks = yticks)

co=contourf!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, (phL.trans_scal[:,:,1]'.-c0_H2)./c0_H2, 
# colormap=:dense, 
# colorrange=(0.2, 1.0),
levels = 0:500:3000,
)

# ax,co=contourf(gp.x[1,:]./xscale, gp.y[:,1]./yscale, (phL.trans_scal[:,:,1]'.-c0_H2)./c0_H2, 
# # colormap=:dense, 
# # colorrange=(0.2, 1.0),
# # levels = 0:500:3000,
# levels = [0,500,1000,1500,2000,2500,3000],
# )

Colorbar(fH2[1, 2], co, ticks = 0:500:3000)
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, gp.LS[1].u', levels = 0:0, color=:red, linewidth = 5);
contour!(gp.x[1,:]./xscale, gp.y[:,1]./yscale, fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 5, linestyle=:dot);
# limits!(ax, -lim, lim, -lim, lim)
limits!(ax, 0.0, num.L0/xscale, 0.0, num.L0/yscale)
# colsize!(fH2.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fH2.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fH2)

Makie.save(prefix*"fH2.pdf", fH2)




# function splot(u, v)
#     nx, ny = size(u)
#     x, y = 1:nx, 1:ny
#     intu, intv = linear_interpolation((x,y), u), linear_interpolation((x,y), v)
#     f(x) = Point2f(intu(x...), intv(x...))
#     return streamplot(f, x, y, colormap=:magma)
# end

function splot(u, v, gp)
    nx, ny = size(u)
    x, y = gp.x[1,:], gp.y[:,1]
    intu, intv = linear_interpolation((x,y), u), linear_interpolation((x,y), v)
    f(x) = Point2f(intu(x...), intv(x...))
    return streamplot(f, x, y, colormap=:magma)
end


# # Streamplot example function
# struct FitzhughNagumo{T}
#     ϵ::T
#     s::T
#     γ::T
#     β::T
# end
# P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

# x = -1.5:0.1:1.5
# nx = length(x)
# u, v = zeros(nx, nx), zeros(nx, nx)

# for (j, xj) in enumerate(x)
#     for (i, xi) in enumerate(x)
#         u[i, j], v[i, j] = (xi-xj-xi^3+P.s)/P.ϵ, P.γ*xi-xj + P.β
#     end
# end

# Eus, Evs = interpolate_regular_grid(gp,fwdL,phL.Eu, phL.Ev)
Eus, Evs = interpolate_regular_grid(gp,fwdL,fwdL.Eu, fwdL.Ev)

fig, ax, pl = splot(Eus[1,:,:], Evs[1,:,:], gp)
plot_velocity_vector(1,us,vs)

fig
Makie.save(prefix*"field_lines.pdf",fig)

function make_video_ele(
    num, grid, field_u, field = nothing,
     vecu = nothing, vecv=nothing;
    title_prefix = "video",
    title_suffix = "",
    xlabel = L"x",
    ylabel = L"y",
    colormap = :viridis,
    sz = (1600, 1000),
    minv = 0.0,
    maxv = 0.0,
    limitsx = false,
    limitsy = false,
    var = 1,
    framerate = 24,
    step = 1,
    step0 = 1,
    stepf = size(field_u, 2),
    xscale = 1, 
    yscale = 1, 
    xticks = nothing,
    yticks = nothing,
    scalscale = 1, 
    scalelabel = nothing,
    scalticks = nothing,
    )

    x = grid.x[1,:] ./xscale
    y = grid.y[:,1] ./yscale

    u = field_u[:,step0:stepf,:,:]
    plot_hmap = true
    if isnothing(field)
        plot_hmap = false
    else
        if length(size(field)) == 2
            if var == 1
                z = reshape(
                    field[step0:stepf, 1:grid.ny*grid.nx], 
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            else
                z = reshape(
                    field[step0:stepf, grid.ny*grid.nx+1:2*grid.ny*grid.nx],
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            end
        else
            z = field[step0:stepf,:,:]
        end

        z = z ./scalscale
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64, Float64}) || isa(limitsx, Tuple{Int, Int})
        lx = limitsx
    else
        lx = (min(x...) - grid.dx[1,1] / 2, max(x...) + grid.dx[1,end] / 2)
    end
    if isa(limitsy, Tuple{Float64, Float64}) || isa(limitsy, Tuple{Int, Int})
        ly = limitsy
    else
        ly = (min(y...) - grid.dy[1,1] / 2, max(y...) + grid.dy[end,1] / 2)
    end

    obs = Observable{Int32}(1)
    iterator = range(0, size(u, 2) - 1, step=step)

    # fig = Figure(size = sz)
    fig = Figure()
    ax  = Axis(fig[1,1], 
    # aspect=DataAspect(), 
    aspect=1,
    xlabel = xlabel, ylabel = ylabel,
    xticks=xticks, yticks=yticks,
        xtickalign = 0,  ytickalign = 0)
    if plot_hmap
        print("scalticks", scalticks, "xticks", xticks, "yticks", yticks,"\n")

        if !var_colorrange
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap, 
                colorrange = (minv, maxv))
        else
            # if !isnothing(scalticks)
            #     print("scalticks", scalticks)

            #     hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap = colormap, ticks=scalticks)
            # else
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap)
            # end
        end
    end
    if !plot_hmap
        contour!(x, y, @lift(u[1,$obs,:,:]'), levels = -10:grid.dx[1,1]:10, linewidth = 2.0)
    end
    for iLS in 1:num.nLS
        contour!(x, y, @lift(u[iLS,$obs,:,:]'), levels = [0.0], color = :red, linewidth = 3.0)
    end

    #TODO
    if !isnothing(vecu)
        arrows!(grid.x[1,:], grid.y[:,1], vecu[1,:,:], vecv[1,:,:], 
        arrowsize = 10, 
        lengthscale = 0.3,
        # arrowcolor = strength, 
        # linecolor = strength,
        )
    end

    # contourf!(x, y, u[2,1,:,:]', levels = 0:0, color = :red, linewidth = 3, extendlow = :gray);
    if plot_hmap
        if !isnothing(scalticks)
            print("scalticks", scalticks)
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0, ticks=scalticks)
        else
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0)
        end

    end

    if !isnothing(scalelabel)
        # Label(fig[1, 1, Top()], halign = :right, scalelabel)
        Label(fig[1, 2, Top()], halign = :center, scalelabel)
    end

    # limits!(ax, lx[1], lx[2], ly[1], ly[2])
    # colgap!(fig.layout, 10)
    # rowgap!(fig.layout, 10)
    # colsize!(fig.layout, 1, widths(ax.scene.viewport[])[1])
    # rowsize!(fig.layout, 1, widths(ax.scene.viewport[])[2])
    # resize_to_layout!(fig)

    vid = record(
            fig, title_prefix*title_suffix*".mp4", iterator;
            framerate = framerate
        ) do it
        obs[] = it+1
    end

    return vid
end

make_video_elec(num, gu, fwd; title_prefix=prefix*"elec",
        title_suffix="", framerate=240, 
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=velscale, scalelabel=L"\times 10^{-4}", xticks=xticks, yticks=yticks, scalticks=0:2:10)