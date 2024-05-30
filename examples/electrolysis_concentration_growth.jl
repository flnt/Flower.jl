using Revise
using Flower
using Printf
using PrettyTables
using Interpolations
# using PyPlot

using PyCall
@pyimport matplotlib.animation as anim
using PyPlot


pygui(false) #do not show figures

#From channel.jl and heat_convection.jl
#Khalighi 2023: "Hydrogen bubble growth in alkaline water electrolysis: An immersed boundary simulation study"

plot_Makie=false
if plot_Makie
    fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
    set_theme!(fontsize_theme)
end

prefix="/local/home/pr277828/flower/"

folder="electrolysis_growth"

prefix *= "/"*folder*"/"

isdir(prefix) || mkdir(prefix)



L0 = 1e-4 
n = 64
max_iter=1000

# n=96

# n=128

n=256

max_iter=100000

max_iter=100

# max_iter=1
max_iter=2
save_every=1

# save_every=10

# max_iter=1

# max_iter=1


# n=64
# max_iter=1

#n=128
#max_iter=2


# max_iter=5
# n = 10
# max_iter=1


xlabel = L"x \left(\mu m \right)"
ylabel = L"y \left(\mu m \right)"

xscale = 1e-6
yscale = xscale

# xticks = 0:20:num.L0/xscale
# yticks = 0:20:num.L0/yscale

xticks = 0:20:100
yticks = 0:20:100

velscale = 1e-4 

concentrationcontour=true

concentrationcontour=false

plot_levelset=true

# cmap = plt.cm.bone
cmap = plt.cm.viridis

adapt_timestep_mode=2


h0 = 0.25*L0 #TODO h0

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)

function fmax(x,v_inlet_max,L0)
    return 4*v_inlet_max*x/L0*(1-x/L0)
end

function favg(x,v_inlet_moy,L0)
    return 6*v_inlet_moy*x/L0*(1-x/L0)
end

function zerovel(x,v_inlet_moy,L0)
    return 0.0
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

rho1=1258.0 #liquid
#TODO need for 80°C, check with other study
rho2=0.069575 #"0.7016E-01" in \citet{cohnTABLETHERMODYNAMICPROPERTIES1963} H2
#Linear interpolation between 350 and 360
# 350	360	353		B	0.13841	353	350	360
# 7.02E-02	6.82E-02		-0.00194999999999999	A	-0.000194999999999999	0.069575	0.07016	0.06821


# radius=2.5e-5 
# radius=1.25e-5 
radius=3.0e-6 


ref_thickness_2d = 4.0 / 3.0 *radius 

mode_2d = 1 #use equivalent cylinder


# xcoord = 0.0
# ycoord = 0.0
# xcoord = radius + 1e-6
# ycoord=5e-5
xcoord = radius + 2e-6 #or 3e-6? not written in the article
ycoord=5e-5

xcoord = 5e-5
ycoord=5e-5

xcoord= 1e-5

xcoord = -xcoord
ycoord = -ycoord

mu=6.7e-7
mu1=mu
mu2=mu #TODO
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
CFL= 0.5 #0.01 #0.5
Re=1.0
TEND=7.3#s

# elec_cond=1 #TODO

elec_cond=2*Faraday^2*c0_KOH*DKOH/(Ru*temperature0)


#print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")




v_inlet=6.7e-4
Re=rho1*v_inlet*L0/mu
printstyled(color=:green, @sprintf "\n Re : %.2e %.2e %.2e %.2e\n" Re rho1/mu1 rho1 mu1)

# save_every=max_iter
# save_every=1

# rho1=1.0
# Re=100.0
# mu1=rho1/Re
# mu=mu1

# Re=1.0
# rho1=1.0
# rho2=1.0
# mu1=1.0
# mu2=1.0
# mu=1.0

Re=rho1/mu1

printstyled(color=:green, @sprintf "\n Re : %.2e %.2e %.2e %.2e\n" Re rho1/mu1 rho1 mu1)





#Concentration
diffusion_coeff=[DH2, DKOH, DH2O]
concentration0=[0.16, 6700, 49000]
nb_transported_scalars=3

nb_saved_scalars=2


if length(concentration0)!=nb_transported_scalars
    print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
    @error ("nb_transported_scalars")
end

if length(diffusion_coeff)!=nb_transported_scalars
    print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)
    @error ("nb_transported_scalars")
end

# print(@sprintf "Re = %.2e\n" Re)




print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)

# pretty_table(concentration0'; header = ["cH2", "cKOH", "cH2O"])
# pretty_table(diffusion_coeff'; header = ["DH2", "DKOH", "DH2O"])
# pretty_table(vcat(hcat("D",diffusion_coeff'),hcat("c",concentration0')); header = ["","H2", "KOH", "H2O"])
# hl = Highlighter((d,i,j)->d[i,j][1]*d[i,j][2] < 0, crayon"red")

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold blue")

pretty_table(vcat(hcat("Diffusion coef",diffusion_coeff'),hcat("Concentration",concentration0')); 
header = ["","H2", "KOH", "H2O"], highlighters=hl)

num = Numerical(
    CFL = CFL,
    Re = Re,
    TEND=TEND,
    x = x,
    y = y,
    shifted = xcoord,
    shifted_y = ycoord,
    case = "Cylinder",#"Cylinder", #"Sphere", #"Nothing",
    R = radius,
    max_iterations = max_iter,
    save_every = save_every,#10,#save_every,#10,
    ϵ = 0.05, 
    nb_transported_scalars=nb_transported_scalars,
    nb_saved_scalars=nb_saved_scalars,
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
    mu1=mu,
    mu2=mu,
    rho1=rho1,
    rho2=rho2,
    u_inf = 0.0,
    pres0=pres0,
    sigma=sigma,
    ref_thickness_2d = ref_thickness_2d,
    )
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
# gp.LS[1].u .= 1.0
# figname0="no_intfc"

@unpack x, nx, ny, ind = gv

vPoiseuille = zeros(gv)
vPoiseuille=favg.(x,v_inlet,L0)
# # vPoiseuille=zerovel.(x,v_inlet,L0)


phL.v .=vPoiseuille
# phL.v .= zeros(gv)
# phL.v .= 0.0


phL.u .= 0.0
phL.p .= 0.0




su = sqrt.((gv.x .+ xcoord).^2 .+ (gv.y .+ ycoord).^2)
R1 = radius + 3.0*num.Δ

bl = 4.0
for II in gv.ind.all_indices
    if su[II] <= R1
        phL.v[II] = 0.0
    # elseif su[II] > R1
    #     uL[II] = tanh(bl*(su[II]-R1))
    end
end


vecb_L(phL.uD, gu) .= 0.0
vecb_B(phL.uD, gu) .= 0.0
vecb_R(phL.uD, gu) .= 0.0
vecb_T(phL.uD, gu) .= 0.0

vecb_L(phL.uD, gu) .= 0.0 #TODO
vecb_B(phL.uD, gu) .= 0.0
vecb_R(phL.uD, gu) .= 0.0
vecb_T(phL.uD, gu) .= 0.0


printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL CFL*L0/n/v_inlet)


xscale = 1e-6
yscale = xscale

x_array=gp.x[1,:]/xscale
y_array=gp.y[:,1]/yscale

us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

# fig, ax = plt.subplots()
# q = ax.quiver(x_array,y_array,us,vs)
# # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
# #              label='Quiver key, length = 10', labelpos='E')

# # plt.show()
# plt.axis("equal")

# plt.savefig(prefix*"vector0.pdf")




xu=gu.x[1,:]/xscale
yu=gu.y[:,1]/yscale

# u_array = phL.u 
u_array = phL.u ./velscale


fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xu,yu,u_array, 10, cmap=cmap)

# Note that in the following, we explicitly pass in a subset of the contour
# levels used for the filled contours.  Alternatively, we could pass in
# additional levels to provide extra resolution, or leave out the *levels*
# keyword argument to use all of the original levels.

# CS2 = ax2.contour(CS, 
# # levels=CS.levels[::2], 
# # levels=
# colors="r")

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("u")
# Add the contour line levels to the colorbar
# cbar.add_lines(CS2)


# cbar.formatter.set_powerlimits((0, 0))
# # to get 10^3 instead of 1e3
# # cbar.formatter.set_useMathText(True)
# cbar.formatter.set_useMathText(1)

cbar.ax.set_title(L"$10^{-4}$")

plt.axis("equal")


plt.savefig(prefix*"u0.pdf")
plt.close(fig1)

######################################################################################################

xv=gv.x[1,:]/xscale
yv=gv.y[:,1]/yscale

# v_array = phL.v 
v_array = phL.v ./velscale

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xv,yv,phL.v ./velscale, 10, cmap=cmap)

# Note that in the following, we explicitly pass in a subset of the contour
# levels used for the filled contours.  Alternatively, we could pass in
# additional levels to provide extra resolution, or leave out the *levels*
# keyword argument to use all of the original levels.

# CS2 = ax2.contour(CS, 
# # levels=CS.levels[::2], 
# # levels=
# colors="r")

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("v")
# Add the contour line levels to the colorbar
# cbar.add_lines(CS2)

# cbar.formatter.set_powerlimits((0, 0))
# # to get 10^3 instead of 1e3
# # cbar.formatter.set_useMathText(True)
# cbar.formatter.set_useMathText(1)

cbar.ax.set_title(L"$10^{-4}$")
plt.axis("equal")

plt.savefig(prefix*"v0.pdf")

plt.close(fig1)






vPoiseuilleb=favg.(gv.x[1,:],v_inlet,L0)

for iscal=1:nb_transported_scalars
    phL.trans_scal[:,:,iscal] .= concentration0[iscal]
end

phL.phi_ele .= phi_ele0
# phL.phi_eleD .= 0.0 #TODO  phL.phi_eleD .= 0.0 ?

printstyled(color=:green, @sprintf "\n Initialisation \n")

print_electrolysis_statistics(nb_transported_scalars,gp,phL)


printstyled(color=:green, @sprintf "\n TODO timestep CFL scal, and print \n")


@unpack τ,CFL,Δ,Re,θd=num
# print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ CFL CFL*Δ CFL*Δ^2*Re Re θd)
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

@time current_i=run_forward(
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
    # BC_int = [FreeSurface()], #[WallNoSlip()],
    BC_int = [WallNoSlip()],

    # BC_TL  = Boundaries(
    # left   = Dirichlet(val = temperature0),
    # right  = Dirichlet(val = temperature0),
    # bottom = Dirichlet(val = temperature0),
    # top    = Dirichlet(val = temperature0),
    # ),

    BC_trans_scal = (
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[1]),
        top    = Neumann(),
        left   = Neumann(val=-i_current/(2*Faraday*DH2)),
        right  = Neumann(),
        int    = Dirichlet(val = concentration0[1])), #H2
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=-i_current/(2*Faraday*DKOH)),
        right  = Neumann(),
        int    = Neumann(val=0.0)), #KOH
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=i_current/(Faraday*DH2O)),
        right  = Neumann(),
        int    = Neumann(val=0.0)) #H2O
    ),

    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=-i_current/elec_cond),
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
    ),


    time_scheme = FE,#CN, #or FE?
    electrolysis = true,

    navier_stokes = true,
    ns_advection=true, #false
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    electrolysis_convection = true,  
    electrolysis_liquid_phase = true,
    electrolysis_phase_change = true,
    electrolysis_phase_change_case = "Khalighi",
    electrolysis_reaction = "Butler_no_concentration",
    # electrolysis_reaction = "nothing",
    adapt_timestep_mode = adapt_timestep_mode,#1,
    non_dimensionalize=0,
    mode_2d = mode_2d,

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

printstyled(color=:green, @sprintf "\n max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(phL.u)) maximum(abs.(phL.v)))
# printstyled(color=:green, @sprintf "\n eps : %.2e \n" eps(0.01))


if plot_Makie

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
end

######################################################################################################
# Matplotlib
######################################################################################################
last_it=10

x_array=gp.x[1,:]/xscale
y_array=gp.y[:,1]/yscale

size_frame=size(fwd.trans_scal[:,:,:,1],1)

size_frame=current_i

function strtitlefunc(isnap)
    # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[i+1] fwd.radius[i+1]
    strtitle = @sprintf "t %.2e (ms) radius %.2e (mm)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6
    return strtitle
end

# phi_array=reshape(vec1(fwdL.phi_eleD[end,:], gp), gp)'

phi_array=phL.phi_ele #do not transpose since python row major

# print("test",phL.phi_ele[1,:])

# print("test",phL.phi_ele[:,1])

# Eus, Evs = interpolate_regular_grid(gp,fwdL,fwdL.Eu, fwdL.Ev)
Eus,Evs = interpolate_grid_liquid(gp,gu,gv,phL.Eu, phL.Ev)


#https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_demo.html

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(x_array,y_array,phi_array, 10, cmap=cmap)

# Note that in the following, we explicitly pass in a subset of the contour
# levels used for the filled contours.  Alternatively, we could pass in
# additional levels to provide extra resolution, or leave out the *levels*
# keyword argument to use all of the original levels.

CS2 = ax2.contour(CS, 
# levels=CS.levels[::2], 
# levels=
colors="r")

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("Electrical potential")
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

plt.streamplot(x_array,y_array, -Eus,-Evs,color="w")
# Eus[last_it,:,:], Evs[last_it,:,:])#, color=(.75,.90,.93)) #do no transpose, python row major
plt.axis("equal")

plt.savefig(prefix*"streamlines.pdf")
plt.close(fig1)


######################################################################################################

us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

# fig, ax = plt.subplots()
# q = ax.quiver(x_array,y_array,us,vs)
# # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
# #              label='Quiver key, length = 10', labelpos='E')
# plt.axis("equal")

# # plt.show()

# plt.savefig(prefix*"vector.pdf")
# plt.close(fig)


# print("\n test u ", vecb_L(phL.uD, gu))
# print("\n test u ", phL.u[1,:])
# print("\n test u ", phL.u[:,1])



# vecb_L(phL.uD, gu) .=0.0
# vecb_R(phL.uD, gu) .=0.0
# vecb_T(phL.uD, gu) .=0.0
# vecb_B(phL.uD, gu) .=0.0
# vec1(phL.uD, gu) .=0.0


# phL.u .= reshape(vec1(phL.uD,gu), gu)
# phL.v .= reshape(vec1(phL.vD,gv), gv)

# us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

# fig, ax = plt.subplots()
# q = ax.quiver(x_array,y_array,us,vs)
# # ax.quiverkey(q, X=0.3, Y=1.1, U=10,
# #              label='Quiver key, length = 10', labelpos='E')

# plt.axis("equal")

# plt.show()

# plt.savefig(prefix*"vector_test1.pdf")

# # print("\n test u ", vecb_L(phL.uD, gu))

######################################################################################################


xu=gu.x[1,:]/xscale
yu=gu.y[:,1]/yscale

# u_array = phL.u 
u_array = phL.u ./velscale


fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xu,yu,u_array, 10, cmap=cmap)

# Note that in the following, we explicitly pass in a subset of the contour
# levels used for the filled contours.  Alternatively, we could pass in
# additional levels to provide extra resolution, or leave out the *levels*
# keyword argument to use all of the original levels.

# CS2 = ax2.contour(CS, 
# # levels=CS.levels[::2], 
# # levels=
# colors="r")

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("u")
# Add the contour line levels to the colorbar
# cbar.add_lines(CS2)


# cbar.formatter.set_powerlimits((0, 0))
# # to get 10^3 instead of 1e3
# # cbar.formatter.set_useMathText(True)
# cbar.formatter.set_useMathText(1)

cbar.ax.set_title(L"$10^{-4}$")

plt.axis("equal")


plt.savefig(prefix*"u.pdf")
plt.close(fig1)

######################################################################################################

xv=gv.x[1,:]/xscale
yv=gv.y[:,1]/yscale

# v_array = phL.v 
v_array = phL.v ./velscale

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xv,yv,phL.v ./velscale, 10, cmap=cmap)

# Note that in the following, we explicitly pass in a subset of the contour
# levels used for the filled contours.  Alternatively, we could pass in
# additional levels to provide extra resolution, or leave out the *levels*
# keyword argument to use all of the original levels.

# CS2 = ax2.contour(CS, 
# # levels=CS.levels[::2], 
# # levels=
# colors="r")

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("v")
# Add the contour line levels to the colorbar
# cbar.add_lines(CS2)

# cbar.formatter.set_powerlimits((0, 0))
# # to get 10^3 instead of 1e3
# # cbar.formatter.set_useMathText(True)
# cbar.formatter.set_useMathText(1)

cbar.ax.set_title(L"$10^{-4}$")
plt.axis("equal")

plt.savefig(prefix*"v.pdf")

plt.close(fig1)


######################################################################################################

fig1, ax2 = plt.subplots(layout="constrained")
# CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
CS = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, 10, cmap=cmap)

# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("concentration")
# Add the contour line levels to the colorbar
if concentrationcontour
    CS2 = ax2.contour(CS, 
    # levels=CS.levels[::2], 
    # levels=
    colors="r")
    cbar.add_lines(CS2)
end
if plot_levelset
    # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
    # CS2 = ax2.contour(CSlvl, 
    # # levels=CS.levels[::2], 
    # # levels=
    # colors="r")
    # cbar.add_lines(CS2)
    CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
end

# if plot_levelset
#     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

#     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
# end

plt.axis("equal")

plt.savefig(prefix*"H2.pdf")
plt.close(fig1)

#TODO python vs julia storage order
#OK if table sent
#if plot i j : careful

for isnap in 1:size_frame

    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
    # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[isnap,:,:,1] .-c0_H2)./c0_H2, 10, cmap=cmap)
    CS = ax2.contourf(x_array,y_array, fwd.trans_scal[isnap,:,:,1], 
    10, 
    cmap=cmap)
    # print("nx ny ", nx,ny)
    nplot=5
    nplotx=nx ÷ nplot 
    nploty=ny ÷ nplot
    fontsize=5
    startmod=1
    ms=2
    mcolor="w"

    # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[isnap] fwd.radius[isnap]
    # strtitle = @sprintf "t %.2e (ms) radius %.2e (mm)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6
    strtitle=strtitlefunc(isnap)

    plt.title(strtitle)

    for i in 1:gp.nx
        for j in 1:gp.ny
            if (i%nplotx==startmod) && (j%nploty==startmod)

                # print("\nplot i ",i," j ",j)

                str=@sprintf "%.2e" fwd.trans_scal[isnap,j,i,1]
                # str=@sprintf "%.2e %.2e" x_array[i] y_array[j]
                # str=@sprintf "%.2e" x_array[i]

                # str=@sprintf "%.5i %.5i" i j 


                ax2.annotate(str,(x_array[i],y_array[j]),fontsize=fontsize)
                plt.scatter(x_array[i],y_array[j],s=ms,c=mcolor)

            end
        end
    end
    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel("concentration")
    # Add the contour line levels to the colorbar
    if concentrationcontour
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)
    end
    if plot_levelset
        # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
        # CS2 = ax2.contour(CSlvl, 
        # # levels=CS.levels[::2], 
        # # levels=
        # colors="r")
        # cbar.add_lines(CS2)
        CSlvl = ax2.contour(x_array,y_array, fwd.u[1,isnap,:,:], [0.0],colors="r")
    end

    # if plot_levelset
    #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

    #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
    # end

    plt.axis("equal")

    plt.savefig(prefix*"H2_"*string(isnap)*".pdf")
    plt.close(fig1)
end




for isnap in 1:size_frame

    fig1, ax2 = plt.subplots(layout="constrained")
    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)
    # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[isnap,:,:,1] .-c0_H2)./c0_H2, 10, cmap=cmap)
    CS = ax2.contourf(x_array,y_array, fwd.mass_flux[isnap,:,:], 10, cmap=cmap)
    # print("nx ny ", nx,ny)
    nplot=5
    nplotx=nx ÷ nplot 
    nploty=ny ÷ nplot
    fontsize=5
    startmod=1
    ms=2
    mcolor="w"

    # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[isnap] fwd.radius[isnap]
    # strtitle = @sprintf "t %.2e (ms) radius %.2e (mm)" fwd.t[isnap]*1e3 fwd.radius[isnap]*1e6

    strtitle=strtitlefunc(isnap)

    plt.title(strtitle)

    for i in 1:gp.nx
        for j in 1:gp.ny
            if (i%nplotx==startmod) && (j%nploty==startmod)

                # print("\nplot i ",i," j ",j)

                str=@sprintf "%.2e" fwd.mass_flux[isnap,j,i]
                # str=@sprintf "%.2e %.2e" x_array[i] y_array[j]
                # str=@sprintf "%.2e" x_array[i]

                # str=@sprintf "%.5i %.5i" i j 


                ax2.annotate(str,(x_array[i],y_array[j]),fontsize=fontsize)
                plt.scatter(x_array[i],y_array[j],s=ms,c=mcolor)

            end
        end
    end
    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig1.colorbar(CS)
    cbar.ax.set_ylabel("flux")
    # Add the contour line levels to the colorbar
    if concentrationcontour
        CS2 = ax2.contour(CS, 
        # levels=CS.levels[::2], 
        # levels=
        colors="r")
        cbar.add_lines(CS2)
    end
    # if plot_levelset
    #     # CSlvl = ax2.contourf(x_array,y_array,(phL.trans_scal[:,:,1] .-c0_H2)./c0_H2, levels=0.0, cmap=cmap)
    #     # CS2 = ax2.contour(CSlvl, 
    #     # # levels=CS.levels[::2], 
    #     # # levels=
    #     # colors="r")
    #     # cbar.add_lines(CS2)
    #     CSlvl = ax2.contour(x_array,y_array, fwd.u[1,isnap,:,:], [0.0],colors="r")
    # end

    # if plot_levelset
    #     gp.LS[1].u .= sqrt.((gp.x .+ xcoord).^ 2 + (gp.y .+ ycoord) .^ 2) - (radius) * ones(gp);

    #     CSlvl = ax2.contour(x_array,y_array, gp.LS[1].u, [0.0],colors="r")
    # end

    plt.axis("equal")

    plt.savefig(prefix*"H2flux_"*string(isnap)*".pdf")
    plt.close(fig1)
end
######################################################################################################
fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,2] .-c0_H2O)./c0_H2O,0.0), 10, cmap=cmap)


# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("concentration")
# Add the contour line levels to the colorbar

if concentrationcontour
    CS2 = ax2.contour(CS, 
    # levels=CS.levels[::2], 
    # levels=
    colors="r")
    cbar.add_lines(CS2)
end

plt.axis("equal")

plt.savefig(prefix*"KOH.pdf")
plt.close(fig1)

######################################################################################################

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,2] .-c0_KOH)./c0_KOH,0.0), 10, cmap=cmap)


# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("concentration")
# Add the contour line levels to the colorbar

if concentrationcontour
    CS2 = ax2.contour(CS, 
    # levels=CS.levels[::2], 
    # levels=
    colors="r")
    cbar.add_lines(CS2)
end

plt.axis("equal")

plt.savefig(prefix*"H2O.pdf")
plt.close(fig1)

######################################################################################################


fig1, ax2 = plt.subplots(layout="constrained")

# print("t",fwd.t)
# print("radius",fwd.radius.*1.e6)
# print("current_i", current_i)
# print("radius ",fwd.radius[1:current_i+1])
# print("\nradius ",fwd.radius)

plt.plot(fwd.t[1:current_i],fwd.radius[1:current_i].*1.e6)
# ax2.set_title("Title")
ax2.set_xlabel(L"$t (s)$")
# ax2.set_ylabel(L"$R (m)$")
ax2.set_ylabel(L"$R (\mu m)$")

ax2.set_xscale("log")
ax2.set_yscale("log")

plt.axis("equal")

plt.savefig(prefix*"R.pdf")
plt.close(fig1)


#H2O
isca=3

fig1, ax2 = plt.subplots(layout="constrained")

CS = ax2.contourf(x_array,y_array,fwd.trans_scal[1,:,:,isca], 
levels=10, 
# levels=range(0,1400,step=200),
cmap=cmap)

 # Make a colorbar for the ContourSet returned by the contourf call.
 cbar = fig1.colorbar(CS)
 cbar.ax.set_ylabel("concentration")

 cbarax=cbar.ax

function make_frame(i)
    # ax1.clear()
    ax2.clear()
    # ax1.imshow(A[:,:,i+1, 1])

    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

    CS = ax2.contourf(x_array,y_array,fwd.trans_scal[i+1,:,:,isca], 
    levels=10, 
    # levels=range(0,1400,step=200),
    cmap=cmap)


    strtitle=strtitlefunc(i+1)
    plt.title(strtitle)

    # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
    # # levels=10, 
    # levels=range(0,1400,step=200),
    # cmap=cmap)

    plt.axis("equal")


    # ax2.imshow(A[:,:,i+1])

    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
     cbar = fig1.colorbar(CS,cax=cbarax)
    #  cbar.ax.set_ylabel("concentration")

    #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

    if (i==0) #(i+1==1) python starts at 0
        # # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel("concentration")
        
        # Add the contour line levels to the colorbar
        if concentrationcontour
            CS2 = ax2.contour(CS, 
            # levels=CS.levels[::2], 
            # levels=
            colors="r")
            cbar.add_lines(CS2)
        end
    end
 
    if plot_levelset
        CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
    end


    # plt.savefig(prefix*"H2.pdf")
    # plt.close(fig1)
end

myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

# myanim.save(prefix*"test.gif")
myanim.save(prefix*"H2O.mp4")

plt.close("all")

######################################################################################################
# function pymovie(fwdL,iscal,name)
    
#     c0ref = concentration0[iscal]

#     fig1, ax2 = plt.subplots(layout="constrained")

#     lim2=200
#     lim1=-lim2
#     step1=10

#     CS = ax2.contourf(x_array,y_array,max.((fwdL.trans_scal[1,:,:,isca].-c0ref)./c0ref,0.0), 
#     #levels=10, 
#     levels=range(lim1,lim2,step=step1),
#     cmap=cmap)

#     # Make a colorbar for the ContourSet returned by the contourf call.
#     cbar = fig1.colorbar(CS)
#     cbar.ax.set_ylabel("concentration")

#     cbarax=cbar.ax

#     function make_frame_movie(i)
#         # ax1.clear()
#         ax2.clear()
#         # ax1.imshow(A[:,:,i+1, 1])

#         # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

#         CS = ax2.contourf(x_array,y_array,max.((fwdL.trans_scal[i+1,:,:,isca].-c0ref)./c0ref,0.0), 
#         #levels=10, 
#         levels=range(lim1,lim2,step=step1),
#         cmap=cmap)


#         strtitle=strtitlefunc(i+1)
#         plt.title(strtitle)

#         # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
#         # # levels=10, 
#         # levels=range(0,1400,step=200),
#         # cmap=cmap)

#         plt.axis("equal")


#         # ax2.imshow(A[:,:,i+1])

#         # ax2.set_title("Title")
#         ax2.set_xlabel(L"$x (\mu m)$")
#         ax2.set_ylabel(L"$y (\mu m)$")

#         # Make a colorbar for the ContourSet returned by the contourf call.
#         cbar = fig1.colorbar(CS,cax=cbarax)
#         #  cbar.ax.set_ylabel("concentration")

#         #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

#         if (i==0) #(i+1==1) python starts at 0
#             # # Make a colorbar for the ContourSet returned by the contourf call.
#             # cbar = fig1.colorbar(CS)
#             # cbar.ax.set_ylabel("concentration")
            
#             # Add the contour line levels to the colorbar
#             if concentrationcontour
#                 CS2 = ax2.contour(CS, 
#                 # levels=CS.levels[::2], 
#                 # levels=
#                 colors="r")
#                 cbar.add_lines(CS2)
#             end
#         end
    
#         #if plot_levelset
#         #   CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
#         # end


#         # plt.savefig(prefix*"H2.pdf")
#         # plt.close(fig1)
#     end

#     myanim = anim.FuncAnimation(fig1, make_frame_movie, frames=size_frame, interval=size_frame, blit=false)

#     # myanim.save(prefix*"test.gif")
#     myanim.save(prefix*name*".mp4")

#     plt.close("all")


# end
# ######################################################################################################

# pymovie(fwdL, iscal=2, name="KOH")

######################################################################################################

#H2O
isca=3

fig1, ax2 = plt.subplots(layout="constrained")

lim2=1
lim1=-lim2
step1=10

CS = ax2.contourf(x_array,y_array,max.((fwdL.trans_scal[1,:,:,isca].-c0_H2O)./c0_H2O,0.0), 
#levels=10, 
levels=range(lim1,lim2,step=step1),
cmap=cmap)

 # Make a colorbar for the ContourSet returned by the contourf call.
 cbar = fig1.colorbar(CS)
 cbar.ax.set_ylabel("concentration")

 cbarax=cbar.ax

function make_frame(i)
    # ax1.clear()
    ax2.clear()
    # ax1.imshow(A[:,:,i+1, 1])

    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

    CS = ax2.contourf(x_array,y_array,max.((fwdL.trans_scal[i+1,:,:,isca].-c0_H2O)./c0_H2O,0.0), 
    #levels=10, 
    levels=range(lim1,lim2,step=step1),
    cmap=cmap)


    strtitle=strtitlefunc(i+1)
    plt.title(strtitle)

    # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
    # # levels=10, 
    # levels=range(0,1400,step=200),
    # cmap=cmap)

    plt.axis("equal")


    # ax2.imshow(A[:,:,i+1])

    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    # Make a colorbar for the ContourSet returned by the contourf call.
     cbar = fig1.colorbar(CS,cax=cbarax)
    #  cbar.ax.set_ylabel("concentration")

    #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

    if (i==0) #(i+1==1) python starts at 0
        # # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel("concentration")
        
        # Add the contour line levels to the colorbar
        if concentrationcontour
            CS2 = ax2.contour(CS, 
            # levels=CS.levels[::2], 
            # levels=
            colors="r")
            cbar.add_lines(CS2)
        end
    end
 
    #if plot_levelset
    #   CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
    # end


    # plt.savefig(prefix*"H2.pdf")
    # plt.close(fig1)
end

myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

# myanim.save(prefix*"test.gif")
myanim.save(prefix*"H2O_liq_norm.mp4")

plt.close("all")

######################################################################################################

# A = (fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2

fig1, ax2 = plt.subplots(layout="constrained")

CS = ax2.contourf(x_array,y_array,max.((fwd.trans_scal[1,:,:,1] .-c0_H2)./c0_H2,0.0), 
# levels=10, 
levels=range(0,1400,step=200),
cmap=cmap)

 # Make a colorbar for the ContourSet returned by the contourf call.
 cbar = fig1.colorbar(CS)
 cbar.ax.set_ylabel("concentration")

function make_frame(i)
    # ax1.clear()
    ax2.clear()
    # ax1.imshow(A[:,:,i+1, 1])

    # CS = ax2.contourf(x_array,y_array,max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), 10, cmap=cmap)

    CS = ax2.contourf(x_array,y_array,max.((fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2,0.0), 
    # levels=10, 
    levels=range(0,1400,step=200),
    cmap=cmap)


    strtitle=strtitlefunc(i+1)
    plt.title(strtitle)

    # CS = ax2.contourf(x_array,y_array,(fwd.trans_scal[i+1,:,:,1] .-c0_H2)./c0_H2, 
    # # levels=10, 
    # levels=range(0,1400,step=200),
    # cmap=cmap)

    plt.axis("equal")


    # ax2.imshow(A[:,:,i+1])

    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

     # Make a colorbar for the ContourSet returned by the contourf call.
    #  cbar = fig1.colorbar(CS)
    #  cbar.ax.set_ylabel("concentration")

    #https://stackoverflow.com/questions/5180518/duplicated-colorbars-when-creating-an-animation

    if (i==0) #(i+1==1) python starts at 0
        # # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel("concentration")
        
        # Add the contour line levels to the colorbar
        if concentrationcontour
            CS2 = ax2.contour(CS, 
            # levels=CS.levels[::2], 
            # levels=
            colors="r")
            cbar.add_lines(CS2)
        end
    end
 
    if plot_levelset
        CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
    end


    # plt.savefig(prefix*"H2.pdf")
    # plt.close(fig1)
end

myanim = anim.FuncAnimation(fig1, make_frame, frames=size_frame, interval=size_frame, blit=false)

# myanim.save(prefix*"test.gif")
myanim.save(prefix*"H2_norm.mp4")

plt.close("all")

######################################################################################################

######################################################################################################

# size_frame=size(fwd.trans_scal[:,:,:,1],1)

# fig1, ax2 = plt.subplots(layout="constrained")

# CS = ax2.contourf(x_array,y_array,fwd.trans_scal[1,:,:,1], 
# levels=10, 
# # levels=range(0,1400,step=200),
# cmap=cmap)


#  # Make a colorbar for the ContourSet returned by the contourf call.
#  cbar = fig1.colorbar(CS)
#  cbar.ax.set_ylabel("concentration")

# function make_frame_2(i)
#     # ax1.clear()
#     ax2.clear()
#     # ax1.imshow(A[:,:,i+1, 1])

#     CS = ax2.contourf(x_array,y_array,fwd.trans_scal[i+1,:,:,1], 
#     levels=10, 
#     # levels=range(0,1400,step=200),
#     cmap=cmap)

#     plt.axis("equal")


#     # ax2.imshow(A[:,:,i+1])

#     # ax2.set_title("Title")
#     ax2.set_xlabel(L"$x (\mu m)$")
#     ax2.set_ylabel(L"$y (\mu m)$")

#     if (i==0) #(i+1==1) python starts at 0
#         # # Make a colorbar for the ContourSet returned by the contourf call.
#         # cbar = fig1.colorbar(CS)
#         # cbar.ax.set_ylabel("concentration")
        
#         # Add the contour line levels to the colorbar
#         if concentrationcontour
#             CS2 = ax2.contour(CS, 
#             # levels=CS.levels[::2], 
#             # levels=
#             colors="r")
#             cbar.add_lines(CS2)
#         end
#     end
 
#     if plot_levelset
#         CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
#     end


#     # plt.savefig(prefix*"H2.pdf")
#     # plt.close(fig1)
# end

# myanim = anim.FuncAnimation(fig1, make_frame_2, frames=size_frame, interval=size_frame, blit=false)

# # myanim.save(prefix*"test.gif")
# myanim.save(prefix*"H2.mp4")

# plt.close("all")

######################################################################################################

######################################################################################################

var = fwd.mass_flux
# size_frame=size(var,1)

fig1, ax2 = plt.subplots(layout="constrained")
lvl=1e-5 #0.1
step = 2*lvl/10
levels=range(-lvl,lvl,step=step)

levels=10

CS = ax2.contourf(x_array,y_array,var[1,:,:], 
# levels=10, 
levels=levels,
cmap=cmap)


 # Make a colorbar for the ContourSet returned by the contourf call.
 cbar = fig1.colorbar(CS)
 cbar.ax.set_ylabel("mass_flux")

function make_frame_3(i)
    # ax1.clear()
    ax2.clear()
    # ax1.imshow(A[:,:,i+1, 1])

    CS = ax2.contourf(x_array,y_array,var[i+1,:,:], 
    levels=levels, 
    cmap=cmap)

    # cbar = fig1.colorbar(CS)
    # cbar.ax.set_ylabel("concentration")

    plt.axis("equal")

    # strtitle = @sprintf "t %.2e radius %.2e" fwd.t[i+1] fwd.radius[i+1]
    strtitle=strtitlefunc(i+1)
    plt.title(strtitle)


    # ax2.imshow(A[:,:,i+1])

    # ax2.set_title("Title")
    ax2.set_xlabel(L"$x (\mu m)$")
    ax2.set_ylabel(L"$y (\mu m)$")

    if (i==0) #(i+1==1) python starts at 0
        # # Make a colorbar for the ContourSet returned by the contourf call.
        # cbar = fig1.colorbar(CS)
        # cbar.ax.set_ylabel("concentration")
        
        # Add the contour line levels to the colorbar
        if concentrationcontour
            CS2 = ax2.contour(CS, 
            # levels=CS.levels[::2], 
            # levels=
            colors="r")
            cbar.add_lines(CS2)
        end
    end
 
    # if plot_levelset
    #     CSlvl = ax2.contour(x_array,y_array, fwd.u[1,i+1,:,:], [0.0],colors="r")
    # end


    # plt.savefig(prefix*"H2.pdf")
    # plt.close(fig1)
end

myanim = anim.FuncAnimation(fig1, make_frame_3, frames=size_frame, interval=size_frame, blit=false)

# myanim.save(prefix*"test.gif")
myanim.save(prefix*"H2flux.mp4")

plt.close("all")

######################################################################################################

# https://genkuroki.github.io/documents/Jupyter/20170624%20Examples%20of%20animations%20in%20Julia%20by%20PyPlot%20and%20matplotlib.animation.html

######################################################################################################

# using PyCall
# @pyimport matplotlib.animation as anim
# using PyPlot

# function showmp4(filename)
#     open(filename) do f
#         base64_video = base64encode(f)
#         display("text/html", """<video controls src="data:video/x-m4v;base64,$base64_video">""")
#     end
# end

# fig = figure()
# ax = axes()

# x = [-1.5:0.05:1.0;]
# y = [-1.0:0.05:1.0;]
# z = x' .+ im*y
# w0 = sqrt.(z.+1)./z

# # plotw has period 1:
# # plotw(t+1) == plotw(1)
# function plotw(t)
#     clf()
#     title("\$t = $t\$")
#     xlim(-1.4,0.9)
#     ylim(-0.9,0.9)
#     w = exp(pi*im*t)*w0
#     streamplot(x,y,real(w),-imag(w))
#     axes()[:set_aspect]("equal")
#     # quick dirty hack to prevent the following error:
#     # TypeError("'StreamplotSet' object is not iterable",)
#     plot()
# end

# n=10 # frame/loop
# l=2 # loop
# interval=200 # milli seconds

# # k=0,1,...,frames-1
# function animate(k)
#     plotw(k/n)
# end

# function init()
#     plotw(0)
# end

# withfig(fig) do
#     global myanim = anim.FuncAnimation(fig, animate, frames=l*n+1, init_func=init, interval=interval, blit=true)
#     myanim[:save]("test4.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
# end

# showmp4("test4.mp4")

# A = randn(20,20,20,2)
# fig, axesplt = PyPlot.subplots(nrows=1, ncols=2, figsize=(7, 2.5))
# ax1, ax2 = axesplt

# function make_frame(i)
#     ax1.clear()
#     ax2.clear()
#     ax1.imshow(A[:,:,i+1, 1])
#     ax2.imshow(A[:,:,i+1, 2])
# end

# myanim = anim.FuncAnimation(fig, make_frame, frames=size(A,3), interval=20, blit=false)

# myanim.save(prefix*"test.gif")
# myanim.save(prefix*"test.mp4")
