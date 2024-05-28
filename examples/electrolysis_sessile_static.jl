using Revise
using Flower
using PrettyTables
# using Printf
# using Interpolations
# using PyCall
# using PyPlot

prefix="/local/home/pr277828/flower/"

folder="electrolysis_sessile_static"

prefix *= "/"*folder*"/"

isdir(prefix) || mkdir(prefix)

pygui(false) #do not show figures

PyPlot.rc("text", usetex=true)
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.latex.preamble"] = raw"\usepackage{siunitx}"


plot_Makie=false
if plot_Makie
    fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
    set_theme!(fontsize_theme)
end

####################################################################################################
# Sessile from sessile.jl
####################################################################################################
Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)
####################################################################################################

L0 = 1e-4 
n = 64 
# n = 128
max_iter=100
save_every=1

xlabel = L"x \left(\mu m \right)"
ylabel = L"y \left(\mu m \right)"

xscale = 1e-6
yscale = xscale

xticks = 0:20:100
yticks = 0:20:100

velscale = 1e-4 

concentrationcontour=true
concentrationcontour=false

plot_levelset=true

cmap = plt.cm.viridis

adapt_timestep_mode=2

####################################################################################################

# h0=radius ?
# h0 = 0.25*L0 #TODO h0

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)


####################################################################################################

rho1=1258.0 #liquid
#TODO need for 80°C, check with other study
rho2=0.069575 #"0.7016E-01" in \citet{cohnTABLETHERMODYNAMICPROPERTIES1963} H2
#Linear interpolation between 350 and 360
# 350	360	353		B	0.13841	353	350	360
# 7.02E-02	6.82E-02		-0.00194999999999999	A	-0.000194999999999999	0.069575	0.07016	0.06821


# radius=2.5e-5 
# radius=1.25e-5 
radius = 3.0e-6 
radius = 6.0e-6 

h0 = radius


ref_thickness_2d = 4.0 / 3.0 *radius 

# mode_2d = 1 #use equivalent cylinder
mode_2d = 2 #mol/meter


xcoord = 0.0
ycoord = L0/2.0

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


print(@sprintf "TODO elec cond and boundary conditions need to be updated for potential\n")

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

print(@sprintf "nb_transported_scalars = %5i\n" nb_transported_scalars)

# pretty_table(concentration0'; header = ["cH2", "cKOH", "cH2O"])
# pretty_table(diffusion_coeff'; header = ["DH2", "DKOH", "DH2O"])
# pretty_table(vcat(hcat("D",diffusion_coeff'),hcat("c",concentration0')); header = ["","H2", "KOH", "H2O"])
# hl = Highlighter((d,i,j)->d[i,j][1]*d[i,j][2] < 0, crayon"red")

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold blue")

pretty_table(vcat(hcat("Diffusion coef",diffusion_coeff'),hcat("Concentration",concentration0')); 
header = ["","H2", "KOH", "H2O"], highlighters=hl)

# printstyled(color=:green, @sprintf "\n nmol : \n" concentration0[1]*4.0/3.0*pi*radius^3  )

current_radius = radius

p_liq= pres0 #+ mean(veci(phL.pD,grid,2)) #TODO here one bubble
p_g=p_liq + 2 * sigma / current_radius

c0test = p_g / (temperature0 * Ru) 


printstyled(color=:green, @sprintf "\n c0test: %.2e \n" c0test)

# printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(temperature0*num.Ru))

printstyled(color=:green, @sprintf "\n Species diffusion timescales: %.2e %.2e %.2e \n" (radius^2)/DH2 (radius^2)/DKOH (radius^2)/DH2O )



####################################################################################################


# h0 = 0.5
# L0x = 3.0
# L0y = 2.0
# n = 96

# x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
# y = collect(LinRange(-L0y / 2, 0, n ÷ 3 + 1))



####################################################################################################
θe= 90
θe= 145

if θe < 40
    max_its = 35000
    n_ext = 10
    CFL = 0.5
elseif θe < 100
    max_its = 15000
    n_ext = 10
    CFL = 0.5
else
    max_its = 5000
    n_ext = 10
    CFL = 0.5
end
####################################################################################################

max_its=10

# save_every = max_its÷100
save_every = 1
pres0 = 1e5

####################################################################################################
# not imposing the angle exactly at the boundary but displaced a cell because ghost cells are not used. 
# So the exact contact angle cannot be imposed
# TODO should at some point modify the levelset so it works as the other fields that we have in the code, 
# with the boundary values in addition to the bulk field. That way we could impose it exactly at the boundary
# needs a lot of work, not priority
_θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
thetaref = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
println("θe = $(_θe)")
####################################################################################################



num = Numerical(
    CFL = CFL,
    Re = Re,
    TEND=TEND,
    x = x,
    y = y,
    shifted = xcoord,
    shifted_y = ycoord,
    case = "Planar",
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
    v_inf = 0.0,
    pres0=pres0,
    σ=sigma,   
    reinit_every = 10,
    nb_reinit = 2,
    δreinit = 10.0,
    n_ext_cl = n_ext,
    NB = 24,
    )
    # ref_thickness_2d = ref_thickness_2d,



#Initialization
gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

####################################################################################################
# r = 0.5
# gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
# gp.LS[1].u .*= -1.0
####################################################################################################

gp.LS[1].u .= sqrt.((gp.x .- xcoord).^2 + (gp.y .- ycoord).^2) - radius * ones(gp)







@unpack x, nx, ny, ind = gv

vPoiseuille = zeros(gv)
vPoiseuille=Poiseuille_favg.(x,v_inlet,L0)

phL.v .=vPoiseuille
phL.u .= 0.0
phL.p .= pres0 #0.0
phL.T .= temperature0
phS.T .= temperature0

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

v_array = phL.v ./velscale

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xv,yv,phL.v ./velscale, 10, cmap=cmap)


# ax2.set_title("Title")
ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("v")

cbar.ax.set_title(L"$10^{-4}$")
plt.axis("equal")

plt.savefig(prefix*"v0.pdf")

plt.close(fig1)
######################################################################################################


vPoiseuilleb=Poiseuille_favg.(gv.x[1,:],v_inlet,L0)

for iscal=1:nb_transported_scalars
    phL.trans_scal[:,:,iscal] .= concentration0[iscal]
end

phL.phi_ele .= phi_ele0

printstyled(color=:green, @sprintf "\n Initialisation \n")

print_electrolysis_statistics(nb_transported_scalars,gp,phL)

printstyled(color=:green, @sprintf "\n TODO timestep CFL scal, and print \n")


@unpack τ,CFL,Δ,Re,θd=num
# print(@sprintf "dt %.2e %.2e %.2e %.2e %.2e %.2e\n" τ CFL CFL*Δ CFL*Δ^2*Re Re θd)
# τ=CFL*Δ/v_inlet
# num.τ=τ
# print(@sprintf "dt %.2e \n" τ)

#TODO pressure left and right BC not mentioned in the article Khalighi 2023

#TODO need to iterate more for potential since phiele=0 initially?

phi_ele=gv.x[1,:] .*0.0
# eta = phi_ele1 .-phi_ele
#TODO precision: number of digits
# i_current=i0*(exp(alphaa*Faraday*eta/(Ru*temperature0))-exp(-alphac*Faraday*eta/(Ru*temperature0)))
i_current=butler_volmer_no_concentration.(alphaa,alphac,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)

print(@sprintf "Butler-Volmer %.2e \n" i_current[1])



######################################################################################################
BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_cl(θe = _θe * π / 180),
    right = Neumann_inh()
)

@time current_i=run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left   = Navier_cl(λ = 1e-2), #Dirichlet(),
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
    # left = Dirichlet(val=pres0),
    # right = Dirichlet(val=pres0),

    BC_u = BC_u,

    BC_int = [FreeSurface()],

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
        int    = Dirichlet(val = concentration0[3])) #H2O
    ),

    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=-i_current/elec_cond),
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
    ),

    auto_reinit = true,
    save_length = true,
    time_scheme = FE,#CN, #or FE?
    electrolysis = true,
    navier_stokes = true,
    ns_advection=false,
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

printstyled(color=:green, @sprintf "\n max abs(u) : %.2e max abs(v)%.2e\n" maximum(abs.(phL.u)) maximum(abs.(phL.v)))


   



######################################################################################################

######################################################################################################
# V0 = 0.5 * π * 0.5^2
V0 = 0.5 * π * radius^2

Vf = volume(gp.LS[1].geoL)
Vratio = Vf / V0

mean_rad = 1 / abs(mean(gp.LS[1].κ[gp.LS[1].MIXED[5:end-5]]))
RR0_sim = mean_rad / 0.5
RR0_theo = RR0(θe * π / 180)

println("Vratio = $(Vratio)")
println("mean rad = $(mean_rad)")
println("RR0_sim = $(RR0_sim)")
println("RR0_theo = $(RR0_theo)")

suffix = "$(θe)deg_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
# suffix = "$(θe)deg_$(num.max_iterations)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
file = suffix*".jld2"
# save_field(prefix*file, num, gp, phL, fwdL, fwd)

tcks = -num.L0/2:0.5:num.L0
lim = (num.L0 + num.Δ) / 2
lim = 1.0
######################################################################################################


######################################################################################################
x_array=gp.x[1,:]/xscale
y_array=gp.y[:,1]/yscale

plot_levelset=true
isocontour=false#true
cmap = plt.cm.viridis

if isnothing(current_i)
    size_frame=size(fwdL.p,1)
else
    size_frame=current_i
end

######################################################################################################
# Streamlines
######################################################################################################
phi_array=phL.phi_ele #do not transpose since python row major
Eus,Evs = interpolate_grid_liquid(gp,gu,gv,phL.Eu, phL.Ev)

#https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_demo.html

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(x_array,y_array,phi_array, 10, cmap=cmap)

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
# us,vs = interpolate_grid_liquid(gp,gu,gv,phL.u,phL.v)

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
# u
######################################################################################################
xu=gu.x[1,:]/xscale
yu=gu.y[:,1]/yscale

# u_array = phL.u 
u_array = phL.u ./velscale

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xu,yu,u_array, 10, cmap=cmap)

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

######################################################################################################
# v
######################################################################################################
xv=gv.x[1,:]/xscale
yv=gv.y[:,1]/yscale

# v_array = phL.v 
v_array = phL.v ./velscale

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(xv,yv,phL.v ./velscale, 10, cmap=cmap)

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

iLS=1
rhs_LS = fzeros(gp)
#test contact angle
xp,yp,xp0,yp0=BC_LS_test!(gp, gp.LS[iLS].u, gp.LS[iLS].A, gp.LS[iLS].B, rhs_LS, BC_u)
printstyled(color=:green, @sprintf "\n θe : %.2e °\n" _θe*180.0/π)
printstyled(color=:green, @sprintf "\n θe : %.2e °\n" thetaref*180.0/π)

######################################################################################################

fig1, ax2 = plt.subplots(layout="constrained")
CS = ax2.contourf(x_array,y_array,gp.LS[1].u ./xscale, 10, cmap=cmap)

# print("\n xp",x_array)
# print("\n yp",y_array)

# print("\n xp",xp)
# print("\n yp",yp)
# print("\n xp0",xp0)
# print("\n yp0",yp0)

xp  = xp/xscale
yp  = yp/xscale
xp0 = xp0/xscale
yp0 = yp0/xscale



ms=2
mcolor="w"

plt.scatter(xp,yp,s=ms,c=mcolor)

mcolor="black"

plt.scatter(xp0,yp0,s=ms,c=mcolor)



ax2.set_xlabel(L"$x (\mu m)$")
ax2.set_ylabel(L"$y (\mu m)$")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
# cbar.ax.set_ylabel("v")

cbar.ax.set_title(raw"$( \unit{\um})$")

# cbar.ax.set_title(L"$10^{-6}$")

plt.axis("equal")

plt.savefig(prefix*"contact.pdf")

plt.close(fig1)
######################################################################################################

# plot_python_pdf(phL.p, "p",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure")

plot_python_pdf(phL.p, "p",prefix,plot_levelset,isocontour,0,range(pres0*0.9999,pres0*1.0001,length=10),cmap,x_array,y_array,gp,"pressure")

plot_python_pdf((phL.p.-pres0)./pres0, "pnorm",prefix,plot_levelset,isocontour,0,range(-1e-4,1e-4,length=10),cmap,x_array,y_array,gp,"pressure")


# python_movie_zoom(fwdL.p,"p",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
# size_frame,1,gp.nx,1,gp.ny,fwd)


python_movie_zoom((fwdL.p.-pres0)./pres0,"pnorm",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)

plot_python_pdf(max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), "H2",prefix,
plot_levelset,concentrationcontour,0,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration")

plot_python_pdf(max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), "H2lvl",prefix,
plot_levelset,concentrationcontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration")

python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2_norm",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),
cmap,x_array,y_array,gp,"concentration",size_frame,1,gp.nx,1,gp.ny,fwd)

######################################################################################################
#TODO bug when putting plotting instructions in function, range(start,end,step) no longer working as argument so give start,end,length instead



plot_python_pdf(max.((phL.trans_scal[:,:,2] .-c0_KOH)./c0_KOH,0.0), "KOH",prefix,
plot_levelset,concentrationcontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration")

plot_python_pdf(max.((phL.trans_scal[:,:,3] .-c0_H2O)./c0_H2O,0.0), "H2O",prefix,
plot_levelset,concentrationcontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration")

# plot_python_several_pdf(fwd.trans_scal[:,:,:,1],"H2",true,size_frame)
# plot_python_several_pdf(fwd.saved_scal[:,:,:,1],"H2massflux",true,size_frame)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2_norm",plot_levelset,size_frame,1,gp.nx,1,gp.ny,0,
range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,2] .-c0_KOH)./c0_KOH,0.0),"KOH_norm",
plot_levelset,size_frame,1,gp.nx,1,gp.ny,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,3] .-c0_H2O)./c0_H2O,0.0),"H2O_norm",
plot_levelset,size_frame,1,gp.nx,1,gp.ny,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(fwd.saved_scal[:,:,:,1],"flux1",false,size_frame,1,gp.nx,1,gp.ny,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(fwd.saved_scal[:,:,:,2],"flux2",false,size_frame,1,gp.nx,1,gp.ny,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
size_frame,1,gp.nx,1,gp.ny,fwd)

# python_movie_zoom(field,name,plot_levelset,size_frame,i0,i1,j0,j1,lmin,lmax,step)




######################################################################################################


######################################################################################################

fig1, ax2 = plt.subplots(layout="constrained")

# print("t",fwd.t)
# print("radius",fwd.radius.*1.e6)
# print("current_i", current_i)
# print("radius ",fwd.radius[1:current_i+1])
# print("\nradius ",fwd.radius)

plt.plot(fwd.t[1:size_frame],fwd.radius[1:size_frame].*1.e6)
# ax2.set_title("Title")
ax2.set_xlabel(L"$t (s)$")
# ax2.set_ylabel(L"$R (m)$")
ax2.set_ylabel(L"$R (\mu m)$")

ax2.set_xscale("log")
ax2.set_yscale("log")

plt.axis("equal")

plt.savefig(prefix*"R.pdf")
plt.close(fig1)

######################################################################################################


######################################################################################################





######################################################################################################
#TODO
# u,v,p,kappa

# gp.x[1,:], gp.y[:,1], gp.LS[1].κ'

# fLS = Figure(size = (1600, 1000))
# ax = Axis(fLS[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
#     xtickalign=0,  ytickalign=0, yticks = tcks)
# arc!(Point2f(0, -1.0 + center(Rf(θe * π / 180, V0), θe * π / 180)), Rf(θe * π / 180, V0), -π, π,
#     linewidth = 3, linestyle = :dash, color = :red, label = "Teo")
# contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0,
#     color = :black, linewidth = 3, label = "Sim");
# fLS[1,2] = Legend(fLS, ax, framevisible = false)
# limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
# colsize!(fLS.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fLS.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fLS)

# fLS0 = Figure(size = (1600, 1000))
# ax = Axis(fLS0[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y",
#     xtickalign=0,  ytickalign=0, yticks = tcks)
# heatmap!(gp.x[1,:], gp.y[:,1], fwd.u[1,1,:,:]')
# contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,1,:,:]', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, num.x[1], num.x[end], num.y[1], num.y[end])
# colsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[1])
# rowsize!(fLS0.layout, 1, widths(ax.scene.viewport[])[2])
# resize_to_layout!(fLS0)
######################################################################################################

