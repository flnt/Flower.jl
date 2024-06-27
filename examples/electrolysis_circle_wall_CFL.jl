using Revise
using Flower
using PrettyTables
# using Printf
# using Interpolations
# using PyCall
# using PyPlot

prefix="/local/home/pr277828/flower/"

folder="electrolysis_circle_wall_CFL"

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

OkabeIto=["#E69F00", #0 orange clair 230, 159, 0
"#56B4E9", #1 bleu clair 86, 180, 233
"#009E73", #2 vert 0, 158, 115
"#F0E442", #3 jaune 240, 228, 66
"#0072B2", #4 bleu 0, 114, 178
"#D55E00", #5 orange 213, 94, 0
"#CC79A7", #6 rose 204, 121, 171
"#000000"] #7 noir 0 0 0

colors=OkabeIto

colors=["#000000" for color in OkabeIto]

colors[1]="#000000"
colors[2]=OkabeIto[5] #bleu
colors[3]=OkabeIto[6] #orange

####################################################################################################
# Sessile from sessile.jl
####################################################################################################
Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)
####################################################################################################

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

# max_its=10
# save_every = max_its÷100
save_every = 1
# pres0 = 1e5
pres0 = 0.0

L0 = 1e-4 
n = 64 
n = 128
# n = 256
# n = 512
# n=1024

max_iter=100
max_iter=10
max_iter = 2

max_iter = 1

# max_iter=5

# max_iter = 1

xlabel = L"x \left(\mu m \right)"
ylabel = L"y \left(\mu m \right)"

xscale = 1e-6
yscale = xscale

xticks = 0:20:100
yticks = 0:20:100

velscale = 1e-4 

concentrationcontour=true
concentrationcontour=false

plot_grid=true
plot_grid=false

plot_levelset=true

cmap = plt.cm.viridis

adapt_timestep_mode=2
adapt_timestep_mode=3

dt0 = 5e-5
dt0 = 2.5e-5 #CFL >1 for n = 128
dt0 = 1.25e-5 #CFL >0.5 for n = 128
dt0 = 6.25e-6 #CFL >0.5 for n = 128
dt0 = 3.125e-6 #radius CFL: 3.54e-01 -7.17e+00% 
# dt0 = 1e-6 #error dnH2 -2.40e-15



epsilon = 0.05
# epsilon = 0.01


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
# radius = 6.0e-6 

h0 = radius


ref_thickness_2d = 4.0 / 3.0 *radius 

# mode_2d = 1 #use equivalent cylinder
# mode_2d = 2 #mol/meter
mode_2d = 3

xcoord = 0.0
ycoord = L0/2.0
# xcoord = -xcoord
# ycoord = -ycoord

mu=6.7e-7
mu1=mu
mu2=mu #TODO
i0=1.0
temperature0=353.0
pres0= 0.0 #1e5
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
plot_xscale = 1e-6
concentration_check_factor = 1e-2 # 1e-1 #1e-2

# elec_cond=1 #TODO
elec_cond=2*Faraday^2*c0_KOH*DKOH/(Ru*temperature0)


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

# nb_saved_scalars=4 
# nb_saved_scalars=5
# nb_saved_scalars=6

nb_saved_scalars=1


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

hl = Highlighter((d,i,j)->d[i,j] isa String, crayon"bold cyan")

diffusion_t = (radius^2)./diffusion_coeff

pretty_table(vcat(
    hcat("Diffusion time",diffusion_t'),
    hcat("Diffusion coef",diffusion_coeff'),
    hcat("Concentration",concentration0')); 
formatters    = ft_printf("%0.2e", 2:4), #not ecessary , 2:4
header = ["","H2", "KOH", "H2O"], 
highlighters=hl)



# printstyled(color=:green, @sprintf "\n Species diffusion timescales: %.2e %.2e %.2e \n" (radius^2)/DH2 (radius^2)/DKOH (radius^2)/DH2O )


# printstyled(color=:green, @sprintf "\n nmol : \n" concentration0[1]*4.0/3.0*pi*radius^3  )

current_radius = radius

p_liq= pres0 #+ mean(veci(phL.pD,grid,2)) #TODO here one bubble
# p_g=p_liq + 2 * sigma / current_radius
p_g=p_liq + sigma / current_radius


c0test = p_g / (temperature0 * Ru) 


printstyled(color=:green, @sprintf "\n c0test: %.2e \n" c0test)

# printstyled(color=:green, @sprintf "\n Mole test: %.2e %.2e\n" concentration0[1]*4.0/3.0*pi*current_radius^3 p_g*4.0/3.0*pi*current_radius^3/(temperature0*num.Ru))




####################################################################################################


# h0 = 0.5
# L0x = 3.0
# L0y = 2.0
# n = 96

# x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
# y = collect(LinRange(-L0y / 2, 0, n ÷ 3 + 1))





####################################################################################################
# not imposing the angle exactly at the boundary but displaced a cell because ghost cells are not used. 
# So the exact contact angle cannot be imposed
# TODO should at some point modify the levelset so it works as the other fields that we have in the code, 
# with the boundary values in addition to the bulk field. That way we could impose it exactly at the boundary
# needs a lot of work, not priority
_θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
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
    case = "Cylinder",#"Planar",
    R = radius,
    max_iterations = max_iter,
    save_every = save_every,#10,#save_every,#10,
    ϵ = epsilon, 
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
    # σ=sigma,   
    # reinit_every = 10,
    # nb_reinit = 2,
    # δreinit = 10.0,
    # n_ext_cl = n_ext,
    # NB = 24,
    plot_xscale = plot_xscale,
    plot_prefix = prefix,
    dt0 = dt0,
    concentration_check_factor = concentration_check_factor,
    )
    # ref_thickness_2d = ref_thickness_2d,

#Initialization
gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)





@unpack x, nx, ny, ind = gv

velocity = 1.0
velocity = 0.0

vPoiseuille = zeros(gv)
vPoiseuille = Poiseuille_favg.(x,v_inlet,L0) .* velocity
vPoiseuilleb = Poiseuille_favg.(gv.x[1,:],v_inlet,L0) .* velocity



phL.v .=vPoiseuille .* velocity
phL.u .= 0.0
phL.p .= pres0 #0.0
phL.T .= temperature0
phS.T .= temperature0

#Initialize for CRASH detection (cf isnan...)
phS.v .= 0.0
phS.u .= 0.0
phS.vD .= 0.0
phS.uD .= 0.0


####################################################################################################
# r = 0.5
# gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
# gp.LS[1].u .*= -1.0
####################################################################################################

gp.LS[1].u .= sqrt.((gp.x .- xcoord).^2 + (gp.y .- ycoord).^2) - radius * ones(gp)
# gp.LS[1].u .*= -1.0


# gp.LS[1].u .= 1.0

# gp.LS[1].u .*= -1.0


su = sqrt.((gv.x .- xcoord).^2 .+ (gv.y .- ycoord).^2)
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

# vecb_L(phL.uD, gu) .= 0.0 #TODO
# vecb_B(phL.uD, gu) .= 0.0
# vecb_R(phL.uD, gu) .= 0.0
# vecb_T(phL.uD, gu) .= 0.0


printstyled(color=:green, @sprintf "\n CFL : %.2e dt : %.2e\n" CFL CFL*L0/n/v_inlet)


xscale = plot_xscale
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

print(@sprintf "Butler-Volmer %.2e %.2e %.2e %.2e\n" i_current[1] -i_current[1]/(2*Faraday*DH2) c0_H2-i_current[1]/(2*Faraday*DH2)*gp.dx[1,1] c0_H2+i_current[1]/(2*Faraday*DH2)*gp.dx[1,1])



######################################################################################################
BC_u = Boundaries(
    bottom = Neumann_inh(),
    top = Neumann_inh(),
    left = Neumann_cl(θe = _θe * π / 180),
    right = Neumann_inh()
)

BC_uS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

BC_vS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

BC_pS = Boundaries(
    left   = Dirichlet(),
    right  = Dirichlet(),
    bottom = Dirichlet(),
    top    = Dirichlet(),
)

@time current_i=run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left   = Dirichlet(),#Navier_cl(λ = 1e-2), #Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(),
        top    = Neumann(val=0.0),
    ),
    BC_uS=BC_uS,
    BC_vL = Boundaries(
        left   = Dirichlet(),
        right  = Dirichlet(),
        bottom = Dirichlet(val = copy(vPoiseuilleb)),
        top    = Neumann(val=0.0),
    ),
    BC_vS=BC_vS,
    BC_pL = Boundaries(
        left   = Neumann(val=0.0),
        right  = Neumann(val=0.0),
        bottom = Neumann(val=0.0),
        top    = Dirichlet(),
    ),
    BC_pS=BC_pS,
    # left = Dirichlet(val=pres0),
    # right = Dirichlet(val=pres0),

    # BC_u = BC_u,
    # BC_int = [FreeSurface()],
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
        left   = Neumann(val=-i_current/(2*Faraday*DH2)), #Dirichlet(val = concentration0[1]), #
        right  = Dirichlet(val = concentration0[1]),
        int    = Dirichlet(val = concentration0[1])), #H2
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[2]),
        top    = Neumann(),
        left   = Neumann(val=-i_current/(2*Faraday*DKOH)),
        right  = Dirichlet(val = concentration0[2]),
        int    = Neumann(val=0.0)), #KOH
         
        BoundariesInt(
        bottom = Dirichlet(val = concentration0[3]),
        top    = Neumann(),
        left   = Neumann(val=i_current/(Faraday*DH2O)),
        right  = Dirichlet(val = concentration0[3]),
        int    = Neumann(val=0.0)), #Dirichlet(val = concentration0[3])),#Neumann(val=0.0)) 
        #H2O #Dirichlet(val = concentration0[3])
    ),

    BC_phi_ele = BoundariesInt(
        left   = Neumann(val=-i_current/elec_cond),
        right  = Dirichlet(),
        bottom = Neumann(val=0.0),
        top    = Neumann(val=0.0),
        int    = Neumann(val=0.0),
    ),

    # auto_reinit = true,
    # save_length = true,
    time_scheme = FE,#CN, #or FE?
    electrolysis = true,
    navier_stokes = true,
    ns_advection=true,#false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    electrolysis_convection = true,  
    electrolysis_liquid_phase = true,
    electrolysis_phase_change = false, #true,
    # electrolysis_phase_change_case = "levelset",
    electrolysis_phase_change_case = "Khalighi",
    electrolysis_reaction = "Butler_no_concentration", #"nothing", #
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

# tcks = -num.L0/2:0.5:num.L0
# lim = (num.L0 + num.Δ) / 2
# lim = 1.0
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



# iLS=1
# rhs_LS = fzeros(gp)
# #test contact angle
# xp,yp,xp0,yp0=BC_LS_test!(gp, gp.LS[iLS].u, gp.LS[iLS].A, gp.LS[iLS].B, rhs_LS, BC_u)

# fig1, ax2 = plt.subplots(layout="constrained")
# CS = ax2.contourf(x_array,y_array,gp.LS[1].u ./xscale, 10, cmap=cmap)

# # print("\n xp",x_array)
# # print("\n yp",y_array)

# # print("\n xp",xp)
# # print("\n yp",yp)
# # print("\n xp0",xp0)
# # print("\n yp0",yp0)

# xp  = xp/xscale
# yp  = yp/xscale
# xp0 = xp0/xscale
# yp0 = yp0/xscale

# ms=2
# mcolor="w"

# plt.scatter(xp,yp,s=ms,c=mcolor)

# mcolor="black"

# plt.scatter(xp0,yp0,s=ms,c=mcolor)



# ax2.set_xlabel(L"$x (\mu m)$")
# ax2.set_ylabel(L"$y (\mu m)$")

# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = fig1.colorbar(CS)
# # cbar.ax.set_ylabel("v")

# cbar.ax.set_title(raw"$( \unit{\um})$")

# # cbar.ax.set_title(L"$10^{-6}$")

# plt.axis("equal")

# plt.savefig(prefix*"contact.pdf")

# plt.close(fig1)
######################################################################################################

#TODO bug when putting plotting instructions in function, range(start,end,step) no longer working as argument so give start,end,length instead

# gen_name=nx
plot_grid = false

# plt_it = [2]

# plt_list =  range(2,current_i,1)

# plt_list = 2:current_i:1
# print(plt_list)

i0 = 1
i1 = 10
# j0 = gp.ny ÷ 2 +1
# j1 = gp.ny ÷ 2 + 10

# j0=67
# j1=76

j0 = gp.ny ÷ 2 +1
j1 = gp.ny ÷ 2 + 10

fontsize = 2 #2
printmode = "val"


# for plt_it = 2:size_frame+1

# end

plt_it = 1
plot_python_pdf(plt_it,fwd.u[1,:,:,:]./xscale, "LS",prefix,
plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"LS",1,gp.nx,1,gp.ny,fwd)

# plot_python_pdf(plt_it,fwdL.p, "p",prefix,
# true,isocontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",1,gp.nx,1,gp.ny,fwd)


plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1],fwdL.trans_scalD[:,:,1], "H2liqlvlzoomfullfield",prefix,
plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)

plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,2],fwdL.trans_scalD[:,:,2], "KOHliqlvlzoomfullfield",prefix,
plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)

plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,3],fwdL.trans_scalD[:,:,3], "H2Oliqlvlzoomfullfield",prefix,
plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)

fig1, ax2 = plt.subplots(layout="constrained")
plt.plot(vecb_L(phL.trans_scalD[:,3], gp),gp.y/xscale)
ax2.set_xlabel(L"$H2O$")
ax2.set_ylabel(L"$y$")

plt.savefig(prefix*"H2O_electrode.pdf")
plt.close(fig1)


fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

# varx = [0, 1, 2]
varx = gp.y/xscale
# vecb_L(phL.trans_scalD[:,3], gp)
label1 = "c(H2O)"
label2 = "Overpotential: -eta"
label3 = "Current"
alpha = 0.5
ls="-" #"--"

# print("current", phL.i_current_mag[:,1])

# print(colors)

twin1 = ax.twinx()
twin2 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin2.spines.right.set_position(("axes", 1.2))
#colors "C0", "C1", "C2"
p1, = ax.plot(varx, phL.trans_scal[:,1,3],colors[1], label=label1,ls=ls,alpha=alpha)
p2, = twin1.plot(varx, phi_ele0 .- phL.phi_ele[:,1], colors[2], label=label2,ls=ls,alpha=alpha)
p3, = twin2.plot(varx, phL.i_current_mag[:,1], colors[3], label=label3,ls=ls,alpha=alpha)


# for l in fig.gca().lines
#     l.set_alpha(alpha)
# end

# p1.set_alpha(alpha)
# p2.set_alpha(alpha)
# p3.set_alpha(alpha)




ax.set(
    # xlim=(0, 2),
    # ylim=(0, 2),
    xlabel="y", ylabel=label1)
twin1.set(
    # ylim=(0, 4), 
ylabel=label2)
twin2.set(
    # ylim=(1, 65), 
ylabel=label3)

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

ax.tick_params(axis="y", colors=p1.get_color())
twin1.tick_params(axis="y", colors=p2.get_color())
twin2.tick_params(axis="y", colors=p3.get_color())

ax.legend(handles=[p1, p2, p3])


plt.savefig(prefix*"test.pdf")
plt.close(fig1)


#############################################################################################################################################"

# for plt_it = 2:current_i+1
for plt_it = 2:size_frame+1


    # printstyled(color=:green, @sprintf "\n plt_it %.5i" plt_it)


    plot_python_pdf(plt_it,fwdL.p, "p",prefix,plot_levelset,isocontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",1,gp.nx,1,gp.ny,fwd)

  
    # plot_python_pdf(phL.p, "p",prefix,plot_levelset,isocontour,0,range(pres0*0.9999,pres0*1.0001,length=10),cmap,x_array,y_array,gp,"pressure")

    # plot_python_pdf((phL.p.-pres0)./pres0, "pnorm",prefix,plot_levelset,isocontour,0,range(-1e-4,1e-4,length=10),cmap,x_array,y_array,gp,"pressure")


    plot_python_pdf(plt_it,max.((fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2_norm",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",0,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)

    plot_python_pdf(plt_it,max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2lvl",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)

    # python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2_norm",prefix,
    # plot_levelset,isocontour,0,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",size_frame,1,gp.nx,1,gp.ny,fwd)

    plot_python_pdf(plt_it,fwdL.trans_scal[:,:,:,1], "H2",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)


    ######################################################################################################



    plot_python_pdf(plt_it,max.((fwdL.trans_scal[:,:,:,2] .-c0_KOH)./c0_KOH,0.0), "KOH_norm",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)

    plot_python_pdf(plt_it,max.((fwdL.trans_scal[:,:,:,3] .-c0_H2O)./c0_H2O,0.0), "H2O_norm",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)


    plot_python_pdf(plt_it,fwdL.trans_scal[:,:,:,2] , "KOH",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)

    plot_python_pdf(plt_it,fwdL.trans_scal[:,:,:,3], "H2O",prefix,
    plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",1,gp.nx,1,gp.ny,fwd)


    
    # plot_python_several_pdf(fwd.trans_scal[:,:,:,1],"H2",true,size_frame)
    # plot_python_several_pdf(fwd.saved_scal[:,:,:,1],"H2massflux",true,size_frame)


    plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,1], "flux1zoom",prefix,
    true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,2], "flux2zoom",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,3], "flux3zoom",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,4], "flux4zoom",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,5], "flux5zoom",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,6], "flux6zoom",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)


    # plot_python_pdf_full2(plt_it,max.((fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2normliqlvlzoomfullfield",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)


    plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1],fwdL.trans_scalD[:,:,1], "H2liqlvlzoomfullfield",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,3],fwdL.trans_scalD[:,:,3], "H2Ofull",prefix,
    # false,concentrationcontour,true,"pcolormesh",0,range(48950,49050,length=10),cmap,x_array,y_array,gp,"concentration",
    # 1,gp.nx,1,3,fwd,fwdL,xscale,fontsize,printmode)

    plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,2],fwdL.trans_scalD[:,:,2], "KOHliqlvlzoomfullfield",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)


    plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,3],fwdL.trans_scalD[:,:,3], "H2Oliqlvlzoomfullfield",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)


    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield2",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,1,10,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield3",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,10,20,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield4",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,20,30,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield5",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,30,40,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield6",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,40,50,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield7",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,50,60,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf_full2(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoomfullfield8",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,60,70,fwd,fwdL,xscale,fontsize,printmode)

    # plot_python_pdf(1,fwd.saved_scal[:,:,:,1], "flux1zoom1",prefix,
    # true,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",i0,i1,j0,j1,fwd)


    # plot_python_pdf(plt_it,max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2normlvlzoom",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)

    plot_python_pdf(plt_it,max.((fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2normliqlvlzoom",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)

    # plot_python_pdf(1,max.((fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2normliqlvlzoom1",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)


    # plot_python_pdf_full2(plt_it,max.((fwdL.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0), "H2normliqlvlzoomfullfield",prefix,
    # plot_levelset,concentrationcontour,true,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)


    # plot_python_pdf(plt_it,fwd.trans_scal[:,:,:,1], "H2lvlzoom",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)

    plot_python_pdf(plt_it,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoom",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)


    # @views vec = reshape(veci(fwdL.trans_scalD[plt_it,:,1],gp,2),gp)

    vec = zeros(gp)
    # vec = reshape(veci(phL.trans_scalD[:,1],gp,2),gp)   
    vec = reshape(veci(fwdL.trans_scalD[plt_it,:,1],gp,2),gp)


    # print("\n testvec",maximum(vec),size(vec))

    # vec = phL.trans_scal[:,:,1]
    # veci(phL.trans_scalD[:,1],grid,2)

    plot_python_pdf_full2(-plt_it,vec,fwdL.trans_scalD[:,:,1], "H2_int_liqlvlzoomfullfield",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)

    plot_python_pdf(-plt_it,vec, "H2_int_liqlvlzoom",prefix,
    plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)


    # plot_python_pdf_full2(-plt_it,vec, "H2_int_liqlvlzoomfullfield",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
    # i0,i1,j0,j1,fwd,fwdL,xscale,fontsize,printmode)


    # plot_python_pdf(1,fwdL.trans_scal[:,:,:,1], "H2liqlvlzoom1",prefix,
    # plot_levelset,concentrationcontour,true,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",i0,i1,j0,j1,fwd)

    # plot_python_pdf(plt_it,fwd.u[1,:,:,:]./xscale, "LSzoom",prefix,
    # true,concentrationcontour,true,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"LS",i0,i1,j0,j1,fwd)

    # plot_python_pdf(1,fwd.u[1,:,:,:]./xscale, "LSzoom1",prefix,
    # true,concentrationcontour,true,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"LS",i0,i1,j0,j1,fwd)



    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,1], "flux1phL",prefix,
    # plot_levelset,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)


    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,1], "flux1",prefix,
    # plot_levelset,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,2], "flux2",prefix,
    # plot_levelset,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,3], "flux3",prefix,
    # plot_levelset,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,4], "flux4",prefix,
    # plot_levelset,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)


    plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,1], "flux1noLS",prefix,
    false,concentrationcontour,plot_grid,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,2], "flux2noLS",prefix,
    # false,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,3], "flux3noLS",prefix,
    # false,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    # plot_python_pdf(plt_it,fwd.saved_scal[:,:,:,4], "flux4noLS",prefix,
    # false,concentrationcontour,plot_grid,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",1,gp.nx,1,gp.ny,fwd)

    
end

#############################################################################################################################################"


python_movie_zoom(fwdL.v,"v",prefix,false,isocontour,"pcolormesh",10,range(0,1400,length=8),cmap,xv,yv,gv,"v",
size_frame,1,gv.nx,1,gv.ny,fwd)

python_movie_zoom(fwdL.p,"p",prefix,plot_levelset,isocontour,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)
# python_movie_zoom((fwdL.p.-pres0)./pres0,"pnorm",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
# size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2lvl",prefix,
plot_levelset,isocontour,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2norm",prefix,
plot_levelset,isocontour,"pcolormesh",0,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",size_frame,1,gp.nx,1,gp.ny,fwd)



python_movie_zoom(max.((fwd.trans_scal[:,:,:,2] .-c0_KOH)./c0_KOH,0.0),"KOH_norm",prefix,
plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
size_frame,1,gp.nx,1,gp.ny,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,3] .-c0_H2O)./c0_H2O,0.0),"H2O_norm",prefix,
plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
size_frame,1,gp.nx,1,gp.ny,fwd)

# python_movie_zoom(fwd.saved_scal[:,:,:,1],"flux1",prefix,
# plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
# size_frame,1,gp.nx,1,gp.ny,fwd)

# python_movie_zoom(fwd.saved_scal[:,:,:,2],"flux2",prefix,
# plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
# size_frame,1,gp.nx,1,gp.ny,fwd)

# python_movie_zoom(fwd.saved_scal[:,:,:,1],"flux1noLS",prefix,
# false,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
# size_frame,1,gp.nx,1,gp.ny,fwd)



python_movie_zoom(fwd.saved_scal[:,:,:,1],"flux1zoom",prefix,
plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"flux",
size_frame,i0,i1,j0,j1,fwd)

python_movie_zoom(max.((fwd.trans_scal[:,:,:,1] .-c0_H2)./c0_H2,0.0),"H2lvlzoom",prefix,
plot_levelset,false,"pcolormesh",10,range(0,1400,length=8),cmap,x_array,y_array,gp,"concentration",
size_frame,i0,i1,j0,j1,fwd)





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


fig1, ax2 = plt.subplots(layout="constrained")

plt.plot(gp.y,vecb_L(phL.trans_scalD[:,3], gp))
ax2.set_xlabel(L"$y$")
ax2.set_ylabel(L"$H2O$")
plt.axis("equal")

plt.savefig(prefix*"H2O_electrode.pdf")
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

