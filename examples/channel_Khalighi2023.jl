using Revise
using Flower
# using CairoMakie


prefix="/local/home/pr277828/flower/channel3/"

isdir(prefix) || mkdir(prefix)

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

L0 = 1e-4 
L0y=L0
n = 64

x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)

# function f(x,v_inlet,L0)
#     return 8/3*v_inlet*x/L0*(1-x/L0)
# end

function f(x,v_inlet,L0)
    return 6*v_inlet*x/L0*(1-x/L0)
end


v_inlet=6.7e-4

rho=1258.0
mu=6.7e-7
# Re = 1.0
# CFL = 1.0
CFL=0.5
Re=rho*v_inlet*L0/mu

max_iterations=500
max_iterations=1


# v_inlet = 


########################################################################

num = Numerical(
    case = "Nothing",
    x = x,
    y = y,
    Re = Re,
    CFL = CFL,
    max_iterations = max_iterations,
    save_every = 10,
    ϵ = 0.05,
    mu1=mu,
    mu2=mu,
    rho1=rho,
    rho2=rho,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)



gp.LS[1].u .= 1.0

λ = 1e-2
R = L0y / 2.0
########################################
# uPoiseuille = R^2 .- gu.y[:,1].^2

# phL.u .= 1.0
# phL.v .= 0.0
##############################

vPoiseuille = zeros(gv)
@unpack x, nx, ny, ind = gv
vPoiseuille=f.(x,v_inlet,L0)

phL.v .=vPoiseuille
# phL.v .=v_inlet
phL.u .= 0.0
phL.p .= 0.0

vPoiseuilleb=f.(gv.x[1,:],v_inlet,L0)


@unpack τ,CFL,Δ,Re,θd=num
print(@sprintf "dt %.2e CFL %.2e CFL*Δ %.2e CFL*Δ^2*Re %.2e Re %.2e θd %.2e\n" τ CFL CFL*Δ CFL*Δ^2*Re Re θd)
# τ=CFL*Δ/v_inlet
# num.τ=τ
print(@sprintf "dt %.2e \n" τ)

#Re<2000


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
    BC_int = [Wall()],
    time_scheme = CN,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    adapt_timestep_mode = 1,
)

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
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=velscale,
        scalelabel=L"\times 10^{-4}", xticks=xticks, yticks=yticks, scalticks=0:2:10)
make_video_vec(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v",
        title_suffix="", framerate=240, 
        xlabel=xlabel, ylabel=ylabel, xscale=xscale, yscale=yscale, scalscale=velscale, 
        scalelabel=L"\times 10^{-4}", xticks=xticks, yticks=yticks, scalticks=0:2:10)


# function kinetic_energy(fwdL, gp, gu, gv)
#     its = size(fwdL.u)[1]
#     energy = zeros(its)
#     _energy = zeros(gp)

#     for i in 1:its
#         _energy .= (
#             (fwdL.u[i,:,2:end].^2.0 .* gu.geoL.dcap[:,2:end,6] .+ 
#             fwdL.u[i,:,1:end-1].^2.0 .* gu.geoL.dcap[:,1:end-1,6]) ./ 
#             (gu.geoL.dcap[:,1:end-1,6] .+ gu.geoL.dcap[:,2:end,6] .+ 1e-8)
#         )
#         _energy .+= (
#             (fwdL.v[i,2:end,:].^2.0 .* gv.geoL.dcap[2:end,:,7] .+ 
#             fwdL.v[i,1:end-1,:].^2.0 .* gv.geoL.dcap[1:end-1,:,7]) ./
#             (gv.geoL.dcap[1:end-1,:,7] .+ gv.geoL.dcap[2:end,:,7] .+ 1e-8)
#         )
#         energy[i] = sum(_energy .* gp.geoL.dcap[:,:,5])
#     end

#     return energy ./ energy[2]
# end

print(gp.LS[1].geoL.dcap[1,:,5])
print(gp.LS[1].geoL.dcap[:,1,5])

print("Delta",num.Δ, (num.Δ)^2)