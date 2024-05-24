using Revise
using Flower

# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/contact_line/sessile/"
prefix="/local/home/pr277828/flower/sessile2/"

prefix="/local/home/pr277828/flower/"

folder="electrolysis_sessile_1"

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


Rf(θ, V) = sqrt(V / (θ - sin(θ) * cos(θ)))
RR0(θ) = sqrt(π / (2 * (θ - sin(θ) * cos(θ))))
center(r, θ) = r * cos(π - θ)

h0 = 0.5
L0x = 3.0
L0y = 2.0
n = 96

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
y = collect(LinRange(-L0y / 2, 0, n ÷ 3 + 1))


# function run_sessile(θe = 90)

θe= 90

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

max_its=10

# save_every = max_its÷100
save_every = 1
pres0 = 1e5

_θe = acos((0.5 * diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
# _θe = acos((diff(y)[1] + cos(θe * π / 180) * h0) / h0) * 180 / π
println("θe = $(_θe)")
num = Numerical(
    case = "Planar",
    x = x,
    y = y,
    Re = 1.0,
    CFL = CFL,
    max_iterations = max_its,
    u_inf = 0.0,
    v_inf = 0.0,
    save_every = save_every,
    # save_every = 1,
    reinit_every = 10,
    nb_reinit = 2,
    δreinit = 10.0,
    σ = 1.0,
    ϵ = 0.05,
    n_ext_cl = n_ext,
    NB = 24,
    pres0=pres0,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

r = 0.5
gp.LS[1].u .= sqrt.(gp.x.^2 + (gp.y .+ L0y / 2).^2) - r * ones(gp)
# gp.LS[1].u .*= -1.0

phL.u .= 0.0
phL.v .= 0.0

phL.p .= pres0

# plot_python_pdf(phL.p, "p0",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure")


@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        bottom = Navier_cl(λ = 1e-2),
        top = Dirichlet(),
    ),
    BC_vL = Boundaries(
        bottom = Dirichlet(),
        top = Dirichlet(),
    ),
    BC_pL = Boundaries(
        left = Dirichlet(val=pres0),
        right = Dirichlet(val=pres0),
    ),
    BC_u = Boundaries(
        bottom = Neumann_cl(θe = _θe * π / 180),
        top = Neumann_inh(),
        left = Neumann_inh(),
        right = Neumann_inh()
    ),
    BC_int = [FreeSurface()],
    time_scheme = FE,
    auto_reinit = true,
    navier_stokes = true,
    ns_advection = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    save_length = true,
)

V0 = 0.5 * π * 0.5^2
Vf = volume(gp.LS[1].geoL)
Vratio = Vf / V0

mean_rad = 1 / abs(mean(gp.LS[1].κ[gp.LS[1].MIXED[5:end-5]]))
RR0_sim = mean_rad / 0.5
RR0_teo = RR0(θe * π / 180)

println("Vratio = $(Vratio)")
println("mean rad = $(mean_rad)")
println("RR0_sim = $(RR0_sim)")
println("RR0_teo = $(RR0_teo)")

suffix = "$(θe)deg_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
# suffix = "$(θe)deg_$(num.max_iterations)_$(n)_reinit$(num.reinit_every)_nb$(num.nb_reinit)"
file = suffix*".jld2"
# save_field(prefix*file, num, gp, phL, fwdL, fwd)

tcks = -num.L0/2:0.5:num.L0
lim = (num.L0 + num.Δ) / 2
lim = 1.0

x_array=gp.x[1,:]
y_array=gp.y[:,1]

plot_levelset=true
isocontour=false#true
cmap = plt.cm.viridis

# if isnothing(current_i)
#     size_frame=size(fwdL.p,1)
# else
#     size_frame=current_i
# end
size_frame=size(fwdL.p,1)


# plot_python_pdf(phL.p, "p",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure")

plot_python_pdf(phL.p, "p",prefix,plot_levelset,isocontour,0,range(pres0*0.9999,pres0*1.0001,length=10),cmap,x_array,y_array,gp,"pressure")

plot_python_pdf((phL.p.-pres0)./pres0, "pnorm",prefix,plot_levelset,isocontour,0,range(-1e-4,1e-4,length=10),cmap,x_array,y_array,gp,"pressure")


# python_movie_zoom(fwdL.p,"p",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
# size_frame,1,gp.nx,1,gp.ny,fwd)


python_movie_zoom((fwdL.p.-pres0)./pres0,"pnorm",prefix,plot_levelset,isocontour,10,range(0,1400,length=8),cmap,x_array,y_array,gp,"pressure",
size_frame,1,gp.nx,1,gp.ny,fwd)

# plot_python_pdf(max.((phL.trans_scal[:,:,1] .-c0_H2)./c0_H2,0.0), "H2test",prefix,
# plot_levelset,isocontour,0,range(0,1400,length=8),cmap,x_array,y_array,gp)



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


