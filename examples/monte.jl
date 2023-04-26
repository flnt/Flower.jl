using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/conferences/2023/Montestigliano_2023/project/"

fontsize_theme = Theme(fontsize = 70)
set_theme!(fontsize_theme)

Lx = 10
Ly = 1
n = 64
nx = Lx * n
# nx = 7 * n
ny = n

x = LinRange(-Lx/2, Lx/2, nx+1)
y = LinRange(-Ly/2, Ly/2, ny+1)

# d0 = 0.009
# s1 = stretching(nx÷2, d0, 0.075, 0.75)
# x = vcat(-reverse(s1),s1[2:end])

nu = 0.0044

Re = 100.0
refU = 0.0006 * Re
refT = 0.0075 / refU

T = 1 / refT

TEND = 3 * T

num = Numerical(
    case = "Monte",
    x = x,
    y = y,
    Re = Re,
    CFL = 0.5,
    TEND = TEND,
    # max_iterations = 6000,
    u_inf = 1.0,
    v_inf = 0.0,
    save_every = 100,
    ϵ = 0.05,
    NB = 12,
    subdomains = 2,
    overlaps = 4
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

uy = zeros(ny)
uy[ny÷4+2:end] .= 1.0/(4nu) * ((3*Ly/8) ^ 2 .- (gp.y[ny÷4+2:end,1] .- Ly./8.) .^ 2)
uy ./= findmax(uy)[1]

data = load(prefix*"data_100.0_3T.jld")
phL = data["phL"]
fwdL = data["fwdL"]

# phL.u .= num.u_inf
# phL.v .= num.v_inf
# phL.T .= 0.

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd, fwdS, fwdL,
    BC_uL = Boundaries(
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet),
        left = Boundary(t = pul, f = dirichlet, val=uy),
        right = Boundary(t = neu, f = neumann)
    ),
    BC_vL = Boundaries(
        bottom = Boundary(t = dir, f = dirichlet),
        top = Boundary(t = dir, f = dirichlet),
        left = Boundary(t = dir, f = dirichlet),
        right = Boundary(t = neu, f = neumann)
    ),
    BC_pL = Boundaries(
        left = Boundary(t = neu, f = neumann),
        right = Boundary(t = dir, f = dirichlet),
        bottom = Boundary(t = neu, f = neumann),
        top = Boundary(t = neu, f = neumann)
    ),
    BC_TL = Boundaries(
        left = Boundary(t = dir, f = dirichlet, val=1.0),
    ),
    BC_u = Boundaries(
    ),
    stefan = false,
    advection = false,
    heat = true,
    heat_liquid_phase = true,
    heat_convection = true,
    navier_stokes = true,
    ns_advection = true,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    T = T,
    nT0 = 0,
)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xlabel="x", ylabel="y")  # customized as you see fit
heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1)  # customized as you see fit
heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.p')
contour!(gp.x[1,:], gp.y[:,1], fwd.u[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fT = Figure(resolution = (1600, 1000))
colsize!(fT.layout, 1, Aspect(1, 1.0))
ax = Axis(fT[1,1], aspect = 1)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.T')
contour!(gp.x[1,:], gp.y[:,1], fwd.u[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fT)

limx = num.x[end]
limy = num.y[end]

suffix = "_Re$(num.Re)_3T"

make_video(gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field",
        title_suffix=suffix, framerate=24, limitsx=(-limx-num.Δ/2,limx+num.Δ/2), limitsy=(-limy,limy))#, minv=0.0, maxv=1.0)
make_video(gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx,limx), limitsy=(-limy-num.Δ/2,limy+num.Δ/2))
make_video(gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx,limx), limitsy=(-limy,limy))
make_video(gp, fwd.u, fwdL.T; title_prefix=prefix*"T_field",
        title_suffix=suffix, framerate=240, limitsx=(-limx,limx), limitsy=(-limy,limy), minv=0.3, maxv=1.0)

JLD.save(prefix*"data_$(num.Re)_3T.jld", "phL", phL, "fwdL", fwdL)


# RpointRe100 = CartesianIndex(1,353)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], fwdL.u[225,:,:]')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# # scatter!([gp.x[RpointRe100]], [gp.y[RpointRe100]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding=0, height = Relative(2/4))
# resize_to_layout!(fu)
# Makie.save(prefix*"u_rec_Re1000_syst.png", fu)

# RpointRe500 = CartesianIndex(1,418)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# scatter!([gp.x[RpointRe500]], [gp.y[RpointRe500]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# resize_to_layout!(fu)
# Makie.save(prefix*"u_rec_Re500.png", fu)

# RpointRe1000 = CartesianIndex(1,473)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# scatter!([gp.x[RpointRe1000]], [gp.y[RpointRe1000]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# resize_to_layout!(fu)
# Makie.save(prefix*"u_rec_Re1000.png", fu)

# RpointRe1000 = CartesianIndex(1,235)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# scatter!([gp.x[RpointRe1000]], [gp.y[RpointRe1000]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# resize_to_layout!(fu)
# Makie.save(prefix*"u_rec_Re1000_mesh2.png", fu)

# RpointRe3500 = CartesianIndex(1,561)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# scatter!([gp.x[RpointRe3500]], [gp.y[RpointRe3500]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# resize_to_layout!(fu)
# Makie.save(prefix*"c_Re3500.png", fu)

# RpointRe5000 = CartesianIndex(1,598)
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=DataAspect(), xlabel="x", ylabel="y")  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# scatter!([gp.x[RpointRe5000]], [gp.y[RpointRe5000]], color=:red, markersize=20)
# limits!(ax, -1.0, 5.0, -0.5, 0.5)
# resize_to_layout!(fu)
# Makie.save(prefix*"u_rec_Re5000.png", fu)

# RpointRe5000 = CartesianIndex(1,1)
# for t = 1:num.max_iterations÷num.save_every
#     for i = 322:640
#         if fwdL.u[t,1,i]<0 && fwdL.u[t,1,i+1]>0
#             RpointRe5000 = CartesianIndex(1,i)
#             println("$(t) $(RpointRe5000)")
#         end
#     end
# end

# RpointRe5000 = CartesianIndex(1,1)
# for i = 322:640
#     if fwdL.u[225,1,i]<0 && fwdL.u[225,1,i+1]>0
#         RpointRe5000 = CartesianIndex(1,i)
#         println("$(RpointRe5000)")
#     end
# end

ls2 = [12, 25, 21, 19, 17] ./ 10
fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect=1, xlabel="Re", ylabel="Length")
lines!(Res, ls, linewidth=3, label="DNS")
lines!(Res, ls2, linewidth=3, label="RANS")
scatter!([940], [1.27], markersize=30, color=:red, label="Ref1 thrombus length")
scatter!([666.6], [2.0], markersize=30, color=:black, label="Ref2 Recirc. length")
axislegend(position=:lt)
resize_to_layout!(fu)
Makie.save(prefix*"le_vs_re.png", fu)

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(phL.u[:,300], gu.y[:,300])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_before_Re1000.png", fu)

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(phL.u[:,340], gu.y[:,340])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_after_Re1000.png", fu)

# it = 55
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(fwdL.u[it,:,300], gu.y[:,300])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_before_Re100_3T_dias.png", fu)

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(fwdL.u[it,:,340], gu.y[:,340])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_after_Re100_3T_dias.png", fu)

# it = 290
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(fwdL.u[it,:,300], gu.y[:,300])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_before_Re1000_3T_dias.png", fu)

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(fwdL.u[it,:,340], gu.y[:,340])
# resize_to_layout!(fu)
# Makie.save(prefix*"u_profile_after_Re1000_3T_dias.png", fu)


# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="t", ylabel="Factor")
# lines!([0.0, 0.2], [1.0, 1.0], color=:black)
# lines!([0.2, 0.25], [1.0, 0.0], color=:black)
# lines!([0.25, 0.95], [0.0, 0.0], color=:black)
# lines!([0.95, 1.0], [0.0, 1.0], color=:black)
# resize_to_layout!(fu)
# Makie.save(prefix*"pulsating.png", fu)


# uy = zeros(ny)
# uy[ny÷4+2:end] .= 1.0/(4nu) * ((3*Ly/8) ^ 2 .- (gp.y[ny÷4+2:end,1] .- Ly./8.) .^ 2)
# uy ./= findmax(uy)[1]
# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect=1, xlabel="Vel", ylabel="y")
# lines!(uy, gu.y[:,1])
# resize_to_layout!(fu)
# Makie.save(prefix*"inlet_profile.png", fu)