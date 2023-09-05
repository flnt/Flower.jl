using Revise
using Flower
using Peaks

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/stretching/grid3/"

L0 = 10.
n = 128
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

##### Grid 1 #####
# d0 = 0.03
# n = 128
# s1 = stretching(2n÷3, d0, 0.1, 0.55)
# s2 = stretching(n+2n÷3, d0, 0.07, 0.55)

##### Grid 2 #####
# d0 = 0.015
# n = 256
# s1 = stretching(n÷2, d0, 0.075, 0.55)
# s2 = stretching(n+2n÷2, d0, 0.035, 0.55)

##### Grid 3 #####
# d0 = 0.03
# n = 256
# s1 = stretching(2n÷3, d0, 0.1, 0.55)
# s2 = stretching(n+2n÷3, d0, 0.07, 0.55)

##### Grid 4 #####
# d0 = 0.015
# n = 256
# s1 = stretching(6n÷6, d0, 0.075, 0.55)
# s2 = stretching(n+15n÷6, d0, 0.035, 0.55)

##### Grid 5 #####
# d0 = 0.06
# n = 94
# s1 = stretching(6n÷6, d0, 0.2, 0.55)
# s2 = stretching(n+8n÷6, d0, 0.15, 0.55)
# x = vcat(-reverse(s1),s2[2:end])
# y = vcat(-reverse(s1),s1[2:end])

num = Numerical(case = "Cylinder",
    Re = 100.0,
    CFL = 0.5,
    # TEND = 150.0,
    max_iterations = 100,
    x = x,
    y = y,
    R = 0.5,
    u_inf = 1.0,
    save_every = 1,
    shifted = 0.0,
    ϵ = 0.05)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)
fgrid = plot_grid(gp);
lim = 0.75
fgrid2 = plot_grid(gp; limitsx = (-lim, lim), limitsy = (-lim, lim));

# Makie.save(prefix*"grid_close_G1.pdf", fgrid2)

@time MIXED, SOLID, LIQUID = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left = Dirichlet(val = num.u_inf),
    ),
    BC_vL = Boundaries(
        left = Dirichlet(),
    ),
    BC_pL = Boundaries(
        bottom = Dirichlet(),
        right =Dirichlet(),
        top = Dirichlet(),
    ),
    time_scheme = CN,
    navier_stokes = true,
    ns_advection = true,
    ns_solid_phase = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1
)

# make_video(gu, fwd.ux, fwdL.u; title_prefix="u_field",
#         title_suffix="", framerate=240)
# make_video(gv, fwd.uy, fwdL.v; title_prefix="v_field",
#         title_suffix="", framerate=240)
# make_video(gp, fwd.u, fwdL.p; title_prefix="p_field",
#         title_suffix="", framerate=240)

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = DataAspect())
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
resize_to_layout!(fu)

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = DataAspect())
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
resize_to_layout!(fv)

pavg = mean(phL.p[LIQUID].*num.τ)
pstd = std(phL.p[LIQUID].*num.τ)*2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = DataAspect())
heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
resize_to_layout!(fp)

# fCd = Figure(resolution = (1600, 1000))
# colsize!(fCd.layout, 1, Aspect(1, 1.0))
# ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")
# lines!(fwd.Cd)
# limits!(ax, 0, size(fwd.Cd, 1), 0, 2)
# resize_to_layout!(fCd)

# fCl = Figure(resolution = (1600, 1000))
# colsize!(fCl.layout, 1, Aspect(1, 1.0))
# ax = Axis(fCl[1,1], xlabel="it", ylabel="Cl")
# lines!(fwd.Cl)
# limits!(ax, 0, size(fwd.Cl, 1), -1, 1)
# resize_to_layout!(fCl)

# f_grid = plot_grid(gu);

# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/stretching/grid3/"
# suffix = "_Re$(num.Re)_d$(d0)_n$(n).jld"

# save_field(prefix*"data"*suffix, num, phL)

# suffix = "_Re$(num.Re)_d$(d0)_n$(n)"
# make_video(num, fwd, gu, "u"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, minv=0.0, maxv=1.2)

# make_video(num, fwd, gv, "v"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, minv=-0.4, maxv=0.4)

# make_video(num, fwd, gp, "p"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100)

# pks, vals = findmaxima(fwd.Cl)
# f = 1 / ((pks[end-1]-pks[end-2])*num.τ*num.save_every)
# # f = 1 / ((pks[end]-pks[end-1])*num.τ*num.save_every)

# rms_Cl = sqrt(1/(pks[end-1]-pks[end-2]) * sum(fwd.Cl[pks[end-2]:pks[end-1]].^2))
# # rms_Cl = sqrt(1/(pks[end]-pks[end-1]) * sum(fwd.Cl[pks[end-1]:pks[end]].^2))

# pks, vals = findmaxima(fwd.Cd)
# # mean_Cd = mean(fwd.Cd[pks[end-2]:pks[end-1]])
# mean_Cd = mean(fwd.Cd[pks[end-1]:pks[end]])
# # pm_Cd = abs.(fwd.Cd[pks[end-2]] - mean_Cd)
# pm_Cd = abs.(fwd.Cd[pks[end-1]] - mean_Cd)

# # x_cen = getproperty.(gp.geoL.centroid[MIXED], :x)
# # y_cen = getproperty.(gp.geoL.centroid[MIXED], :y)
# # x_mix = gp.x[MIXED] .+ x_cen .* gp.dx[MIXED]
# # y_mix = gp.y[MIXED] .+ y_cen .* gp.dy[MIXED]
# # θ = atan.(y_mix, x_mix) .* 180 ./ pi

# # perm = sortperm(θ)
# # θp = θ[perm]
# # pp = phL.p[MIXED][perm]
# # Vp = gp.geoL.cap[:,:,5][MIXED][perm]

# # θp = θp[pp .!= 0.0]
# # Vp = Vp[pp .!= 0.0]
# # pp = pp[pp .!= 0.0]

# # fpt = Figure(resolution = (1600, 1000))
# # ax = Axis(fpt[1,1], xlabel="θ", ylabel="p", xticks = -180:45:180)
# # lines!(θp, pp, linewidth = 3)
# # lines!(θp, Vp, linewidth = 3)
# # resize_to_layout!(fpt)

# CairoMakie.save(prefix*"grid_cylinder.pdf", fgrid)
# CairoMakie.save(prefix*"Cd.png", fCd)
# CairoMakie.save(prefix*"Cl.png", fCl)
# # CairoMakie.save(prefix*"p_trace.png", fpt)

# open(prefix*"coeffs.txt", "w") do file
#     println(file, "St = "*"$f")
#     println(file, "rms Cl = "*"$rms_Cl")
#     println(file, "mean Cd = "*"$mean_Cd")
#     println(file, "pm Cd = "*"$pm_Cd")
#     close(file)
# end





# # data = load("/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/stretching/grid4/data_Re100.0_d0.015_n256.jld")
# # load_phase!(data, phL)

# # ω = vorticity(gp, phL)

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect = DataAspect())
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# poly!(Circle(Point2f(0, 0), 0.5), color = :grey)
# cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding=0)
# limits!(ax, -2, 10, -6, 6)
# resize_to_layout!(fu)

# fv = Figure(resolution = (1600, 1000))
# colsize!(fv.layout, 1, Aspect(1, 1.0))
# ax = Axis(fv[1,1], aspect = DataAspect())
# hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
# contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# poly!(Circle(Point2f(0, 0), 0.5), color = :grey)
# cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding=0)
# limits!(ax, -2, 10, -6, 6)
# resize_to_layout!(fv)

# # fp = Figure(resolution = (1600, 1000))
# # colsize!(fp.layout, 1, Aspect(1, 1.0))
# # ax = Axis(fp[1,1], aspect = DataAspect())
# # hmap = heatmap!(gp.x[1,:], gp.y[:,1], ω', colorrange=(-2,2))
# # contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# # poly!(Circle(Point2f(0, 0), 0.5), color = :grey)
# # cbar = fp[1,2] = Colorbar(fp, hmap, labelpadding=0)
# # limits!(ax, -2, 10, -6, 6)
# # resize_to_layout!(fp)

# Makie.save(prefix*"u_cylinder.pdf", fu)
# Makie.save(prefix*"v_cylinder.pdf", fv)
# Makie.save(prefix*"vort_cylinder.pdf", fp)


# e = [1.3875575017061683, 1.375320201758595, 1.3713]#, 1.3561427711809102]
# npts = [1/0.06, 1/0.03, 1/0.015]
# r = npts[end]/npts[end-1]
# e .= Richardson_extrapolation(e, r)

# mutable struct conv_regression
#     yreg::Array{Float64, 1}
#     coef0::Float64
#     coef1::Float64
# end
# function regression(x, y, x_reg)
#     coeffs = fit(log.(x), log.(y), 1).coeffs
#     y_reg = exp.(coeffs[2].*log.(x_reg) .+ coeffs[1])
#     result = conv_regression(y_reg, -coeffs[1], -coeffs[2])

#     return result
# end

# conv = regression(npts, e, npts)

# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/cylinder/stretching/"

# fconv = Figure(resolution = (1300, 1000))
# colsize!(fconv.layout, 1, Aspect(1, 1.0))
# ax = Axis(fconv[1,1], aspect = 1, xscale=log10, yscale=log10,
#             xlabel="Points per diameter", ylabel="Error"
# )
# ax.xticks = [16,33,66]
# lines!(npts, conv.yreg, label="Conv. order = $(@sprintf("%.3f", conv.coef1))", color=:black, linewidth=3.0)
# scatter!(npts, e, marker=:rect, color=:black, markersize=20)
# axislegend(position = :lb)
# fconv = current_figure()

# Makie.save(prefix*"conv_cd.pdf", fconv)
