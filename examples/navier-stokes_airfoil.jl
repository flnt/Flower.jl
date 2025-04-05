using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/airfoil/Re500/"

function suction_side(x, t)
    y =  5*t/100*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
end

function pressure_side(x, t)
    y = -5*t/100*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
end

function naca(n, t)
    x = LinRange(0.0, 1.0, n)

    y_suc = suction_side.(x, t)
    y_pres = pressure_side.(reverse(x), t)

    x_airfoil = vcat(x, reverse(x)) .- 0.5
    y_airfoil = vcat(y_suc, y_pres)

    return x_airfoil, y_airfoil
end

fontsize_theme = Theme(fontsize = 60)
set_theme!(fontsize_theme)

# NACA 0010
x_airfoil, y_airfoil = naca(1000, 10)

θ = -30*π/180
rot_mat = [[cos(θ), sin(θ)] [-sin(θ), cos(θ)]]

for (i,(xi, yi)) in enumerate(zip(x_airfoil, y_airfoil))
    x_airfoil[i], y_airfoil[i] = rot_mat * [x_airfoil[i], y_airfoil[i]]
end

# d0 = 0.01
# n = 240
# s1 = stretching(4n÷6, d0, 1.0, 0.55)
# s2 = stretching(n+n, d0, 0.07, 0.55)
# x = vcat(-reverse(s1),s2[2:end])
# y = vcat(-reverse(s1),s1[2:end])

d0 = 0.01
n = 256
s1 = stretching(7n÷6, d0, 0.075, 0.55)
s2 = stretching(n+8n÷4, d0, 0.035, 0.55)
x = vcat(-reverse(s1),s2[2:end])
y = vcat(-reverse(s1),s1[2:end])

num = Numerical(case = "Airfoil",
    Re = 500.0,
    CFL = 0.25,
    TEND = 80.0,
    # TEND = 0.1,
    max_iterations = 0,
    x = x,
    y = y,
    R = 0.5,
    u_inf = 1.0,
    save_every = 30,
    x_airfoil = x_airfoil,
    y_airfoil = y_airfoil,
    ϵ = 0.05
)

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
fgrid = plot_grid(gp; stepx = gp.nx, stepy=gp.ny);
Makie.save(prefix*"grid_airfoil.pdf", fgrid)

# lim = 0.75
# fgrid2 = plot_grid(gp; limitsx = (-lim, lim), limitsy = (-lim, lim), stepx = 2, stepy = 2)
# Makie.save(prefix*"grid_close_airfoil.pdf", fgrid2)

# println("$(min(gp.dx...))$(max(gp.dx...)) $(x[end,end]) $(x[1,1])")

# @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
#     opS, opL, phS, phL, fwd,
#     BC_uL = Boundaries(
#         left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
#     ),
#     BC_vL = Boundaries(
#         left = Boundary(t = dir, f = dirichlet, val = 0.0),
#     ),
#     BC_pL = Boundaries(
#         bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
#         right = Boundary(t = dir, f = dirichlet, val = 0.0),
#         top = Boundary(t = dir, f = dirichlet, val = 0.0),
#     ),
#     stefan = false,
#     advection = false,
#     heat = false,
#     navier_stokes = true,
#     ns_advection = true,
#     ns_solid_phase = false,
#     ns_liquid_phase = true,
#     verbose = true,
#     show_every = 1
# )

# # # tcks = -20.0:2:20.0

# # # fu = Figure(resolution = (1600, 1000))
# # # colsize!(fu.layout, 1, Aspect(1, 1.0))
# # # ax = Axis(fu[1,1], aspect = DataAspect())
# # # hmap = heatmap!(gu.x[1,:], gu.y[:,1], fwd.Uxsave[end-1,:,:]')
# # # contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# # # resize_to_layout!(fu)

# # # fv = Figure(resolution = (1600, 1000))
# # # colsize!(fv.layout, 1, Aspect(1, 1.0))
# # # ax = Axis(fv[1,1], aspect = DataAspect())
# # # hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
# # # contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# # # resize_to_layout!(fv)

# # # pavg = mean(phL.p[LIQUID].*num.τ)
# # # pstd = std(phL.p[LIQUID].*num.τ)*2

# # # fp = Figure(resolution = (1600, 1000))
# # # colsize!(fp.layout, 1, Aspect(1, 1.0))
# # # ax = Axis(fp[1,1], aspect = DataAspect())
# # # heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
# # # contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# # # resize_to_layout!(fp)

# # # fCd = Figure(resolution = (1600, 1000))
# # # colsize!(fCd.layout, 1, Aspect(1, 1.0))
# # # ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")
# # # lines!(fwd.Cd)
# # # limits!(ax, 0, size(fwd.Cd, 1), 0, 2)
# # # resize_to_layout!(fCd)

# # # fCl = Figure(resolution = (1600, 1000))
# # # colsize!(fCl.layout, 1, Aspect(1, 1.0))
# # # ax = Axis(fCl[1,1], xlabel="it", ylabel="Cl")
# # # lines!(fwd.Cl)
# # # limits!(ax, 0, size(fwd.Cl, 1), 0, 2)
# # # resize_to_layout!(fCl)

# # suffix = "_Re$(num.Re).jld"

# # # save_field(prefix*"data"*suffix, num, phL)

# # # suffix = "_Re$(num.Re).jld"
# # # lim = 6.0
# # # make_video(num, fwd, gu, "u"; title_prefix=prefix,
# # #         title_suffix=suffix, framerate=100, minv=0.0, maxv=1.2, limitsx = (-lim, lim), limitsy = (-lim, lim))
# # # make_video(num, fwd, gv, "v"; title_prefix=prefix,
# # #         title_suffix=suffix, framerate=100, minv=-0.4, maxv=0.4, limitsx = (-lim, lim), limitsy = (-lim, lim))
# # # make_video(num, fwd, gp, "p"; title_prefix=prefix,
# # #         title_suffix=suffix, framerate=100, limitsx = (-lim, lim), limitsy = (-lim, lim))

# # # using Peaks
# # # pks, vals = findmaxima(fwd.Cl)
# # # f = 1 / ((pks[end-1]-pks[end-2])*num.τ*num.save_every)
# # # # rms_Cl = sqrt(1/(pks[end-1]-pks[end-2]) * sum(fwd.Cl[pks[end-2]:pks[end-1]].^2))
# # # mean_Cl = mean(fwd.Cl[pks[end-2]:pks[end-1]])
# # # # mean_Cl = mean(fwd.Cl)

# # # pks, vals = findmaxima(fwd.Cd)
# # # # f = 1 / ((pks[end-1]-pks[end-2])*num.τ*num.save_every)
# # # # rms_Cd = sqrt(1/(pks[end-1]-pks[end-2]) * sum(fwd.Cd[pks[end-2]:pks[end-1]].^2))
# # # mean_Cd = mean(fwd.Cd[pks[end-2]:pks[end-1]])
# # # # mean_Cd = mean(fwd.Cd)

# # # CairoMakie.save(prefix*"u.pdf", fu)
# # # CairoMakie.save(prefix*"v.pdf", fv)
# # # CairoMakie.save(prefix*"p.pdf", fp)

# # # CairoMakie.save(prefix*"grid.pdf", fgrid)
# # # CairoMakie.save(prefix*"Cd.pdf", fCd)
# # # CairoMakie.save(prefix*"Cl.pdf", fCl)

# # # open(prefix*"coeffs.txt", "w") do file
# # #     println(file, "St = "*"$f")
# # #     println(file, "mean Cl = "*"$mean_Cl")
# # #     println(file, "mean Cd = "*"$mean_Cd")
# # #     close(file)
# # # end

# data = load("/Users/alex/Documents/PhD/Cutcell/New_ops/navier-stokes/airfoil/Re500/data_Re500.0.jld")
# load_phase!(data, phL)

# minu = min(phL.u...)
# minv = min(phL.v...)

# # phL.u[SOLID_u] .= NaN
# # phL.v[SOLID_v] .= NaN

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect = DataAspect())
# hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u', colorrange=(minu, 1.2))
# contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding=0)
# limits!(ax, -2, 10, -6, 6)
# resize_to_layout!(fu)

# fv = Figure(resolution = (1600, 1000))
# colsize!(fv.layout, 1, Aspect(1, 1.0))
# ax = Axis(fv[1,1], aspect = DataAspect())
# hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v', colorrange=(-1.0, 1.0))
# contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding=0)
# limits!(ax, -2, 10, -6, 6)
# resize_to_layout!(fv)

# Makie.save(prefix*"u_airfoil.pdf", fu)
# Makie.save(prefix*"v_airfoil.pdf", fv)