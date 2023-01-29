using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

L0 = 6
n = 128
d0 = L0 / n

# n_const = ceil(Int, 1.3/d0)
# x_const = d0 * n_const / 2

# s0 = range(-x_const, x_const, n_const)
# s1 = stretching((n-n_const)÷2+1, d0, 3d0, 0.0)

# x = vcat(-reverse(s1[2:end]) .- x_const, collect(s0), x_const .+ s1[2:end])
# # y = vcat(-reverse(s1),s1[2:end])

x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

num = Numerical(case = "Cylinder",
    CFL = 1.0,
    Re = 1.0,
    TEND = 0.5,
    x = x,
    y = y,
    max_iterations = 1,
    R = 0.5,
    u_inf = 1,
    ϵ = 0.02)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
# @profview MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
phS, phL, fwd,
BC_uL = Boundaries(
left = Boundary(t = dir, f = dirichlet, val = num.u_inf),
bottom = Boundary(t = dir, f = dirichlet, val = num.u_inf),
top = Boundary(t = dir, f = dirichlet, val = num.u_inf)
),
BC_vL = Boundaries(
left = Boundary(t = dir, f = dirichlet, val = num.v_inf),
bottom = Boundary(t = dir, f = dirichlet, val = num.v_inf),
top = Boundary(t = dir, f = dirichlet, val = num.v_inf)
),
BC_pL = Boundaries(
        # bottom = Boundary(t = dir, f = dirichlet, val = 0.0),
        right = Boundary(t = dir, f = dirichlet, val = 0.0),
        # top = Boundary(t = dir, f = dirichlet, val = 0.0)
        ),
stefan = false,
advection = false,
heat = false,
navier_stokes = true,
ns_advection = false,
ns_solid_phase = false,
ns_liquid_phase = true,
verbose = true,
show_every = 1
)

# tcks = -num.L0/2:2:num.L0
# lim = num.L0 / 2

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
# contour!(gu.x[1,:], gu.y[:,1], gu.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fu)

# fv = Figure(resolution = (1600, 1000))
# colsize!(fv.layout, 1, Aspect(1, 1.0))
# ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
# contour!(gv.x[1,:], gv.y[:,1], gv.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fv)

# pavg = mean(phL.p[LIQUID].*num.τ)
# pstd = std(phL.p[LIQUID].*num.τ)*2

# fp = Figure(resolution = (1600, 1000))
# colsize!(fp.layout, 1, Aspect(1, 1.0))
# ax = Axis(fp[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fp)

# fphi = Figure(resolution = (1600, 1000))
# colsize!(fphi.layout, 1, Aspect(1, 1.0))
# ax = Axis(fphi[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gp.x[1,:], gp.y[:,1], phL.ϕ')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fphi)

# fCd = Figure(resolution = (1600, 1000))
# colsize!(fCd.layout, 1, Aspect(1, 1.0))
# ax = Axis(fCd[1,1], xlabel="it", ylabel="Cd")  # customized as you see fit
# lines!(fwd.Cd)
# limits!(ax, 0, size(fwd.Cd, 1), 0, 100)
# resize_to_layout!(fCd)

# # f_grid = plot_grid(gu; limitsx=(-1.0, 1.0), limitsy=(-1.0, 1.0));


# x_cen = getproperty.(gp.geoL.centroid[MIXED], :x)
# y_cen = getproperty.(gp.geoL.centroid[MIXED], :y)
# x_mix = gp.x[MIXED] .+ x_cen .* gp.dx[MIXED]
# y_mix = gp.y[MIXED] .+ y_cen .* gp.dy[MIXED]
# θ = atan.(y_mix, x_mix) .* 180 ./ pi

# pD = reshape(phL.pD.data[2], (gp.ny, gp.nx))
# perm = sortperm(θ)
# θp = θ[perm]
# pp = pD[MIXED][perm]
# Vp = gp.geoL.cap[:,:,5][MIXED][perm]

# θp = θp[pp .!= 0.0]
# Vp = Vp[pp .!= 0.0]
# pp = pp[pp .!= 0.0]

# fpt = Figure(resolution = (1600, 1000))
# ax = Axis(fpt[1,1], xlabel="θ", ylabel="p", xticks = -180:45:180)
# lines!(θp, pp, linewidth = 3)
# lines!(θp, Vp, linewidth = 3)
# resize_to_layout!(fpt)


# prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/navier_stokes/cylinder/"
# suffix = ""
# make_video(num, fwd, gu, "u"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100)

# # # make_video(num, fwd, gu, "ucorr"; title_prefix=prefix,
# # #         title_suffix=suffix, framerate=100)

# make_video(num, fwd, gv, "v"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100)

# make_video(num, fwd, gp, "ϕ"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100)

# make_video(num, fwd, gp, "p"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100)
