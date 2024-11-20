using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ratio = 8
L0 = 1.
n = 64
nx = ratio * n
ny = n

x = collect(LinRange(-ratio*L0 / 2., ratio*L0 / 2., nx + 1))
y = collect(LinRange(-L0 / 2., L0 / 2., ny + 1))

num = Numerical(
    x = x,
    y = y,
    Re = 1.0,
    CFL = 0.5,
    # TEND = 1.0,
    max_iterations = 1000,
    u_inf = 0.0,
    v_inf = 0.0,
    save_every = 10,
    reinit_every = 1,
    nb_reinit = 10,
    δreinit = 0.65,
    g = 0.0,
    β = 0.0,
    σ = 0.0,
    ϵ = 0.01,
    NB = 24,
    nLS = 1,
)

Ra = 1e6
λ = 1.0
H0 = 0.05
T1 = 1.0
T2 = 0.0

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

@. gp.LS[1].u = -gp.y - L0/2 + H0 + 0.0001


# for II in CartesianIndices(gp.x)
#     ystar = gp.y[II] + L0/2
#     if (ystar <= H0)
#         phL.T[II] = 1 +(num.θd-1)*ystar/H0;
#     elseif (ystar > H0)
#         phS.T[II] = num.θd*(ystar-1)/(H0-1)
#     end
# end

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = true,
    BC_TS = Boundaries(left = Periodic(),
                        right = Periodic(),
                        top = Dirichlet(val = T2)),
    BC_TL = Boundaries(left = Periodic(),
                        right = Periodic(),
                        bottom = Dirichlet(val = T1)),
    BC_uL = Boundaries(left = Periodic(),
                        right = Periodic()),
    BC_vL = Boundaries(left = Periodic(),
                        right = Periodic()),
    BC_pL = Boundaries(left = Periodic(),
                        right = Periodic()),
    BC_u = Boundaries(left = Periodic(),
                        right = Periodic()),
    BC_int = [Stefan()],
    time_scheme = FE,
    auto_reinit = true,
    heat = true,
    heat_convection = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    Ra = Ra,
    λ = λ
)

xtcks = -ratio*L0/2:0.5:ratio*L0/2
ytcks = -L0/2:0.5:L0/2
lim = L0 / 2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = ratio, xticks = xtcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.T', colorrange=(num.θd, T1), colormap=Reverse(:ice))
contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, x[1], x[end], y[1], y[end])
resize_to_layout!(fp)

fp = current_figure()

# fu = Figure(resolution = (1600, 1000))
# colsize!(fu.layout, 1, Aspect(1, 1.0))
# ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gu.x[1,:], gu.y[:,1], phL.u', colorrange=(minimum(phL.u), maximum(phL.u)))
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fu)

# fu = current_figure()

# fv = Figure(resolution = (1600, 1000))
# colsize!(fv.layout, 1, Aspect(1, 1.0))
# ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gv.x[1,:], gv.y[:,1], phL.v', colorrange=(minimum(phL.v), maximum(phL.v)))
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fv)

# fv = current_figure()


# # prefix = "/Users/alex/codes/Flower.jl/figures/rayleigh/"
# # suffix = "_Ra$(abs(Ra))_λ$(λ)_θ$(num.θd)"

# # make_video(num, fwd, gu, "u"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
# # make_video(num, fwd, gv, "v"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
# # make_video(num, fwd, gp, "p"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
# # make_video(num, fwd, gp, "T"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))

# # make_video(num, fwd, gu, "u"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every)
# # make_video(num, fwd, gv, "v"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every)
# # make_video(num, fwd, gp, "p"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every)
# # make_video(num, fwd, gp, "T"; title_prefix=prefix,
# #         title_suffix=suffix, framerate=500÷num.save_every)

# # using Interpolations
# # height = zeros(size(fwd.usave)[1:2])

# # for i in 1:size(fwd.usave,1)
# #     for j in 1:size(fwd.usave,2)
# #         itp = LinearInterpolation(reverse(fwd.usave[i,j,2:end-1]), reverse(gp.x[1,2:end-1]))
# #         height[i,j] = itp(0.0)
# #     end
# # end

# # av_height1 = mean(height, dims=2)[:,1] .+ 0.5 .- H0

# # eRa = av_height1.^3 .* Ra .* (1.0 .- num.θd)

# # fh = Figure(resolution = (1600, 1000))
# # colsize!(fh.layout, 1, Aspect(1, 1.0))
# # ax = Axis(fh[1,1], aspect = 1)
# # lines!(fwd.time, av_height1)
# # # lines!(fwd.time, av_height2)
# # # lines!(fwd.time, av_height3)
# # resize_to_layout!(fh)
