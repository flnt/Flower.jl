using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ratio = 1
L0 = 1.
nx = 64
ny = ratio * 64

x = LinRange(-L0/2, L0/2, nx+1)
y = LinRange(-ratio*L0/2, ratio*L0/2, ny+1)

num = Numerical(case = "Mullins",
    Re = 1.0,
    CFL = 0.1,
    TEND = 1.0,
    x = x,
    y = y,
    R = 0.5,
    u_inf = 0.0,
    save_every = 5,
    shifted = 0.000,
    N = 1,
    max_iterations = 2000,
    A = 0.0,
    ϵ_κ = 0.00000,
    θd = 0.0,
    NB = 2,
    reinit_every = 1,
    ϵ = 0.05
)

Ra = 1e5
λ = 1.0
H0 = 0.05
T1 = 1.0
T2 = 0.0

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@. gp.u = -gp.x - L0/2 + H0 + 0.0001

# for II in CartesianIndices(gp.x)
#     y = gp.x[II] + L0/2
#     if (y <= H0)
#         phL.T[II] = 1 +(num.θd-1)*y/H0;
#     elseif (y > H0)
#         phS.T[II] = num.θd*(y-1)/(H0-1)
#     end
# end


@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    periodic_y = true,
    BC_TL = Boundaries(
        left = Boundary(t = dir, f = dirichlet, val = T1),
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
    ),
    BC_TS = Boundaries(
        right = Boundary(t = dir, f = dirichlet, val = T2),
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
    ),
    BC_u = Boundaries(
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
    ),
    BC_uL = Boundaries(
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
        left = Boundary(t = dir, f = dirichlet, val = 0.0),
        right = Boundary(t = dir, f = dirichlet, val = 0.0),
    ),
    BC_vL = Boundaries(
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
        left = Boundary(t = dir, f = dirichlet, val = 0.0),
        right = Boundary(t = dir, f = dirichlet, val = 0.0),
    ),
    BC_pL = Boundaries(
        top = Boundary(t = per, f = periodic),
        bottom = Boundary(t = per, f = periodic),
    ),
    stefan = true,
    advection = true,

    heat = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    
    heat_convection = true,

    navier_stokes = true,
    ns_advection = false,
    ns_liquid_phase = true,
    ns_solid_phase = false,

    verbose = true,
    show_every = 1,

    adaptative_CFL = true,

    Ra = Ra,
    λ = λ,
)

tcks = -ratio*L0/2:2:ratio*L0/2
lim = L0 / 2

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, 1.0))
ax = Axis(fp[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.T', colorrange=(num.θd, T1))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fp)

fp = current_figure()

fu = Figure(resolution = (1600, 1000))
colsize!(fu.layout, 1, Aspect(1, 1.0))
ax = Axis(fu[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gu.x[1,:], gu.y[:,1], phL.u', colorrange=(minimum(phL.u), maximum(phL.u)))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fu)

fu = current_figure()

fv = Figure(resolution = (1600, 1000))
colsize!(fv.layout, 1, Aspect(1, 1.0))
ax = Axis(fv[1,1], aspect = 1, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gv.x[1,:], gv.y[:,1], phL.v', colorrange=(minimum(phL.v), maximum(phL.v)))
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fv)

fv = current_figure()


prefix = "/Users/alex/codes/Flower.jl/figures/rayleigh/"
suffix = "_Ra$(abs(Ra))_λ$(λ)"

make_video(num, fwd, gu, "u"; title_prefix=prefix,
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
make_video(num, fwd, gv, "v"; title_prefix=prefix,
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
make_video(num, fwd, gp, "p"; title_prefix=prefix,
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
make_video(num, fwd, gp, "T"; title_prefix=prefix,
        title_suffix=suffix, framerate=500÷num.save_every, limitsx=(-lim,lim), limitsy=(-lim,lim))
# make_video(num, fwd, gu, "ucorr"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, limitsx=(-lim,lim), limitsy=(-lim,lim))
# make_video(num, fwd, gp, "ϕ"; title_prefix=prefix,
#         title_suffix=suffix, framerate=100, limitsx=(-lim,lim), limitsy=(-lim,lim))


using Interpolations
height = zeros(size(fwd.usave)[1:2])

for i in 1:size(fwd.usave,1)
    for j in 1:size(fwd.usave,2)
        itp = LinearInterpolation(reverse(fwd.usave[i,j,2:end-1]), reverse(gp.x[1,2:end-1]))
        height[i,j] = itp(0.0)
    end
end

av_height1 = mean(height, dims=2)[:,1] .+ 0.5 .- H0

eRa = av_height1.^3 .* Ra .* (1.0 .- num.θd)

fh = Figure(resolution = (1600, 1000))
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1)
lines!(fwd.time, av_height1)
lines!(fwd.time, av_height2)
# lines!(fwd.time, av_height3)
resize_to_layout!(fh)
