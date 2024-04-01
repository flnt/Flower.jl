using Revise
using Flower
using Peaks

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

prefix = "/Users/alex/Documents/PhD/EnVar/flapping/"

L0 = 10.
n = 128
# x = LinRange(-L0+2.0, 2.0, n+1)
x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

KC = 5.0
# KC = 20.0
β = 1.0
free = false

max_its = 500

num = Numerical{Float64, Int64}(case = "Ellipse",
    Re = β * KC / π,
    CFL = 0.5,
    max_iterations = max_its,
    x = x,
    y = y,
    R = 0.5,
    A = 0.5, # Aspect ratio
    save_every = max_its÷100,
    ϵ = 0.05,
    reinit_every = max_its+1,
)

gp, gu, gv = init_meshes(num);
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv);

phL.u .= 0.0
phL.v .= 0.0

@time xc, yc = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(),
    BC_vL = Boundaries(),
    BC_pL = Boundaries(
        left = Dirichlet(),
        bottom = Dirichlet(),
        right = Dirichlet(),
        top = Dirichlet(),
    ),
    BC_int = [Flapping(KC = KC, β = β, free = free)],
    time_scheme = FE,
    navier_stokes = true,
    ns_advection = true,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1,
    adaptative_t = true
)

case = "nits$(max_its)_KC$(KC)_β$(β)_AR$(num.A)_CFL$(num.CFL)_free$(free)"
file = case * ".jld2"
# save_field(prefix*file, num, gp, phL, fwdL, fwd)

fu = Figure(size = (1000, 600))
ax = Axis(fu[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.ucorr')
cbar = fu[1,2] = Colorbar(fu, hmap, labelpadding = 0)
colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fu)

fv = Figure(size = (1000, 600))
ax = Axis(fv[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.vcorr')
cbar = fv[1,2] = Colorbar(fv, hmap, labelpadding = 0)
colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fv)

pavg = mean(phL.p[gp.LS[1].LIQUID].*num.τ)
pstd = std(phL.p[gp.LS[1].LIQUID].*num.τ)*2

fp = Figure(size = (1000, 600))
ax = Axis(fp[1,1], aspect = DataAspect(), xlabel = L"x", ylabel = L"y")
hmap = heatmap!(gp.x[1,:], gp.y[:,1], (phL.p.*num.τ)', colorrange=(pavg-pstd, pavg+pstd))
cbar = fv[1,2] = Colorbar(fp, hmap, labelpadding = 0)
colsize!(fp.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fp.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fp)

make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field_"*case,
    title_suffix="", framerate=30, colormap = :RdBu)
make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_"*case,
    title_suffix="", framerate=30, colormap = :RdBu)
make_video(num, gp, fwd.u, fwdL.p; title_prefix=prefix*"p_field_"*case,
    title_suffix="", framerate=30, colormap = :RdBu)