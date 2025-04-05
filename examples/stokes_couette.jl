using Revise
using Flower

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/navier_slip/"

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 30)
set_theme!(fontsize_theme)

L0 = 2.0
n = 128
d0 = L0 / n

x = LinRange(-L0/2, L0/2, n+1)
y = LinRange(-L0/2, L0/2, n+1)

max_its = 500

num = Numerical(case = "Nothing",
    Re = 1.0,
    CFL = 0.5,
    x = x,
    y = y,
    max_iterations = max_its,
    save_every = max_its÷100,
    ϵ = 0.03,
    nLS = 1,
    nNavier = 1
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

ω = 4.0
R1 = 0.25
R2 = 0.5

phL.u .= 0.0
phL.v .= 0.0

#########################
######## COUETTE ########
#########################
u1 = sqrt.(gp.x .^ 2 + gp.y .^ 2) .- R1 .* ones(gp)
u2 = sqrt.(gp.x .^ 2 + gp.y .^ 2) .- R2 .* ones(gp)
u2 .*= -1.0

r = sqrt.(gp.x .^ 2 + gp.y .^ 2)

λ = zeros(gp)
λ[r .> 0.35] .= 1e10

val = zeros(gp)
val[r .< 0.35] .= ω*R1
val[r .>= 0.35] .= 0.0

#########################
######## CHANNEL ########
#########################
# uPoiseuille = R1^2 .- gu.y[:,1].^2

# u1 = gp.y .+ R1 .+ 1e-8
# u2 = gp.y .- R1 .+ 1e-8
# u2 .*= -1.0

# λ = zeros(gp)
# λ[gp.y .> 0.0] .= 0.0
# λ[gp.y .<= 0.0] .= 1e-5

# val = zeros(gp)
# val[gp.y .> 0.0] .= 0.0
# val[gp.y .<= 0.0] .= 0.0

gp.LS[1].u .= combine_levelsets(gp, u1, u2)

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    BC_uL = Boundaries(
        left = Dirichlet(val = copy(uPoiseuille)),
    ),
    BC_vL = Boundaries(
        left = Dirichlet(),
    ),
    BC_pL = Boundaries(
        right = Dirichlet(),
    ),
    BC_int = [Navier(val = val, λ = λ)],
    time_scheme = FE,
    navier_stokes = true,
    ns_advection = false,
    ns_liquid_phase = true,
    verbose = true,
    show_every = 1
)

fu = Figure(size = (1000, 800))
ax = Axis(fu[1,1], aspect = 1)
hmap = heatmap!(gu.x[1,:], gu.y[:,1], phL.u')
for iLS in 1:num.nLS
    contour!(gu.x[1,:], gu.y[:,1], gu.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3); 
end
cbar = fu[1,2] = Colorbar(fu, hmap)
colsize!(fu.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fu.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fu)

fv = Figure(size = (1000, 800))
ax = Axis(fv[1,1], aspect = 1)
hmap = heatmap!(gv.x[1,:], gv.y[:,1], phL.v')
for iLS in 1:num.nLS
    contour!(gu.x[1,:], gu.y[:,1], gu.LS[iLS].u', levels = 0:0, color=:red, linewidth = 3); 
end
cbar = fv[1,2] = Colorbar(fv, hmap)
colsize!(fv.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fv.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fv)

ft = Figure(size = (1000, 800))
ax = Axis(ft[1,1], aspect = 1)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], reshape(phL.uT, gp)')
cbar = ft[1,2] = Colorbar(ft, hmap)
colsize!(ft.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(ft.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(ft)

# fpr = Figure(size = (1000, 800))
# colsize!(fpr.layout, 1, Aspect(1, 1.0))
# ax  = Axis(fpr[1,1], aspect = 1, xlabel = L"U _ x", ylabel = L"y")
# lines!(phL.u[:,end-10], gu.y[:,1], linewidth = 3)
# limits!(ax, 0.0, 1.1 * maximum(phL.u[:,end-10]), -0.55, 0.55)
# resize_to_layout!(fpr)

suffix = "n$(n)_nits$(num.max_iterations)"

make_video(num, gu, fwd.ux, fwdL.u; title_prefix=prefix*"u_field_",
            title_suffix=suffix, framerate=30)
make_video(num, gv, fwd.uy, fwdL.v; title_prefix=prefix*"v_field_",
        title_suffix=suffix, framerate=30)