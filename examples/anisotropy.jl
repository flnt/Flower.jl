using Revise
using Flower

L0 = 4.
n = 128

x = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)
y = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)

num = Numerical(T_inf = -0.8,
    case = "Crystal",
    θd = 0.,
    ϵ_κ = 0.003,
    ϵ_V = 0.000,
    x = x,
    y = y,
    CFL = 0.5,
    TEND = 0.09,
    aniso = true,
    θ₀ = 0.25*pi,
    m = 6,
    A = -0.2,
    N = 6,
    R = 0.1,
    # max_iterations = 1,
    ϵ = 0.05,
    subdomains = 2,
    overlaps = 1
)

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

phL.T .= num.T_inf;

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd, fwdS, fwdL,
    stefan = true,
    heat = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    verbose = true,
    advection = true,
    show_every = 1
);


fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

ax = Axis(f[1,1])#, title = "θ₀ = $(num.θ₀/π)π, mode = $(Int(num.m))")

colsize!(f.layout, 1, Aspect(1, 1))

resize_to_layout!(f)
hidedecorations!(ax)
f = heatmap!(gp.x[1,:], gp.y[:,1], phL.T', colormap=:ice)


step = num.max_iterations÷20
f = contour!(gp.x[1,:], gp.y[:,1], fwd.u[end-1,:,:]', levels = 0:0, color=:red, linewidth = 3);

if step != 0
for i in 1:step:num.max_iterations
    f = contour!(gp.x[1,:], gp.y[:,1], fwd.u[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
end
f = contour!(gp.x[1,:], gp.y[:,1], fwd.u[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
end

limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)


f = current_figure()

# Makie.save("./figures/paper_figures/aniso_theta_pi_4.png", f)
