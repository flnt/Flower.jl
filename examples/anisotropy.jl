using Revise
using Flower

L0 = 4.
n = 300

x = LinRange(-L0/2, L0/2, n)
y = LinRange(-L0/2, L0/2, n)

num = Numerical(T_inf = -0.8,
    case = "Crystal",
    θd = 0.,
    ϵ_κ = 0.003,
    ϵ_V = 0.000,
    x = x,
    y = y,
    CFL = 0.5,
    # TEND = 0.045,
    max_iterations = 100,
    aniso = true,
    θ₀ = 0.5*pi,
    m = 6,
    A = -0.2,
    N = 6,
    R = 1.0,
    )

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

phL.T .= num.T_inf;

@time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    stefan = false,
    heat = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    verbose = true,
    advection = false,
    show_every = 1
    );


fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

ax = Axis(f[1,1])#, title = "θ₀ = $(num.θ₀/π)π, mode = $(Int(num.m))")

colsize!(f.layout, 1, Aspect(1, 1))

resize_to_layout!(f)
hidedecorations!(ax)
f = heatmap!(gp.x[1,:], gp.y[:,1], phS.T', colormap=:ice)


# step = num.max_iterations÷20
# #f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

# if step != 0
# for i in 1:step:num.max_iterations
#     f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
# end
# f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
# end


f = current_figure()

#Makie.save("./figures/paper_figures/aniso_theta_pi_4.png", f)

make_video(num, fwd, gp, "T"; title_prefix="/Users/alex/Documents/PhD/Cutcell/New_ops/heat/anisotropy/",
        title_suffix="", framerate=24, minv=-0.8, maxv=0.0)
