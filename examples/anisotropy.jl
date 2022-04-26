using Revise
using Flower

num = Numerical(T_inf = -0.8,
    case = "Crystal",
    θd = 0.,
    ϵ_κ = 0.003,
    ϵ_V = 0.000,
    L0 = 4.,
    n = 300,
    CFL = 0.5,
    TEND = 0.09,
    aniso = true,
    θ₀ = 0.25*pi,
    m = 6,
    A = -0.2,
    N = 6,
    R = 0.1,
    #max_iterations = 1
    )

idx = set_indices(num.n)
tmp, fwd = init_fields(num, idx)

fwd.TL .= num.T_inf;

@time MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = true,
    verbose = true,
    advection = true,
    show_every = 50
    );


fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

ax = Axis(f[1,1])#, title = "θ₀ = $(num.θ₀/π)π, mode = $(Int(num.m))")

colsize!(f.layout, 1, Aspect(1, 1))

resize_to_layout!(f)
hidedecorations!(ax)
f = heatmap!(num.H, num.H, (fwd.TL)', colormap=:ice)


step = num.max_iterations÷20
#f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

if step != 0
for i in 1:step:num.max_iterations
    f = contour!(num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
end
f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
end


f = current_figure()

#Makie.save("./figures/paper_figures/aniso_theta_pi_4.png", f)
