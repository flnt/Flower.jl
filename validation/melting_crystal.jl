using Revise
using Flower

eps_array = (0.005, 0.00)
color_array = ("red", "green")

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)

resize_to_layout!(f)

for i in 1:2
    num = Numerical(case = "Crystal",
        T_inf = -0.5,
        L0 = 2.,
        n = 64,
        CFL = 0.5,
        TEND = 0.2,
        R = 0.5,
        ϵ_κ = eps_array[i],
        N = 4,
        A = 0.2
        )

    idx = set_indices(num.n)
    tmp, fwd = init_fields(num, idx)

    MIXED = run_forward(num, idx, tmp, fwd,
    stefan = true,
    advection = true,
    heat = true,
    solid_phase = false,
    liquid_phase = true,
    verbose = false
    )


    if i == 1
        f = heatmap!(num.H, num.H, (fwd.TL+fwd.TS)', colormap= Reverse(:ice))
        f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:black, linewidth = 3);
    end
    for j = num.max_iterations÷5:num.max_iterations÷5:num.max_iterations
        f = contour!(num.H, num.H, fwd.usave[j,:,:]', levels = 0:0, color=color_array[i], linewidth = 3);
    end
end



f = current_figure()
