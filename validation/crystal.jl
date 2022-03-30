using Revise
using Flower

N_array = (256, 128, 64, 32)
color_array = ("red", "orange", "blue", "green")

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)

Colorbar(f[1, 2], ticks = -0.5:0.1:0, colormap=:ice, limits = (-0.5,0))
resize_to_layout!(f)

for i in 1:4
    num = Numerical(T_inf = -0.7,
        θd = 0.0,
        ϵ_κ = 0.001,
        ϵ_V = 0.000,
        R = 0.3,
        case = "Crystal",
        L0 = 2.,
        n = N_array[i],
        CFL = 0.5,
        TEND = 0.2,
        )

    idx = set_indices(num.n)
    tmp, fwd = init_fields(num, idx)

    fwd.u = -num.R*ones(num.n, num.n) + sqrt.(num.X.^2 + transpose(num.X).^2).*(ones(num.n,num.n) + 0.3*sin.(4*atan.(num.X./(transpose(num.X) + 1E-30*ones(num.n, num.n)))))
    fwd.TL[:,:] .= num.T_inf

    MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
    stefan = true,
    advection = true,
    heat = true,
    solid_phase = false,
    liquid_phase = true,
    verbose = true
    )
    if i == 1
        f = heatmap!(num.H, num.H, (fwd.TL+fwd.TS)', colormap= :ice)
        f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:black, linewidth = 2);
    end

    f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=color_array[i], linewidth = 2);
end

f = current_figure()
