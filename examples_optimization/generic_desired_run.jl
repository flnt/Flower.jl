using Revise
using Flower

num = Numerical(T_inf = 3.0,
    case = "Mullins_cos",
    L0 = 2.,
    n = 50,
    CFL = 0.5,
    TEND = 0.2,
    ϵ_κ = 0.00005,
    A = -0.1,
    N = 2,
    )

idx = set_indices(num.n)

tmp, fwd = init_fields(num, idx)

@. model(t, p) =
    p[1]*cos(π*t) + p[2]*cos(π*t)^2 +
    p[3]*cos(π*t)^3 + p[4]*cos(π*t)^4 - p[5];

@. model_desired(t, p) = p[1]*((1 + cos(pi*t))/2)^4

@. gradient(field, opt, x) = (opt.γ[3]*x - field[opt.bc_indices])

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
x_desired = model_desired(num.H[ind], [10.0, 0.0])
p = 0*ones(8)
x_initial = model_desired(num.H[ind], [0.0, 0.0])

#opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-4, 1.0, 1.0], [p], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = fwd.u
initial_temperature = fwd.TL

#fwd.TS .= -1.0
MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = 0.0)),
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = false,
    advection = true,
    verbose = true,
    show_every = 100
    );

#des = Desired(x_desired, fwd.u, fwd.usave, fwd.TL, fwd.TS)

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

ax = Axis(f[1,1])#, title = "θ₀ = $(num.θ₀/π)π, mode = $(Int(num.m))")  # customized as you see fit

colsize!(f.layout, 1, Aspect(1, 1))

resize_to_layout!(f)
hidedecorations!(ax)
#f = heatmap!(num.H, num.H, (fwd.TL)', colormap=:ice)

step = num.max_iterations÷20
f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

if step != 0
for i in 1:step:num.max_iterations
    f = contour!(num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
end
f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
end


f = current_figure()
