using Revise
using Flower

num = Numerical(T_inf = -0.6,
    case = "3Crystals",
    θd = 0.,
    ϵ_κ = 0.0005,
    ϵ_V = 0.002,
    L0 = 4.,
    n = 100,
    CFL = 0.5,
    TEND = 0.45,
    aniso = true,
    θ₀ = 0.5*pi,
    m = 4,
    A = 0.2,
    N = 4,
    R = 0.14,
    #max_iterations = 1
    )

<<<<<<< HEAD
idx = set_indices(num.n)

tmp, fwd = init_fields(num, idx)
=======
idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)
>>>>>>> rayleigh_benard

@. model(t, p) = p[1]*sin(num.N*pi*t/8) + p[2]*cos(num.N*pi*t/8) + p[3]*cos(num.N*pi*t/8)^2 + p[4]*sin(num.N*pi*t/8)^2;

@. model2(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5) + p[2]*(cos(num.N*pi*t/8) - 0.5);

# @. model_desired(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5);
@. model_desired(t, p) = p[1]*((1 + cos(num.N*pi*t/8))/2)^4 + p[2]

@. gradient(field, opt, x) = -(opt.γ[3]*x + field[opt.bc_indices])

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
x_desired = model_desired(num.H[ind], [10., 0.])
p = 0*ones(8)
x_initial = model_desired(num.H[ind], [0., 0.])

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 10.0, 1e-2, 1.0, 1.0], [p], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = fwd.u
initial_temperature = fwd.TL

<<<<<<< HEAD
MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = x_desired),
    bottom = Boundary(f = neumann, val = x_desired[end:-1:1]),
    right = Boundary(f = neumann, val = x_desired),
    left = Boundary(f = neumann, val = x_desired[end:-1:1])),
=======
f = Figure()
ax = Axis(f[1,1])
f = lines!(boundary_values)
f = current_figure()

@time MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = boundary_values),
    bottom = Boundary(f = neumann, val = boundary_values[end:-1:1]),
    right = Boundary(f = neumann, val = boundary_values),
    left = Boundary(f = neumann, val = boundary_values[end:-1:1])),
>>>>>>> rayleigh_benard
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = true,
    advection = true,
    verbose = true,
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
