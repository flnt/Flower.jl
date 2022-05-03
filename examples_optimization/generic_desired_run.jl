using Revise
using Flower

BLAS.set_num_threads(1)


num = Numerical(T_inf = -0.6,
    case = "3Crystals",
    θd = 0.,
    ϵ_κ = 0.0005,
    ϵ_V = 0.0000,
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

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

x_c = 0.2
y_c = 0.2

@. model(t, p) =
    p[1]*sin(π*t) + p[2]*cos(π*t) +
    p[3]*cos(π*t)^2 + p[4]*sin(π*t)^2 +
    p[5]*cos(π*t)^3 + p[6]*sin(π*t)^3 +
    p[7]*cos(π*t)^4 + p[8]*sin(π*t)^4;

@. model_desired(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5);

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]

p = [-0.027580156399409195,
 -40.408381949607566,
  40.57605950959059,
  14.710086551973449]

boundary_values = model_desired(num.H[ind], [0.0])

f = Figure()
ax = Axis(f[1,1])
f = lines!(boundary_values)
f = current_figure()

@time MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = boundary_values),
    bottom = Boundary(f = neumann, val = boundary_values[end:-1:1]),
    right = Boundary(f = neumann, val = boundary_values),
    left = Boundary(f = neumann, val = boundary_values[end:-1:1])),
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = true,
    verbose = true,
    advection = true,
    show_every = 100,
    #periodic_x = true,
    #periodic_y = true
    );


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
