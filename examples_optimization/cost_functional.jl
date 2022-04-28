using Revise
using Flower

@. gradient(field, opt, x) = -(opt.γ[3]*x + field[opt.bc_indices])

@. model(t, p) = p[1]*sin(π*t) + p[2]*cos(π*t) + p[3]*cos(π*t)^2 + p[4]*sin(π*t)^2;

num = Numerical(T_inf = 0.0,
    case = "Sphere",
    ϵ_κ = 0.002,
    θd = 0.,
    L0 = 2.,
    n = 64,
    CFL = 0.5,
    TEND = 0.02,
    R = 0.75
    )

idx = set_indices(num.n)



nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
p1 = [2., 1., 2., 1.]
x_desired = model(num.H[ind], p1)
p2 = [0., 0., 0., 0.]
x_initial = model(num.H[ind], p2)

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-5, 1.0, 1.0], [p2], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = @. sqrt(num.X^ 2 + num.Y^ 2) - (num.R)

res, des = gradient_based_optimization(x_desired, x_initial, opt, num, idx, initial_levelset, model, opt_iter = 1,
method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

rangep = -1:0.5:4
cost = zeros(length(rangep), length(rangep))

for i in enumerate(rangep)
    for j in enumerate(rangep)
        tmp, fwd = init_fields(num, idx)
        p = [i[2], j[2], 2., 1.]
        boundary_values = model(num.H, p)

        MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
            BC_TL = Boundaries(top = Boundary(f = dirichlet, val = boundary_values),
            bottom = Boundary(f = dirichlet, val = boundary_values[end:-1:1]),
            right = Boundary(f = dirichlet, val = boundary_values),
            left = Boundary(f = dirichlet, val = boundary_values[end:-1:1])),
            stefan = true,
            heat = true,
            liquid_phase = true,
            solid_phase = true,
            verbose = true,
            show_every = 50,
            advection = true
            );

        x_boundary = model(num.H[ind], p)
        cost[i[1], j[1]] = cost_functional(fwd, des, opt, idx, num, MIXED)
    end
end

x = y = rangep
z = cost
zmin, zmax = minimum(z), maximum(z)
cmap = :lightrainbow
set_theme!()
fig = Figure(resolution = (1200, 800), fontsize = 22)
ax = Axis3(fig[1, 1], aspect = :data, perspectiveness = 0.5, elevation = π / 7,
    azimuth = π / 4, xlabel = "b₁", ylabel = "a₁", zlabel = "Cost functional",
    xzpanelcolor = (:white, 0.75), yzpanelcolor = (:white, 0.75),
    zgridcolor = :grey, ygridcolor = :grey, xgridcolor = :grey)
sm = surface!(ax, x, y, z; colormap = cmap, colorrange = (zmin, zmax),
    transparency = true)
xm, ym, zm = minimum(ax.finallimits[])
contour!(ax, x, y, z; levels = 20, colormap = cmap, linewidth = 2,
    colorrange = (zmin, zmax), transformation = (:xy, zm),
    transparency = true)
wireframe!(ax, x, y, z; overdraw = true, transparency = true,
    color = (:black, 0.1))
Colorbar(fig[1, 2], sm, height = Relative(0.5))
colsize!(fig.layout, 1, Aspect(1, 1.0))
display(fig)
set_theme!()
