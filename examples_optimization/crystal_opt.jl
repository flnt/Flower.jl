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

idx, idxu, idxv = set_indices(num.n)

tmp, fwd = init_fields(num, idx, idxu, idxv)

@. model(t, p) =
    p[1]*cos(num.N*pi*t/8) + p[2]*cos(num.N*pi*t/8)^2 +
    p[3]*cos(num.N*pi*t/8)^3 + p[4]*cos(num.N*pi*t/8)^4 +
    p[5]*sin(num.N*pi*t/8) + p[6]*sin(num.N*pi*t/8)^2 +
    p[7]*sin(num.N*pi*t/8)^3 + p[8]*sin(num.N*pi*t/8)^4 ;

@. model2(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5) + p[2]*(cos(num.N*pi*t/8) - 0.5);

@. model_desired(t, p) = p[1]*((1 + cos(num.N*pi*t/8))/2)^4

@. gradient(field, opt, x) = (opt.γ[3]*x - field[opt.bc_indices])

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
x_desired = model_desired(num.H[ind], [10., 0.])
p = 0*ones(8)
x_initial = model_desired(num.H[ind], [0., 0.])

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-4, 1.0, 1.0], [p], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = fwd.u
initial_temperature = fwd.TL

MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = x_desired),
    bottom = Boundary(f = neumann, val = x_desired[end:-1:1]),
    right = Boundary(f = neumann, val = x_desired),
    left = Boundary(f = neumann, val = x_desired[end:-1:1])),
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = true,
    advection = true,
    verbose = true,
    );

des = Desired(x_desired, fwd.u, fwd.usave, fwd.TL, fwd.TS)

function fg2!(F, G, x, des, opt, num, idx, idxu, idxv, initial_levelset, initial_temperature, basis)

  tmp, fwd = init_fields(num, idx, idxu, idxv)
  xdata = num.H[opt.ind];
  p = curve_fit(basis, xdata, x, zeros(8))
  @show (p.param)
  boundary_values = basis(num.H, p.param)
  MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
      BC_TL = Boundaries(top = Boundary(f = neumann, val = boundary_values),
      bottom = Boundary(f = neumann, val = boundary_values[end:-1:1]),
      right = Boundary(f = neumann, val = boundary_values),
      left = Boundary(f = neumann, val = boundary_values[end:-1:1])),
      stefan = true,
      heat = true,
      liquid_phase = true,
      solid_phase = true,
      verbose = false,
      advection = true
      );
  s = similar(fwd.u)
  adj = my_Adjoint(s, fwd.u, opt.γ[1].*(des.TS - fwd.TS), opt.γ[1].*(des.TL - fwd.TL), s, s, s, s)
  tmp, fwd = init_fields(num, idx, idxu, idxv)
  run_backward(num, idx, tmp, fwd, adj,
      stefan = true,
      heat = true,
      liquid_phase = true,
      solid_phase = true,
      verbose = false,
      advection = true
      );
  if G != nothing
      G .= gradient(adj.TL, opt, x)
      push!(opt.p, p.param)
      push!(opt.TLsave, fwd.TL)
      push!(opt.TSsave, fwd.TS)
      push!(opt.usave, fwd.usave)
  end
  if F != nothing
    value = cost_functional(fwd, des, opt, idx, num, MIXED)
    @show (value)
    return value
  end
end

function gradient_based_optimization2(x_desired, x_initial, opt, num, idx, idxu, idxv, des, initial_levelset, initial_temperature, basis;
    method_opt = ConjugateGradient(),
    opt_iter = 10)

    res = optimize(Optim.only_fg!((F, G, x)->fg2!(F, G, x, des, opt, num, idx, idxu, idxv, initial_levelset, initial_temperature, basis)), x_initial, method_opt,
    Optim.Options(store_trace = true, show_trace=true, iterations = opt_iter, allow_f_increases = false))

    @show Optim.minimizer(res)

    return res
end

<<<<<<< HEAD
res = gradient_based_optimization2(x_desired, x_initial, opt, num, idx, des, initial_levelset, initial_temperature, model_desired,
    opt_iter = 20,
=======
res = gradient_based_optimization2(x_desired, x_initial, opt, num, idx, idxu, idxv, des, initial_levelset, initial_temperature, model2,
    opt_iter = 5,
>>>>>>> rayleigh_benard
    method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end

# df = DataFrame(iteration = [res.trace[i].iteration for i in 1:length(res.trace)],
#     value = [res.trace[i].value for i in 1:length(res.trace)],
#     g_norm = [res.trace[i].g_norm for i in 1:length(res.trace)],
#     p = opt.p[2:end]);
# CSV.write("examples_optimization/data/opt_crystal.csv", df);
# df2 = DataFrame(res = res);
# CSV.write("examples_optimization/data/res_crystal.csv", df2);


let c = 0
    f = Figure(resolution = (4000, 4000))
    step = num.max_iterations÷10
    Iterations = [2, 17, 19, 22]
    bp = maximum(abs.(des.TL - opt.TLsave[2]))
    bm = minimum(abs.(des.TL - opt.TLsave[2]))
    x = [1:9, 10:18]; y = [1:8, 9:16]; x_s = [2:9, 11:18];
    fontsize_theme = Theme(fontsize = 80)
    set_theme!(fontsize_theme)
    for i in axes(x,1)
        for j in axes(y,1)
            c += 1
            ax = Axis(f[x_s[i], y[j]])
            hidedecorations!(ax)
            hidespines!(ax)
            heatmap!(f[x_s[i], y[j]], abs.(des.TL[1:end-1,1:end-1] - opt.TLsave[Iterations[c]][1:end-1,1:end-1])', colormap=:BuGn_9, colorrange = (bm, bp))
            for ii in 1:step:num.max_iterations
                contour!(f[x_s[i], y[j]], opt.usave[Iterations[c]][ii,:,:]', levels = 0:0, color=:black, linewidth = 5);
            end
            contour!(f[x_s[i], y[j]], des.usave[end,:,:]', levels = 0:0, color=(:blue, 0.7), linewidth = 7);
            contour!(f[x_s[i], y[j]], opt.usave[Iterations[c]][end,:,:]', levels = 0:0, color=:red, linewidth = 7);
            ax2 = Axis(f[x[i][1], y[j]], ylabel = "u", title = @sprintf "Iteration %d" Iterations[c] - 2)
            hidedecorations!(ax2)
            hidespines!(ax2)
            #xlims!(-1, 1)
            ylims!(-1, 15)
            lines!(f[x[i][1], y[j]], num.H, model_desired(num.H, opt.p[Iterations[c]]), linewidth = 7, color=:red)
            lines!(f[x[i][1], y[j]], num.H, x_desired, linewidth = 7, color=(:blue, 0.7))
            hlines!(ax2, [0], color = :black, linestyle =:dash, linewidth = 3)
            Box(f[x[i], y[j]], color = :white, strokewidth = 5)
            Box(f[x[i][1], y[j]], color = :white, strokewidth = 5)
        end
    end
    Colorbar(f[1:18, 17], limits = (bm-0.05, bp), label = "Temperature error", colormap=:BuGn_9)
    resize_to_layout!(f)
    f = current_figure()
    Makie.save("./figures/paper_figures/crystal_opt_heatmap_actuator_new.png", f)
end


f = Figure()
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")

lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)

f = current_figure()

# Makie.save("./figures/paper_figures/crystal_opt_cost.png", f)
