using Revise
using Flower

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

idx = set_indices(num.n)

tmp, fwd = init_fields(num, idx)

@. model(t, p) =
    p[1]*cos(0.5π*t) +
    p[2]*cos(0.5π*t)^2 + p[3]*sin(0.5π*t)^2;
    
@. model2(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5) + p[2]*(cos(num.N*pi*t/8) - 0.5);

@. model_desired(t, p) = p[1]*(cos(num.N*pi*t/8) - 0.5);

@. gradient(field, opt, x) = -(opt.γ[3]*x + field[opt.bc_indices])

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
x_desired = model_desired(num.H[ind], [-1.5])
p = 0*ones(8)
x_initial = model_desired(num.H[ind], [0.5])

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-4, 1.0, 1.0], [p], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

global x_c = 0.2
global y_c = 0.2

initial_levelset = fwd.u
initial_temperature = fwd.TL

MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
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

function fg2!(F, G, x, des, opt, num, idx, initial_levelset, initial_temperature, basis)

  tmp, fwd = init_fields(num, idx)
  xdata = num.H[opt.ind];
  p = curve_fit(basis, xdata, x, zeros(8))
  @show (p.param)
  boundary_values = basis(num.H, p.param)
  MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
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
  tmp, fwd_ = init_fields(num, idx)
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

function gradient_based_optimization2(x_desired, x_initial, opt, num, idx, des, initial_levelset, initial_temperature, basis;
    method_opt = ConjugateGradient(),
    opt_iter = 10)

    res = optimize(Optim.only_fg!((F, G, x)->fg2!(F, G, x, des, opt, num, idx, initial_levelset, initial_temperature, basis)), x_initial, method_opt,
    Optim.Options(store_trace = true, show_trace=true, iterations = opt_iter, allow_f_increases = false))

    @show Optim.minimizer(res)

    return res
end

res = gradient_based_optimization2(x_desired, x_initial, opt, num, idx, des, initial_levelset, initial_temperature, model2,
    opt_iter = 5,
    method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end


step = num.max_iterations÷10
Iterations = [2, 3, 4, 5]
bp = findmax(des.TL - opt.TLsave[Iterations[1]])[1]
bm = findmin(des.TL - opt.TLsave[Iterations[1]])[1]
f = Figure(resolution = (4000, 4000))
fontsize_theme = Theme(fontsize = 80)
set_theme!(fontsize_theme)

for i in 1:4
    j = i
    ix = 1
    if i > 2
        ix = 2
        j = i-2
    end
    ax = Axis(f[ix, j], title = @sprintf "Iteration %d" Iterations[i])#$(N_array[i])
    colsize!(f.layout, j, Aspect(1, 1.0))
    hidedecorations!(ax)
    #heatmap!(f[ix, j], (des.TL - opt.TLsave[Iterations[i]])', colorrange = (bm, bp))
    for ii in 1:step:num.max_iterations
        contour!(f[ix, j], opt.usave[Iterations[i]][ii,:,:]', levels = 0:0, color=:black, linewidth = 5);
    end
    contour!(f[ix, j], des.usave[end,:,:]', levels = 0:0, color=(:red, 0.7), linewidth = 7);
    contour!(f[ix, j], opt.usave[Iterations[i]][end,:,:]', levels = 0:0, color=:blue, linewidth = 7);
end
#Colorbar(f[1:2, 3], limits = (bm, bp), label = "Error")
f = current_figure()

#Makie.save("./figures/paper_figures/crystal_opt_noheatmap.png", f)


f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")

lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)

f = current_figure()

#Makie.save("./figures/paper_figures/crystal_opt_cost.png", f)
