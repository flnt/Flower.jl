using Revise
using Flower

num = Numerical(T_inf = 1.2,
    case = "Mullins_cos",
    L0 = 2.,
    n = 64,
    CFL = 0.5,
    TEND = 0.5,
    A = -0.05,
    N = 2
    )

idx, idxu, idxv = set_indices(num.n)

tmp, fwd = init_fields(num, idx, idxu, idxv)

@. model(t, p) =
    p[1]*sin(2π*t) + p[2]*cos(2π*t) +
    p[3]*cos(2π*t)^2 + p[4]*sin(2π*t)^2 +
    p[5]*cos(2π*t)^3 + p[6]*sin(2π*t)^3 +
    p[7]*cos(2π*t)^4 + p[8]*sin(2π*t)^4;

@. model_desired(t, p) = p[1]*((1 + cos(num.N*pi*t))/2)^4;

@. gradient(field, opt, x) = -(opt.γ[3]*x + field[opt.bc_indices])

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
x_desired = model_desired(num.H[ind], [10.0])
p = 0*ones(8)
x_initial = model_desired(num.H[ind], p)

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-3, 1.0, 1.0], [p], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = fwd.u
initial_temperature = fwd.TL

MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
    BC_TL = Boundaries(top = Boundary(f = neumann, val = x_desired)),
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
      BC_TL = Boundaries(top = Boundary(f = neumann, val = boundary_values)),
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
    Optim.Options(store_trace = true, show_trace=true, iterations = opt_iter, allow_f_increases = true))

    @show Optim.minimizer(res)

    return res
end

res = gradient_based_optimization2(x_desired, x_initial, opt, num, idx, idxu, idxv, des, initial_levelset, initial_temperature, model_desired,
    opt_iter = 25,
    method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end


step = num.max_iterations÷10
Iterations = [2, 5, 8, 12]
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
    ax = Axis(f[ix, j], title = @sprintf "Iteration %d" Iterations[i])
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

#Makie.save("./figures/paper_figures/mullins_opt_noheatmap.png", f)


f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")

lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)

f = current_figure()

#Makie.save("./figures/paper_figures/mullins_opt_cost.png", f)
