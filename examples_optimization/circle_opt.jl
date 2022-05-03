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
    TEND = 0.1,
    R = 0.75
    )

idx, idxu, idxv = set_indices(num.n)

nprobes = num.n
step = num.n÷(nprobes)
ind = [i*step for i in 1:nprobes]
p1 = [2., 1., 2., 1.]
x_desired = model(num.H[ind], p1)
p2 = [0., 0., 0., 0.]
x_initial = model(num.H[ind], p2)

opt = Optim_parameters(nprobes, ind, idx.b_top[1][ind], [1.0, 1.0, 1e-5, 1.0, 1.0], [p2], [zeros(num.n,num.n)], [zeros(num.n,num.n)], [zeros(num.max_iterations+1, num.n,num.n)])

initial_levelset = @. sqrt(num.X^ 2 + num.Y^ 2) - (num.R)

res, des = gradient_based_optimization(x_desired, x_initial, opt, num, idx, idxu, idxv, initial_levelset, model, opt_iter = 20,
method_opt = LBFGS(linesearch = Optim.LineSearches.BackTracking()))

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end

f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")
lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)
f = current_figure()

#Makie.save("./figures/paper_figures/circle_opt_cost.png", f)

step = num.max_iterations÷10
Iterations = [3, 5, 8, 20]
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

#Makie.save("./figures/paper_figures/circle_opt_noheatmap.png", f)
