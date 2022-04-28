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

res, des = gradient_based_optimization(x_desired, x_initial, opt, num, idx, initial_levelset, model, opt_iter = 20,
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


f = Figure(resolution = (4000, 4000))
step = num.max_iterations÷10
bp = maximum(des.TL)
bm = minimum(des.TL)
Iterations = [3, 5, 7, 17]
x = [1:9, 10:18]; y = [1:8, 9:16]; x_s = [2:9, 11:18]; c = 0;
fontsize_theme = Theme(fontsize = 80)
set_theme!(fontsize_theme)
for i in axes(x,1)
    for j in axes(y,1)
        c += 1
        ax = Axis(f[x_s[i], y[j]])
        hidedecorations!(ax)
        hidespines!(ax)
        heatmap!(f[x_s[i], y[j]], opt.TLsave[Iterations[c]][1:end-1,1:end-1]', colormap=:BuGn_9, colorrange = (bm, bp))
        for ii in 1:step:num.max_iterations
            contour!(f[x_s[i], y[j]], opt.usave[Iterations[c]][ii,:,:]', levels = 0:0, color=:black, linewidth = 5);
        end
        contour!(f[x_s[i], y[j]], des.usave[end,:,:]', levels = 0:0, color=(:blue, 0.7), linewidth = 7);
        contour!(f[x_s[i], y[j]], opt.usave[Iterations[c]][end,:,:]', levels = 0:0, color=:red, linewidth = 7);
        ax2 = Axis(f[x[i][1], y[j]], ylabel = "u", title = @sprintf "Iteration %d" Iterations[c] - 2)
        hidedecorations!(ax2)
        hidespines!(ax2)
        xlims!(-1, 1)
        ylims!(-2, 4.5)
        lines!(f[x[i][1], y[j]], num.H, model(num.H, opt.p[Iterations[c]]), linewidth = 7, color=:red)
        lines!(f[x[i][1], y[j]], num.H, model(num.H, p1), linewidth = 7, color=(:blue, 0.7))
        hlines!(ax2, [0], color = :black, linestyle =:dash, linewidth = 3)
        #Box(f[x[i], y[j]], color = :white, strokewidth = 5)
        Box(f[x[i][1], y[j]], color = :white, strokewidth = 5)
    end
end
Colorbar(f[1:18, 17], limits = (bm, bp), label = "Temperature", colormap=:BuGn_9)
resize_to_layout!(f)
f = current_figure()


#Makie.save("./figures/paper_figures/circle_opt_heatmap_actuator.png", f)

f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1])
lines!(f[1,1], num.H, model(num.H, p1), linewidth = 3, linecolor=:black)
for i in Iterations
    lines!(f[1,1], num.H, model(num.H, opt.p[i]), linewidth = 3)
end
f = current_figure()
