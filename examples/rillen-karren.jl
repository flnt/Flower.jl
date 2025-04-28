using Revise
using Flower
using Random

Random.seed!(1234567)

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)


ratio = 8
L0 = 1.0 
n = 32
nx = ratio*n
ny = n

y = LinRange(-L0/2, L0/2, ny+1)
x = LinRange(-ratio*L0/2, ratio*L0/2, nx+1)

CFL = 0.5
max_it = 4000

num = Numerical(
    case = "Planar",
    x = x,
    y = y,
    CFL = CFL,
    max_iterations = max_it,
    save_every = max_it÷100,
    reinit_every = 1,
    nLS = 1,
    NB = 2,
    nb_reinit = 2
);

gp, gu, gv = init_meshes(num);
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv);

@. gp.LS[1].u = gp.y - 0.4*L0;

function f_rillen(α, κ, x, xrand; L0 = ratio*L0, sigma = sigma, c1 = c1, c2 = c2, limiter = limiter, normalized = normalized)
    β = α + pi/2
    x_shifted = mod.(x - xrand + L0/2, L0) .- L0/2  
    if limiter == true
        if abs(κ) > 1. / num.Δ
            if κ > 0. 
                κ = 1. / num.Δ
            else
                κ = -1. / num.Δ
            end
        end
    end
    if normalized == true
        V = -(1/(sigma*sqrt(2pi)))*exp.(-x_shifted.^2 / (2 * sigma^2)) + c1*κ/sqrt(tan(β)^2 + c2^2)
    else
        V = -(sqrt(pi/2.)*sigma)*exp.(-x_shifted.^2 / (2 * sigma^2)) + c1*κ/sqrt(tan(β)^2 + c2^2)
    end
    return min(0, V)
end

sigma = 0.05
c1 = 0.01
c2 = 1.
limiter = true
normalized = true
speed = 1

@time peaky = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    periodic_x = true,
    time_scheme = CN,
    toy_model = true,
    verbose = true,
    f_interface = f_rillen,
    show_every = 1,
    speed = speed
)

f1 = Figure(size = (1600, 1000))

ax = Axis(f1[1, 1],
aspect = ratio,
xlabel = "x", 
ylabel = "y",
xticks = [-ratio*L0/2, 0, ratio*L0/2],
yticks = [-L0/2, 0, L0/2]
)
xlims!(ax, (-ratio*L0/2, ratio*L0/2))
ylims!(ax, (-L0/2, L0/2))

contour!(ax, gp.x[1,:], gp.y[:,1], fwd.u[1,1,:,:]', levels = 0:0, color=:red, linewidth = 3);
for i = size(fwd.u,2)÷10:size(fwd.u,2)÷10:size(fwd.u,2)-size(fwd.u,2)÷10
    contour!(ax, gp.x[1,:], gp.y[:,1], fwd.u[1,i,:,:]', levels = 0:0, color=:black, linewidth = 3);
end
contour!(ax, gp.x[1,:], gp.y[:,1], fwd.u[1,end,:,:]', levels = 0:0, color=:blue, linewidth = 3);
f1