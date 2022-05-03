using Revise
using Flower

num = Numerical(T_inf = 0.0,
    case = "Sphere",
    θd = -0.3,
    ϵ_κ = 0.002,
    ϵ_V = 0.000,
    L0 = 4.,
    n = 64,
    CFL = 0.5,
    TEND = 0.09,
    aniso = false,
    θ₀ = 0*pi,
    m = 6,
    A = -0.2,
    N = 6,
    R = 1.0,
    max_iterations = 200
    )

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)


@time MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
    stefan = true,
    heat = true,
    liquid_phase = true,
    solid_phase = true,
    verbose = true,
    advection = true,
    show_every = 20
    );


fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

ax = Axis(f[1,1])

colsize!(f.layout, 1, Aspect(1, 1))

resize_to_layout!(f)
hidedecorations!(ax)
#f = heatmap!(num.H, num.H, (fwd.TL)', colormap=:ice)


step = num.max_iterations÷10
#f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

if step != 0
for i in 1:step:num.max_iterations
    f = contour!(num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
end
f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
end


f = current_figure()


#=

#a = reshape(fwd.Vsave, (num.max_iterations+1, num.n^2))

r = 0.1
a = .2


xc = 0.
yc = 0.

for II in CartesianIndices(fwd.u)
    x = num.X[II]
    y = num.Y[II]
    u1 = -num.R + sqrt((x+xc)^2 + (y+yc)^2)*(1 + num.A*cos(num.N*atan((x+xc)/(y+yc))))
    u2 = -0.5 + sqrt((x-xc)^2 + (y+yc)^2)*(1 + num.A*cos(num.N*atan((x-xc)/(y+yc))))
    u3 = -num.R + sqrt((x+xc)^2 + (y-yc)^2)*(1 + num.A*cos(num.N*atan((x+xc)/(y-yc))))
    u4 = -num.R + sqrt((x-xc)^2 + (y-yc)^2)*(1 + num.A*cos(num.N*atan((x-xc)/(y-yc))))
    u5 = -num.R + sqrt((x)^2 + (y)^2)*(1 + num.A*cos(num.N*atan((x)/(y))))
    fwd.u[II] = u1#min(u1,u2,u3,u4, u5)
    #fwd.u[II] = min(sqrt((x + a)^2 + y^2)-r, sqrt((x - a)^2 + y^2)-r)
    #fwd.u[II] *= ((x - 1)^2 + (y - 1)^2 + 0.1)
end
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f2 = Figure()

colsize!(f2.layout, 1, Aspect(1, 1.0))
ax = Axis(f2[1,1], aspect = 1)  # customized as you see fit

f2 = scatter!(tmp1)
#f2 = scatter!(tmp2)


f2 = current_figure()
=#
#=
time = Makie.Observable(1)



u = @lift(fwd.usave[$time+1,:,:]')
T = @lift(fwd.Tsave[$time+1,:,:]')

fig = heatmap(num.H, num.H, T, colormap= :ice)
contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);

framerate = 30
timestamps = range(1, num.max_iterations, step=2)

record(fig, "figures/movies/exp_crystal.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end


time = Makie.Observable(1)



u = @lift(fwd.usave[$time,:,:]')
T = @lift(fwd.Tsave[$time,:,:]')
u2 = @lift(fwd2.usave[$time,:,:]')
T2 = @lift(fwd2.Tsave[$time,:,:]')
u3 = @lift(fwd3.usave[$time,:,:]')
T3 = @lift(fwd3.Tsave[$time,:,:]')

fig = contour(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
contour!(num.H, num.H, u2, levels = 0:0, color=:red, linewidth = 2)
contour!(num.H, num.H, u3, levels = 0:0, color=:green, linewidth = 2)

framerate = 30
timestamps = range(1, num.max_iterations, step=10)

record(fig, "figures/movies/1crystal.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end=#
