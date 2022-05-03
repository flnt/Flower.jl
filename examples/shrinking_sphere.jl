using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1000))

colsize!(f.layout, 1, Aspect(1, 1.0))
ax = Axis(f[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)  # customized as you see fit

resize_to_layout!(f)

num = Numerical(case = "Sphere",
    L0 = 2.,
    n = 64,
    CFL = 1.0,
    TEND = 0.5,
    R = 0.7)

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)
fwd.TL .= 0.

MIXED = run_forward(num, idx, idxu, idxv, tmp, fwd,
BC_TL = Boundaries(top = Boundary(t = dir, f = dirichlet, val = 1.0)),
stefan = true,
advection = true,
heat = true,
solid_phase = false,
liquid_phase = true,
verbose = false
)

f = heatmap!(num.H, num.H, (fwd.TL+fwd.TS)', colormap= Reverse(:ice))
f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

for j = num.max_iterations÷5:num.max_iterations÷5:num.max_iterations
    f = contour!(num.H, num.H, fwd.usave[j,:,:]', levels = 0:0, color=:black, linewidth = 3);
end

f = current_figure()
