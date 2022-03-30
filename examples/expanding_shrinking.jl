using Revise
using Flower

BLAS.set_num_threads(1)

num = Numerical(R = 0.5,
    case = "Sphere",
    L0 = 2.,
    n = 80,
    CFL = 1.0,
    TEND = 0.1,
    )

idx = set_indices(num.n)

tmp, fwd = init_fields(num, idx)
MIXED = run_forward(num, idx, tmp, fwd,
    advection = true,
    speed = -3
    )

tmp, fwd2 = init_fields(num, idx)
MIXED = run_forward(num, idx, tmp, fwd2,
    advection = true,
    speed = 3
    )


num = Numerical(R = 0.5,
    case = "Crystal",
    L0 = 2.,
    n = 80,
    CFL = 1.0,
    TEND = 0.1,
    A = 0.2,
    N = 4
    )

tmp, fwd3 = init_fields(num, idx)
MIXED = run_forward(num, idx, tmp, fwd3,
    advection = true,
    speed = -3
    )

tmp, fwd4 = init_fields(num, idx)
MIXED = run_forward(num, idx, tmp, fwd4,
    advection = true,
    speed = 3
    )

fontsize_theme = Theme(fontsize = 40)
f = Figure(resolution = (1600, 1000))
fontsize_theme = Theme(fontsize = 40)
f = Figure(Resolution = (1600, 1600))

ax1 = Axis(f[1, 1], title = "Expanding circle")
ax2 = Axis(f[1, 2], title = "Shrinking circle")
ax3 = Axis(f[1, 3], title = "Expanding crystal")
ax4 = Axis(f[1, 4], title = "Shrinking crystal")
colsize!(f.layout, 1, Aspect(1, 1.0))
colsize!(f.layout, 2, Aspect(1, 1.0))
colsize!(f.layout, 3, Aspect(1, 1.0))
colsize!(f.layout, 4, Aspect(1, 1.0))
hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax3)
hidedecorations!(ax4)

resize_to_layout!(f)
contour!(f[1, 1],num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 1.5);

for i = num.max_iterations÷10:num.max_iterations÷10:num.max_iterations
    contour!(f[1, 1], num.H, num.H, fwd2.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 1.5);
end

contour!(f[1, 2],num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 1.5);

for i = num.max_iterations÷10:num.max_iterations÷10:num.max_iterations
    contour!(f[1, 2], num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 1.5);
end

contour!(f[1, 3],num.H, num.H, fwd3.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 1.5);

for i = num.max_iterations÷10:num.max_iterations÷10:num.max_iterations
    contour!(f[1, 3], num.H, num.H, fwd4.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 1.5);
end

contour!(f[1, 4],num.H, num.H, fwd3.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 1.5);

for i = num.max_iterations÷10:num.max_iterations÷10:num.max_iterations
    contour!(f[1, 4], num.H, num.H, fwd3.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 1.5);
end
#f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 2);
f = current_figure()
