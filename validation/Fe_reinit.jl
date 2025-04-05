using Revise
using Flower

num = Numerical(
    L0 = 4.,
    n = 128
    )

idx, idxu, idxv = set_indices(num.n)
tmp, fwd = init_fields(num, idx, idxu, idxv)

r = 1.
a = 0.7

@inline kink(q,s,d,R) = (q + s)/sqrt((q + s)^2 + d^2) >= abs(q)/abs(R)

BC_u = Boundaries(
    left = Boundary(ind = idx.b_left),
    right = Boundary(ind = idx.b_right),
    bottom = Boundary(ind = idx.b_bottom),
    top = Boundary(ind = idx.b_top))

f = Figure()
fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)

ax1 = Axis(f[1, 1], title = "Iteration 0")
ax2 = Axis(f[1, 2], title = "Iteration 30")
ax3 = Axis(f[1, 3], title = "Iteration 90")
colsize!(f.layout, 1, Aspect(1, 1.0))
colsize!(f.layout, 2, Aspect(1, 1.0))
colsize!(f.layout, 3, Aspect(1, 1.0))
hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax3)

resize_to_layout!(f)
f
for II in CartesianIndices(fwd.u)
    x = num.X[II]
    y = num.Y[II]
    fwd.u[II] = min(sqrt((x + a)^2 + y^2)-r, sqrt((x - a)^2 + y^2)-r)
    fwd.u[II] *= ((x - 1)^2 + (y - 1)^2 + 0.1)
end
contour!(f[1,1],num.H, num.H, fwd.u', levels = -2:.1:1, linewidth = 1.5)
contour!(f[1,1],num.H, num.H, fwd.u', levels = 0:0, linewidth = 3, color=:red)

for II in CartesianIndices(fwd.u)
    x = num.X[II]
    y = num.Y[II]
    fwd.u[II] = min(sqrt((x + a)^2 + y^2)-r, sqrt((x - a)^2 + y^2)-r)
    fwd.u[II] *= ((x - 1)^2 + (y - 1)^2 + 0.1)
end
FE_reinit(fwd.u, num.Δ, num.n, 30, BC_u, idx)

contour!(f[1,2], num.H, num.H, fwd.u', levels = -2:.1:1, linewidth = 1.5)
contour!(f[1,2], num.H, num.H, fwd.u', levels = 0:0, linewidth = 3, color=:red)

for II in CartesianIndices(fwd.u)
    x = num.X[II]
    y = num.Y[II]
    fwd.u[II] = min(sqrt((x + a)^2 + y^2)-r, sqrt((x - a)^2 + y^2)-r)
    fwd.u[II] *= ((x - 1)^2 + (y - 1)^2 + 0.1)
end
FE_reinit(fwd.u, num.Δ, num.n, 90, BC_u, idx)

contour!(f[1,3], num.H, num.H, fwd.u', levels = -2:.1:1, linewidth = 1.5)
contour!(f[1,3], num.H, num.H, fwd.u', levels = 0:0, linewidth = 3, color=:red)

f = current_figure()
