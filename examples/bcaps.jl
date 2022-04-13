using Revise
using Flower

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

num = Numerical(case = "Sphere",
    L0 = 2.,
    n = 16,
    CFL = 1.0,
    max_iterations = 1,
    R = 0.7)

idx = set_indices(num.n)
tmp, fwd = init_fields(num, idx)
fwd.TL .= 0.

MIXED = run_forward(num, idx, tmp, fwd,
BC_TL = Boundaries(top = Boundary(f = dirichlet, val = 1.0)),
stefan = false,
advection = false,
heat = false,
solid_phase = false,
liquid_phase = false,
verbose = false
)

fa1 = Figure(resolution = (1600, 1000))
colsize!(fa1.layout, 1, Aspect(1, 1.0))
ax = Axis(fa1[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fa1)
fa1 = heatmap!(num.H, num.H, tmp.LIQ[:,:,1]')
fa1 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fa1 = current_figure()

fa2 = Figure(resolution = (1600, 1000))
colsize!(fa2.layout, 1, Aspect(1, 1.0))
ax = Axis(fa2[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fa2)
fa2 = heatmap!(num.H, num.H, tmp.LIQ[:,:,2]')
fa2 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fa2 = current_figure()

fa3 = Figure(resolution = (1600, 1000))
colsize!(fa3.layout, 1, Aspect(1, 1.0))
ax = Axis(fa3[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fa3)
fa3 = heatmap!(num.H, num.H, tmp.LIQ[:,:,3]')
fa3 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fa3 = current_figure()

fa4 = Figure(resolution = (1600, 1000))
colsize!(fa4.layout, 1, Aspect(1, 1.0))
ax = Axis(fa4[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fa4)
fa4 = heatmap!(num.H, num.H, tmp.LIQ[:,:,4]')
fa4 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fa4 = current_figure()

fb1 = Figure(resolution = (1600, 1000))
colsize!(fb1.layout, 1, Aspect(1, 1.0))
ax = Axis(fb1[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fb1)
fb1 = heatmap!(num.H, num.H, tmp.LIQ[:,:,6]')
fb1 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fb1 = current_figure()

fb2 = Figure(resolution = (1600, 1000))
colsize!(fb2.layout, 1, Aspect(1, 1.0))
ax = Axis(fb2[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fb2)
fb2 = heatmap!(num.H, num.H, tmp.LIQ[:,:,7]')
fb2 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fb2 = current_figure()

fw1 = Figure(resolution = (1600, 1000))
colsize!(fw1.layout, 1, Aspect(1, 1.0))
ax = Axis(fw1[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fw1)
fw1 = heatmap!(num.H, num.H, tmp.LIQ[:,:,8]')
fw1 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fw1 = current_figure()

fw2 = Figure(resolution = (1600, 1000))
colsize!(fw2.layout, 1, Aspect(1, 1.0))
ax = Axis(fw2[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fw2)
fw2 = heatmap!(num.H, num.H, tmp.LIQ[:,:,9]')
fw2 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fw2 = current_figure()

fw3 = Figure(resolution = (1600, 1000))
colsize!(fw3.layout, 1, Aspect(1, 1.0))
ax = Axis(fw3[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fw3)
fw3 = heatmap!(num.H, num.H, tmp.LIQ[:,:,10]')
fw3 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fw3 = current_figure()

fw4 = Figure(resolution = (1600, 1000))
colsize!(fw4.layout, 1, Aspect(1, 1.0))
ax = Axis(fw4[1,1], aspect = 1, xticks = -1:1:1, yticks = -1:1:1)
resize_to_layout!(fw4)
fw4 = heatmap!(num.H, num.H, tmp.LIQ[:,:,11]')
fw4 = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
fw4 = current_figure()