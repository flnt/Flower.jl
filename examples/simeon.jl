using Revise
using Flower

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)

L0x = 2.0
L0y = 2.0

n = 129
CFL = 0.5
max_it = 400
A = 0.5
N = 0.5

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
dx = diff(x)[1]
y = collect(-L0y/2:dx:L0y/2)

num = Numerical(
    case = "Planar",
    x = x,
    y = y,
    CFL = CFL,
    max_iterations = max_it,
    save_every = 1,
    reinit_every = 1,
    nLS = 1,
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

# @. gp.LS[1].u = gp.y + A*sin(N*2*pi*gp.x);
@. gp.LS[1].u = gp.y^2+gp.x^2-A^2*(gp.y>=0);
# @. gp.LS[1].u = gp.y*(gp.x^2<A^2)*(gp.y-sqrt(3)*(gp.x+A))*(gp.y+sqrt(3)*(gp.x-A))*(gp.y<sqrt(3)*A);

f1 = Figure(size = (1600, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
f1

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    time_scheme = CN,
    auto_reinit = true,
    verbose = true,
    show_every = 1,
    speed = -1
)

f1 = Figure(size = (1600, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
for i = 1:10:100
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,i,:,:]', levels = 0:0, color=:red, linewidth = 3);
end
f1