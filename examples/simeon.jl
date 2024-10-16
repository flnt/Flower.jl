using Revise
using Flower

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)

L0x = 2.
L0y = 2.

n = 129
CFL = 0.5
max_it = 100
A = 0.5
N = 1
alpha = pi/9

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
dx = diff(x)[1]
y = collect(-L0y/2:dx:L0y/2+dx)

num = Numerical(
    case = "Planar",
    x = x,
    y = y,
    CFL = CFL,
    max_iterations = max_it,
    save_every = 1,
    reinit_every = 1,
    nLS = 1,
    NB = 2,
    nb_reinit = 2
)

gp, gu, gv = init_meshes(num)
op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)

# @. gp.LS[1].u = -gp.y;
# @. gp.LS[1].u = -(gp.y + A*sin(N*2*pi*gp.x));
@. gp.LS[1].u = -(sqrt(gp.y^2 + gp.x^2) - A); # circle of radius A
# @. gp.LS[1].u = -(gp.y^2 + gp.x^2 - A^2*(gp.y>=-dx)); # upper semicircle of radius A
# @. gp.LS[1].u = -(gp.y*(gp.x^2<A^2)*(gp.y-sqrt(3)*(gp.x+A))*(gp.y+sqrt(3)*(gp.x-A))*(gp.y<=sqrt(3)*A)); # ??
# @. gp.LS[1].u = -(gp.y + abs(gp.x)); # ??
# @. gp.LS[1].u = -(-(abs(gp.x)<=A/2)*(gp.y<=tan(alpha/2-pi/4)*(abs(gp.x)-A/2))*gp.y + (gp.y>tan(alpha/2-pi/4)*(abs(gp.x)-A/2))*(gp.y<=tan(alpha)*abs(gp.x)+A/(2*tan(alpha)))*(gp.y>=tan(alpha)*(abs(gp.x)-A/2))*(cos(alpha)*abs(gp.x)+sin(alpha)*gp.y-cos(alpha)*A/2) + (abs(gp.x)>A/2)*(gp.y<tan(alpha)*(abs(gp.x)-A/2))*sqrt((abs(gp.x)-A/2)^2+gp.y^2) + (gp.y>tan(alpha)*abs(gp.x)+A/(2*tan(alpha)))*sqrt(gp.x^2+(gp.y-A/(2*tan(alpha)))^2)); # isosceles triangle of base A and opposite angle 2*alpha

f1 = Figure(size = (1600, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
contour!(gp.x[1,:], gp.y[:,1], gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
f1

function f_interface(α, x, y)
    β = α + pi/2
    r = 0*sqrt(x^2 + y^2)
    V = cos(β)/(1. - cos.(β)*cbrt(r)/(cbrt(sin(β)) + 1e-12)) # simple dissolution model normal velocity
    return V
end

@time run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    time_scheme = CN,
    toy_model = true,
    verbose = true,
    f_interface = f_interface,
    show_every = 1,
    speed = 10
)

f1 = Figure(size = (1600, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 3);
for i = 10:10:max_it
    contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,i,:,:]', levels = 0:0, color=:red, linewidth = 3);
end
f1