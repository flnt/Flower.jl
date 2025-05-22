using Revise
using Flower

using MeshGrid

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 50)
set_theme!(fontsize_theme)

function sdf_line(x, y, x1, y1, x2, y2)
    # Calculate vector from line start to point
    vx = x - x1
    vy = y - y1

    # Calculate vector of line segment
    wx = x2 - x1
    wy = y2 - y1

    # Calculate dot product of vx and wx
    dot_product = vx * wx + vy * wy

    # Calculate length squared of line segment
    length_squared = wx^2 + wy^2

    # Calculate parameter t
    t = clamp(dot_product / length_squared, 0, 1)

    # Calculate closest point on line segment
    closest_x = x1 + t * wx
    closest_y = y1 + t * wy

    # Calculate distance to closest point
    ddx = x - closest_x
    ddy = y - closest_y

    # Return distance
    return sqrt(ddx^2 + ddy^2)
end

n = 64
CFL = 0.5
max_it = 600
K = 3
A = 12
L = 3*A
N = 1
alpha = atan(24*(K*A/2)/(4+(27*(K*A/2)-sqrt(81*(K*A/2)^2+12))*cbrt(sqrt((K*A/2)^2/4+1/27)+(K*A/2)/2)-(27*(K*A/2)+sqrt(81*(K*A/2)^2+12))*cbrt(sqrt((K*A/2)^2/4+1/27)-(K*A/2)/2)))

L0x = 2*(A+1) # 13
L0y = A+L+2 # ceil(A/2/tan(alpha))

peaky0 = A # A/2/tan(alpha)

x = collect(LinRange(-L0x / 2, L0x / 2, n + 1))
dx = diff(x)[1]
y = collect(-(L+1):dx:(A+1)+dx)

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
# @. gp.LS[1].u = -(gp.y^2 + gp.x^2 - A^2*(gp.y>=-dx)); # upper semicircle of radius A
# @. gp.LS[1].u = -(gp.y*(gp.x^2<A^2)*(gp.y-sqrt(3)*(gp.x+A))*(gp.y+sqrt(3)*(gp.x-A))*(gp.y<=sqrt(3)*A)); # ??
# @. gp.LS[1].u = -(gp.y + abs(gp.x)); # ??

# @. gp.LS[1].u = -(sqrt(gp.y^2 + gp.x^2) - A); # circle of radius A
# @. gp.LS[1].u = -((gp.y>0)*(sqrt(gp.y^2 + gp.x^2) - A)+(gp.y<=0)*(abs(gp.x)-A)); # infinite capped rod of radius A
@. gp.LS[1].u = -((gp.y>=0)*(sqrt(gp.y^2 + gp.x^2) - A)-(gp.y<=0)*((gp.y>=-L)*(abs(gp.x)<A)*min(sdf_line(abs(gp.x),gp.y,A,0,A,-L),sdf_line(abs(gp.x),gp.y,0,-L,A,-L))-(1-(gp.y>=-L)*(abs(gp.x)<A))*min(sdf_line(abs(gp.x),gp.y,A,0,A,-L),sdf_line(abs(gp.x),gp.y,0,-L,A,-L)))); # finite capped rod of radius A and height L
# @. gp.LS[1].u = -(-(abs(gp.x)<=A/2)*(gp.y<=tan(alpha/2-pi/4)*(abs(gp.x)-A/2))*gp.y + (gp.y>tan(alpha/2-pi/4)*(abs(gp.x)-A/2))*(gp.y<=tan(alpha)*abs(gp.x)+A/(2*tan(alpha)))*(gp.y>=tan(alpha)*(abs(gp.x)-A/2))*(cos(alpha)*abs(gp.x)+sin(alpha)*gp.y-cos(alpha)*A/2) + (abs(gp.x)>A/2)*(gp.y<tan(alpha)*(abs(gp.x)-A/2))*sqrt((abs(gp.x)-A/2)^2+gp.y^2) + (gp.y>tan(alpha)*abs(gp.x)+A/(2*tan(alpha)))*sqrt(gp.x^2+(gp.y-A/(2*tan(alpha)))^2)); # isosceles triangle of base A and opposite angle 2*alpha
# @. gp.LS[1].u = -((gp.y>=-13)*(gp.y+gp.x^2-6)+(gp.y<-13)*(13-gp.y));


function eta_f(x, k)
    return -1/(24*k).*(4 .+ (27*(k.*x).-sqrt.(81*(k.*x).^2 .+ 12)).*cbrt.(sqrt.((k.*x).^2/4 .+ 1/27).+(k.*x)./2).-(27*(k.*x).+sqrt.(81*(k.*x).^2 .+ 12)).*cbrt.(sqrt.((k.*x).^2/4 .+ 1/27).-(k.*x)/2))
end

function deta_f(x)
    return cbrt.(sqrt.(x.^2/4 .+ 1/27) .- x/2) .- cbrt.(sqrt.(x.^2/4 .+ 1/27) .+ x/2)
end

function vXi(x)
    return 3*(eta_f(x,1) .- x.*deta_f(x))
end

function IC(r)
    phi = sqrt.(A^2 .- r.^2)
    kappa = A^2 ./ sqrt.(A^2 .- r.^2).^3
    return phi,kappa
end


f1 = Figure(size = (1600, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
contour!(gp.x[1,:], gp.y[:,1].-0*peaky0.*ones(length(y)-1), gp.LS[1].u', levels = 0:0, color=:red, linewidth = 3);
# lines!(ax, x, eta_f(x,K), color=:green, linestyle=:dash, linewidth = 3)
f1

function f_interface(α, κ, x, y)
    β = α + pi/2
    # r = 0*sqrt(x^2 + y^2)
    # V = (α<-pi/6)*y
    # if abs(x)>0.01
    # V = (α<0)*min(cos(β)/(1-cos(β)*cbrt(x/sin(β))),1.) # simple dissolution model normal velocity
    if x*β < -0.01
        V = (α<0)*cos(β)/(1-cos(β)*cbrt(x/sin(β))) # simple dissolution model normal velocity
    else
        V = (α<0)/(1-cbrt(cos(β)^2/κ)) # simple dissolution model normal velocity near a pole
    end
    # V = -κ
    return V
end

@time peakx, peaky, tt = run_forward(
    num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL;
    time_scheme = CN,
    toy_model = true,
    verbose = true,
    f_interface = f_interface,
    show_every = 1,
    speed = 1
)


r0 = collect(LinRange(1e-9, A-1e-3, 1001)) #1001))
RR,R0 = meshgrid(r0,r0)
H0, K0 = IC(R0)
HH = H0 .+ eta_f(K0.*RR,1)./K0 .- eta_f(K0.*R0,1)./K0 .- (1 .+ 1 ./ cbrt.(K0)).*(vXi(K0.*R0) .- vXi(K0.*RR))./cbrt.(K0).^2
TT = (1 .+ 1 ./ cbrt.(K0)).^2 .* (vXi(K0.*R0) .- vXi(K0.*RR))./cbrt.(K0).^2



peakx[1], peaky[1] = 0, peaky0

f1 = Figure(size = (1000, 1000))
ax = Axis(f1[1,1], aspect=DataAspect(), xlabel=L"x", ylabel=L"y", xtickalign=0,  ytickalign=0)
contour!(gp.x[1,:], gp.y[:,1] .-0*peaky0.*ones(length(y)-1), fwd.u[1,1,:,:]', levels = 0:0, color=:black, linewidth = 3);
contour!(RR, HH, TT, levels = 0:0, color=:green, linestyle=:dash, linewidth = 3)
for i = 100:100:max_it
    contour!(gp.x[1,:], gp.y[:,1].-0*peaky[i].*ones(length(y)-1), fwd.u[1,i,:,:]', levels = 0:0, color=:red, linewidth = 3);
    contour!(RR, HH, TT, levels = [tt[i]], color=:green, linestyle=:dash, linewidth = 3)
end
contour!(gp.x[1,:], gp.y[:,1] .-0*peaky[end].*ones(length(y)-1), fwd.u[1,end,:,:]', levels = 0:0, color=:black, linewidth = 3);
contour!(RR, HH, TT, levels = [tt[end]], color=:green, linestyle=:dash, linewidth = 3)
# lines!(ax, x, eta_f(x,K), color=:green, linestyle=:dash, linewidth = 3)
ylims!(-L,A)
f1