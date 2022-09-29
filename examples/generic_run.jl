using Revise
using Flower

ratio = 1
L0 = 2.
nx = 64
ny = ratio * nx

x = LinRange(0., L0/2, nx+1)
y = LinRange(0., ratio*L0/2, ny+1)

num = Numerical(case = "Planar",
    CFL = 1.0,
    TEND = 0.00390625, #0.0009765625,
    # max_iterations = 1,
    x = x,
    y = y,
    R = 0.8,
    save_every = 1,
    T_inf = 0.0,
    NB = 2
)

initial_temp(y, t, lambda, L) = 1 - erf((L - y)/(2*sqrt(t)))/erf(lambda)

initial_pos = 0.125 + eps(1.0)
lambda = 0.9

t = (initial_pos/(2*lambda))^2

gp, gu, gv = init_meshes(num)
opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

@. gp.u = gp.y - L0/2 + initial_pos

H = gp.x[1,:];
for j = 1:size(phL.T, 1)
    for i = 1:size(phL.T, 1)
        y = H[i]
        if y >= (1 - initial_pos)
            phL.T[i,j] = initial_temp(y, t, lambda, 1)
        end
    end
end


@time MIXED, SOLID, LIQUID, radius = run_forward(num, gp, gu, gv,
    opS, opL, phS, phL, fwd,
    BC_TL = Boundaries(
        top = Boundary(t = dir, f = dirichlet, val = 1.0),
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_TS = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    stefan = true,
    advection = true,
    heat = true,
    heat_solid_phase = false,
    heat_liquid_phase = true,
    verbose = false,
    show_every = 1,
    hill = true,
    λ = 1/2.85
)

TLanalytical = similar(phL.T)
TLanalytical .= 0
yfinal = 1 - 2lambda*sqrt(t+num.TEND)

H = gp.x[1,:];
for j = 1:size(phL.T, 1)
    for i = 1:size(phL.T, 1)
        y = H[i]
        if y >= yfinal
            TLanalytical[i,j] = initial_temp(y, (t+num.TEND), lambda, 1)
        end
    end
end

e = phL.T - TLanalytical
@show (radius[end] - yfinal)
@show (normf(e, vcat(LIQUID, MIXED), gp.geoL.cap[:,:,5], num.Δ))

tcks = -ratio*L0/2:2:ratio*L0/2
lim = L0 / 2

fsolid = Figure(resolution = (1600, 1000))  
colsize!(fsolid.layout, 1, Aspect(1, 1.0))
ax = Axis(fsolid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phS.T')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fsolid)

#fsolid = current_figure()

fliquid = Figure(resolution = (1600, 1000))
colsize!(fliquid.layout, 1, Aspect(1, 1.0))
ax = Axis(fliquid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
heatmap!(gp.x[1,:], gp.y[:,1], phL.T')
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# limits!(ax, -lim, lim, -lim, lim)
resize_to_layout!(fliquid)

fliquid = current_figure()