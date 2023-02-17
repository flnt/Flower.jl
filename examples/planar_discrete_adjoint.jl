using Revise
using Flower

function J_TS(grid, TDS, TDL, u)
    @unpack nx, ny = grid
    its, npts = size(TDS)
    res = zeros(npts)

    id = nx-3
    s = 0
    for i = 1:its
        TS = reshape(TDS[i,1:ny*nx], (ny, nx))
        s += 2 / nx * (mean(TS[id,:]) + 0.8)
    end

    res[id:ny:ny*nx] .= s

    return res
end

function J_TL(grid, TDS, TDL, u)
    its, npts = size(TDL)
    zeros(npts)
end

function J_u(grid, TDS, TDL, u)
    vec(zero(u))
end

L0 = 2.
n = 10

x = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)
y = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)

num = Numerical(T_inf = -1.0,
    case = "Planar",
    θd = 0.,
    x = x,
    y = y,
    CFL = 0.1,
    shifted = 1e-3,
    max_iterations = 10,
    save_every = 1,
    NB = 0,
    ϵ = 1e-8,
    subdomains = 2,
    overlaps = 1
    )

gp, gu, gv = init_meshes(num)
opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd = init_fields(num, gp, gu, gv)

gp.u .*= -1.
phL.T .= num.T_inf;
phS.T .= num.T_inf;

# @profview MIXED, SOLID, LIQUID, A = run_forward(num, gp, gu, gv,
@time MIXED, _, _, SOLID, LIQUID = run_forward(num, gp, gu, gv,
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
    phS, phL, fwd,
    periodic_x = true,
    BC_TL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = dir, f = dirichlet, val = 1.0)
    ),
    BC_TS = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        top = Boundary(t = dir, f = dirichlet, val = -1.0)
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    stefan = true,
    heat = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    verbose = true,
    advection = true,
    show_every = 1
);

TDS = zeros(num.max_iterations + 1, 2*gp.ny*gp.nx)
TDL = zeros(num.max_iterations + 1, 2*gp.ny*gp.nx)
u = zeros(num.max_iterations + 1, gp.ny*gp.nx)

adj = discrete_adjoint(TDS, TDL, u)

@time MIXED, _, _, SOLID, LIQUID = run_backward_discrete(num, gp, gu, gv,
    fwd, adj, phS, phL,
    opC_TS, opC_TL,
    J_TS, J_TL, J_u,
    periodic_x = true,
    BC_TL = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        bottom = Boundary(t = dir, f = dirichlet, val = 1.0)
    ),
    BC_TS = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
        top = Boundary(t = dir, f = dirichlet, val = -1.0)
    ),
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    stefan = true,
    heat = true,
    heat_liquid_phase = true,
    heat_solid_phase = true,
    verbose = true,
    show_every = 1,
    ϵ_adj = 1e-8
);

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

f = Figure(resolution = (1600, 1600))
ax = Axis(f[1,1])
colsize!(f.layout, 1, Aspect(1, 1))
resize_to_layout!(f)
hidedecorations!(ax)
f = heatmap!(gp.x[1,:], gp.y[:,1], fwd.TLsave[end,:,:]', colormap=:ice)
f = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
f = current_figure()

fS = Figure(resolution = (1600, 1600))
ax = Axis(fS[1,1])
colsize!(fS.layout, 1, Aspect(1, 1))
resize_to_layout!(fS)
hidedecorations!(ax)
fS = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.TDS[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fS = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fS = current_figure();

fL = Figure(resolution = (1600, 1600))
ax = Axis(fL[1,1])
colsize!(fL.layout, 1, Aspect(1, 1))
resize_to_layout!(fL)
hidedecorations!(ax)
fL = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.TDL[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fL = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fL = current_figure();

fu = Figure(resolution = (1600, 1600))
ax = Axis(fu[1,1])
colsize!(fu.layout, 1, Aspect(1, 1))
resize_to_layout!(fu)
hidedecorations!(ax)
fu = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.u[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fu = contour!(gp.x[1,:], gp.y[:,1], fwd.usave[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fu = current_figure();