using Revise
using Flower

fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)

prefix = "/Users/alex/Documents/PhD/Cutcell/New_ops/robin/adjoint/"
suffix = "_planar_motion"

function J(grid, TDS, TDL, u, x, α=0.)
    @unpack nx, ny = grid
    its, npts = size(TDS)

    id = ny-3
    TS = reshape(TDS[end,1:ny*nx], (ny, nx))
    res = (mean(TS[id,:]) + 0.8) ^ 2

    res += α * x ^ 2

    return res
end

function J_TS(num, grid, TDS, TDL, u, i)
    @unpack nx, ny = grid
    its, npts = size(TDS)
    res = zeros(npts)

    id = ny-3
    s = 0
    if i == num.max_iterations
        TS = reshape(TDS[end,1:ny*nx], (ny, nx))
        s += 2 / nx * (mean(TS[id,:]) + 0.8)
    end

    res[id:ny:ny*nx] .= s

    return res
end

function J_TL(num, grid, TDS, TDL, u, i)
    its, npts = size(TDL)
    zeros(npts)
end

function J_u(num, grid, TDS, TDL, u, i)
    vec(zero(u))
end

function J_x(x, α=0.0)
    2.0 * α * x
end

function R_x(grid, opC_TL)
    @unpack nx, ny = grid
    r2 = zeros(ny*nx)

    for II in grid.ind.b_bottom[1]
        pII = lexicographic(II, ny)
        r2[pII] = -1.0
    end

    return vcat(zeros(ny*nx), opC_TL.χ * r2)
end

function dJ_dx(grid, opC_TL, x, adj)
    its, npts = size(adj.TDS)

    rx = R_x(grid, opC_TL)

    res = J_x(x)
    for i = 2:its
        res -= sum(adj.TDL[i,:] .* rx)
    end

    return res
end

L0 = 2.
n = 10

x = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)
y = LinRange(-L0/2.0 - L0/2.0/n, L0/2.0 + L0/2.0/n, n+1)

num = Numerical(T_inf = 0.0,
    case = "Planar",
    θd = 0.,
    x = x,
    y = y,
    CFL = 0.5,
    shifted = 1e-1,
    max_iterations = 20,
    save_every = 1,
    NB = 1,
    ϵ = 5e-3,
    subdomains = 2,
    overlaps = 1
    )

gp, gu, gv = init_meshes(num)

function planar_motion(num, gp, gu, gv, eps)
    opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)
    
    gp.u .*= -1.
    phL.T .= num.T_inf;
    phS.T .-= num.T_inf;
    
    @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
        opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL,
        phS, phL, fwd, fwdS, fwdL,
        periodic_x = true,
        BC_TL = Boundaries(
            left = Boundary(t = per, f = periodic),
            right = Boundary(t = per, f = periodic),
            bottom = Boundary(t = dir, f = dirichlet, val = 1.0+eps),
        ),
        BC_TS = Boundaries(
            left = Boundary(t = per, f = periodic),
            right = Boundary(t = per, f = periodic),
            top = Boundary(t = dir, f = dirichlet, val = -2.0)
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

    return fwd, fwdS, fwdL, phS, phL, opC_TS, opC_TL
end

fwd0, fwdS0, fwdL0, phS, phL, opC_TS, opC_TL = planar_motion(num, gp, gu, gv, 0.0)
J0 = J(gp, fwdS0.TD, fwdL0.TD, fwd0.u, 1.0, 0.0)

make_video(gp, fwd0.u, fwd0.T; title_prefix=prefix*"T_field", title_suffix=suffix, colormap=Reverse(:ice))

eps_v = [1e-1 / (10 ^ (i/2.0)) for i = 0:5]
Jv = zero(eps_v)
grads = zero(eps_v)
for (i, epsi) in enumerate(eps_v)
    fwd, fwdS, fwdL = planar_motion(num, gp, gu, gv, epsi)[1:3]
    Jv[i] = J(gp, fwdS.TD, fwdL.TD, fwd.u, 1.0 + epsi, 0.0)
    grads[i] = (Jv[i] - J0) / epsi
end

fdf = Figure(resolution = (1600, 1600))
ax = Axis(fdf[1,1], title="Finite differences", xlabel="eps", ylabel="grad"; xscale=log10)#, yscale=log10)
fdf = lines!(eps_v, grads)
fdf = current_figure()

TDS = zeros(num.max_iterations + 1, 2*gp.ny*gp.nx)
TDL = zeros(num.max_iterations + 1, 2*gp.ny*gp.nx)
u = zeros(num.max_iterations + 1, gp.ny*gp.nx)

adj = discrete_adjoint(gp, TDS, TDL, u)

@time MIXED, SOLID, LIQUID = run_backward_discrete(num, gp, gu, gv,
    fwd0, fwdS0, fwdL0, adj, phS, phL,
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
        top = Boundary(t = dir, f = dirichlet, val = -2.0)
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

g_adj = dJ_dx(gp, opC_TL, 1.0, adj)

f = Figure(resolution = (1600, 1600))
ax = Axis(f[1,1])
colsize!(f.layout, 1, Aspect(1, 1))
resize_to_layout!(f)
hidedecorations!(ax)
f = heatmap!(gp.x[1,:], gp.y[:,1], fwdL0.T[end,:,:]', colormap=:ice)
f = contour!(gp.x[1,:], gp.y[:,1], fwd0.u[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
f = current_figure()

fS = Figure(resolution = (1600, 1600))
ax = Axis(fS[1,1])
colsize!(fS.layout, 1, Aspect(1, 1))
resize_to_layout!(fS)
hidedecorations!(ax)
fS = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.TDS[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fS = contour!(gp.x[1,:], gp.y[:,1], fwd0.u[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fS = current_figure();

fL = Figure(resolution = (1600, 1600))
ax = Axis(fL[1,1])
colsize!(fL.layout, 1, Aspect(1, 1))
resize_to_layout!(fL)
hidedecorations!(ax)
fL = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.TDL[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fL = contour!(gp.x[1,:], gp.y[:,1], fwd0.u[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fL = current_figure();

fu = Figure(resolution = (1600, 1600))
ax = Axis(fu[1,1])
colsize!(fu.layout, 1, Aspect(1, 1))
resize_to_layout!(fu)
hidedecorations!(ax)
fu = heatmap!(gp.x[1,:], gp.y[:,1], reshape(veci(adj.u[2,:], gp, 1), (gp.ny, gp.nx))', colormap=:ice)
fu = contour!(gp.x[1,:], gp.y[:,1], fwd0.u[2,:,:]', levels = 0:0, color=:red, linewidth = 3);
limits!(ax, -L0/2., L0/2., -L0/2., L0/2.)
fu = current_figure();

make_video(gp, fwd0.u, adj.TDS; title_prefix=prefix*"adj_TS", title_suffix=suffix, step0=2)
make_video(gp, fwd0.u, adj.TDL; title_prefix=prefix*"adj_TL", title_suffix=suffix, step0=2)
make_video(gp, fwd0.u, adj.u; title_prefix=prefix*"adj_u", title_suffix=suffix, step0=2)

println("Rel error: $(100*abs((grads[4] - g_adj) / grads[4]))%")

nothing