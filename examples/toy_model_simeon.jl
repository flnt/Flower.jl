using Revise
using Flower
using Random

Random.seed!(1234)

ratio = 1
L0 = 1.0 #1e-2 #metres
n = 100
nx = ratio*n
ny = n

y = LinRange(-L0/2, L0/2, ny+1)
x = LinRange(-ratio*L0/2, ratio*L0/2, nx+1)

num = Numerical(case = "Mullins",
    x = x,
    y = y,
    max_iterations = 1000,
    shift = 0.0,
    nb_reinit = 10
)

gp, gu, gv = init_meshes(num)
# opS, opL, opC_TS, opC_TL, opC_pS, opC_pL, opC_uS, opC_uL, opC_vS, opC_vL, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv);

gp.u .= gp.y .+ num.shifted #.- 0.15*L0/2*sin.(2pi*gp.x/L0) 
gp.u .-= 0.4*L0 

u_inital = copy(gp.u)
# gp.u .= u_final

function periodic_gaussian(sigma, mu, x; L0 = ratio*L0)
    x_shifted = mod.(x - mu + L0/2, L0) .- L0/2  
    return (1/(sigma*sqrt(2pi)))*exp.(-x_shifted.^2 / (2 * sigma^2))
end

function periodic_gaussian_loop(sigma, mu, x; L0 = ratio*L0)
    gaus = zeros(size(x))
    for i in axes(mu,1)
        gaus += periodic_gaussian.(sigma, mu[i], x)
    end
    return gaus
end


function toy_karren(num, grid;
    periodic_x = false,
    periodic_y = false,
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    amp = 2e-7,
    sigma = 0.5*1e-3,
    nb_drops = 10
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, σ, T_inf, τ, L0, NB, Δ, CFL, Re,
        max_iterations, current_i, save_every, reinit_every, nb_reinit, ϵ, m, θ₀, aniso = num
    @unpack x, y, nx, ny, dx, dy, ind, u, iso, faces, geoS, geoL, V, κ, LSA, LSB = grid

    if periodic_x
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
        BC_u.left.f = BC_u.right.f = periodic
    else
        BC_u.left.ind = ind.b_left;
        BC_u.right.ind = ind.b_right;
    end

    if periodic_y
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
        BC_u.bottom.f = BC_u.top.f = periodic
    else
        BC_u.bottom.ind = ind.b_bottom;
        BC_u.top.ind = ind.b_top;
    end

    # @views fwd.u[1,:,:] .= u

    while current_i < max_iterations + 1

        # gaussian_drop = similar(V)
        # x_drop = rand(x)
        # gaussian_drop = periodic_gaussian.(sigma, x_drop, x);
        
        # x_drop = [rand(x), rand(x)]
        x_drop = [rand(num.x) for i = 1:nb_drops]
        # x_drop = [-0.235]
        # x_drop = [-0.235]

        gaussian_drop = periodic_gaussian_loop(sigma, x_drop, x)

        # @show (gaussian_drop)
        V .= -amp*gaussian_drop * (1 / CFL) * (1/Δ^2)

        # V .= -amp * (1/Δ^2)

        IIOE(grid, LSA, LSB, u, V, CFL, periodic_x, periodic_y)

        u .= reshape(gmres(LSA,(LSB*vec(u))), (ny,nx))

        # u .-= gp.V* Δ^2 * CFL^2
        # if iszero(current_i%save_every) || current_i==max_iterations
        #     snap = current_i÷save_every+1
        #     if current_i==max_iterations
        #         snap = size(fwd.T,1)
        #     end
        #     @views fwd.V[snap,:,:] .= V
        #     @views fwd.u[snap,:,:] .= u
        # end

        FE_reinit(grid, ind, u, nb_reinit, BC_u, periodic_x, periodic_y)

        @show (current_i)

        current_i += 1

    end

    return gp.u
end

@time u = toy_karren(num, gp,
    periodic_x = true,
    BC_u = Boundaries(
        left = Boundary(t = per, f = periodic),
        right = Boundary(t = per, f = periodic),
    ),
    amp = 2e-5,
    sigma = 2e-3,
    nb_drops = 10
    );

tcks_x = -ratio*L0/2:ratio*L0/10:ratio*L0/2
tcks_y = -L0/2:L0/10:L0/2


# contour!(gp.x[1,:], gp.y[:,1], fwd.u[1,:,:]', levels = 0:0, color=:black, linewidth = 1);


fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, ratio))
ax = Axis(fp[1,1], aspect = ratio, xticks = tcks_x, yticks = tcks_y)  # customized as you see fit


contour!(gp.x[1,:], gp.y[:,1], u_inital', levels = 0:0, color=:blue, linewidth = 3);
contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidth = 3);

# limits!(ax, -lim, lim, -lim, lim)
ϵ = 2e-4
# ylims!(ax, 0.5*L0/2-ϵ, 0.5*L0/2+ϵ)
resize_to_layout!(fp)

fp = current_figure()

