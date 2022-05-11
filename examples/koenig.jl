using Revise
using Flower

# Main function solving the level set advection equation
function run_koenig(num, idx, tmp, fwd, s_l, u_x, u_y;
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary())
    )

    @unpack L0, A, N, θd, ϵ_κ, ϵ_V, T_inf, τ, L0, NB, n, Δ, CFL, max_iterations, current_i, reinit_every, nb_reinit, ϵ, H, B, BT, m, θ₀, aniso = num
    @unpack inside, b_left, b_bottom, b_right, b_top = idx
    @unpack SCUT, LCUT, AS, AL, BS, BL, LSA, LSB, SOL, LIQ, sol_projection, liq_projection = tmp
    @unpack iso, u, TS, TL, Tall, V, κ, usave, TSsave, TLsave, Tsave, Vsave, κsave, lengthsave = fwd

    local t = 0.; local faces = zeros(n,n,4);
    local Vsave = zeros(max_iterations+1,n,n); local Vxsave = zeros(max_iterations+1,n,n)
    local Vysave = zeros(max_iterations+1,n,n); local MIXED; local SOLID; local LIQUID;

    usave[1, :, :] .= u[:,:]

    while current_i < max_iterations + 1

        t += τ
        marching_squares!(H, iso, u, TS, TL, κ, SOL, LIQ, sol_projection, liq_projection, Δ, L0, B, BT, inside, ϵ, n, faces)

        bcs!(faces, BC_u.left, Δ); bcs!(faces, BC_u.right, Δ)
        bcs!(faces, BC_u.bottom, Δ); bcs!(faces, BC_u.top, Δ)

        MIXED, SOLID, LIQUID = get_cells_indices(iso, inside)

        Vx = u_x(t)
        Vy = u_y(t)

        V = zeros(n,n)
        V[MIXED] = s_l(t, MIXED)
        velocity_extension!(V, u, vcat(SOLID,LIQUID), n, Δ, NB, BC_u)

        level_update_koenig!(LSA, LSB, u, V, Vx, Vy, inside, CFL, Δ, n)
        u .= reshape(gmres(LSA,(LSB*vec(u))), (n,n))

        FE_reinit(u, Δ, n, nb_reinit, BC_u, idx)

        usave[current_i+1, :, :] .= u[:,:]; Vsave[current_i+1,:,:] .= V
        Vxsave[current_i+1,:,:] .= Vx; Vysave[current_i+1,:,:] .= Vy

        current_i += 1
    end
    return Vsave, Vxsave, Vysave
end

# Numerical setup
num = Numerical(
    L0 = 2., # Box size
    n = 64, # Points per dimension
    TEND = 0.5, # Final time
    R = 0.2 # Initial radius
    )

# Initialize indices and fields
idx = set_indices(num.n)
tmp, fwd = init_fields(num, idx)

# Initialize level set function
@. fwd.u = -num.R + sqrt(num.X^2 + (num.Y+1)^2) *(0.83 - 0.7cos(2*atan(num.X/(num.Y+1))))

# Define time and space dependent velocities
@. s_l(t, ind) = 1*cos(6π*t/num.TEND)*(num.Y[ind] + num.L0/2) # The normal component is only defined in interfacial (ind) cells
@. u_x(t) = 3*cos(4π*t/num.TEND)*(num.Y + num.L0/2) # The x-component of the velocity is defined everywhere
@. u_y(t) = 1*cos(4π*t/num.TEND)*(num.Y + num.L0/2) # The y-component of the velocity is defined everywhere

# Run the simulation
@time Vsave, Vxsave, Vysave = run_koenig(num, idx, tmp, fwd, s_l, u_x, u_y)


# Uncomment for movie

# Vsave .= replace!(Vsave, 0 => NaN);
#
# TIME = Makie.Observable(1)
#
# u = @lift(fwd.usave[$TIME+1,:,:]')
# V = @lift(Vsave[$TIME+1,:,:]')
# Vx = @lift(Vxsave[$TIME+1,:,:]')
# Vy = @lift(Vysave[$TIME+1,:,:]')
#
# fontsize_theme = Theme(fontsize = 15)
#
# # fig = heatmap(num.H, num.H, V)
# # contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
# fig = contour(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
#
# framerate = 24
# TIMEstamps = range(1, num.max_iterations, step=5)
#
# record(fig, "test_koenig.mp4", TIMEstamps;
#         framerate = framerate) do t
#     TIME[] = t
# end


# Uncomment for plot

# fontsize_theme = Theme(fontsize = 30)
# set_theme!(fontsize_theme)
# f = Figure(resolution = (1600, 1600))
# ax = Axis(f[1,1])
# colsize!(f.layout, 1, Aspect(1, 1))
# resize_to_layout!(f)
# hidedecorations!(ax)
#
# step = num.max_iterations÷10
#
# f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);
#
# if step != 0
#     for i in 1:step:num.max_iterations
#         f = contour!(num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
#     end
# end
#
# f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:red, linewidth = 3);
# f = current_figure()
