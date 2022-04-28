using Revise
using Flower

num = Numerical(
    case = "Sphere",
    L0 = 2.,
    n = 64,
    TEND = 1.0,
    max_iterations = 10,
    R = 0.4
    )


idx = set_indices(num.n)
tmp, fwd = init_fields(num, idx)

#fwd.u .= sqrt.(num.X .^ 2 + (num.Y.+0.7) .^ 2) - (num.R) * ones(num.n, num.n);

s_l = -4.0
u_x = 1.0
u_y = 1.0

function run_koening(num, idx, tmp, fwd, s_l::Number, u_x::Number, u_y::Number;
    BC_u = Boundaries(
        left = Boundary(),
        right = Boundary(),
        bottom = Boundary(),
        top = Boundary()),
    verbose = false,
    show_every = 100
    )

    @unpack L0, A, τ, n, Δ, CFL, max_iterations, current_i, nb_reinit = num
    @unpack inside, b_left, b_bottom, b_right, b_top = idx
    @unpack LSA, LSB = tmp
    @unpack iso, u, usave = fwd

    usave[1, :, :] .= u[:,:]

    while current_i < max_iterations + 1

        V = s_l*ones(n,n)
        Vx = u_x*ones(n,n)
        Vy = u_y*ones(n,n)

        @sync begin
            @spawn bcs!(u, BC_u.left, Δ)
            @spawn bcs!(u, BC_u.right, Δ)
            @spawn bcs!(u, BC_u.bottom, Δ)
            @spawn bcs!(u, BC_u.top, Δ)
        end

        Flower.level_update_koening!(LSA, LSB, u, V, Vx, Vy, inside, CFL, Δ, n)
        try
            u .= reshape(gmres(LSA,(LSB*vec(u))), (n,n))
        catch
            @error ("Inadequate level set function, iteration $current_i")
            break
        end
        usave[current_i+1, :, :] .= u[:,:]

        current_i += 1
    end
end

run_koening(num, idx, tmp, fwd, s_l, u_x, u_y)

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
f = Figure(resolution = (1600, 1600))
ax = Axis(f[1,1])
colsize!(f.layout, 1, Aspect(1, 1))
resize_to_layout!(f)
hidedecorations!(ax)

step = num.max_iterations÷10

f = contour!(num.H, num.H, fwd.usave[1,:,:]', levels = 0:0, color=:red, linewidth = 3);

if step != 0
    for i in 1:step:num.max_iterations
        f = contour!(num.H, num.H, fwd.usave[i,:,:]', levels = 0:0, color=:black, linewidth = 3);
    end
f = contour!(num.H, num.H, fwd.usave[end,:,:]', levels = 0:0, color=:black, linewidth = 3);
end

f = current_figure()
