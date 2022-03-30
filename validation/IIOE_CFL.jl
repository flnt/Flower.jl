using Revise
using Flower

N_array = 2 .^[4, 5, 6]

CFL_array = [1.0, 4.0, 16.0]

plt2 = ("1", "2", "∞")

function conv_IIOE_CFL(x, y;
    case = "Sphere",
    speed = 1,
    L0 = 2.,
    ii = 1,
    R = 0.5,
    verbose = true)

    N = length(x)
    C = length(y)
    r = (x[end]+1)/(x[end-1]+1)
    ERR = zeros(3,N,C)
    u_data = Vector{Matrix{Float64}}(undef, 0)
    ua_data = Vector{Matrix{Float64}}(undef, 0)
    e_data = Vector{Matrix{Float64}}(undef, 0)
    grid = Vector{Vector{Float64}}(undef, 0)
    R += 1/x[1]

    for j in 1:C
        for i in 1:N
            num = Numerical(case = case,
                L0 = L0,
                n = x[i] + 1,
                CFL = y[j],
                max_iterations = Int(ii*4^(i)/y[j]),
                R = R,
                nb_reinit = 0);

            idx = set_indices(num.n);
            tmp, fwd = init_fields(num, idx);

            MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
            advection = true,
            analytical = true,
            speed = speed,
            verbose = verbose
            )

            @show (R + speed*num.max_iterations*num.τ, num.CFL)

            ua = sqrt.(num.X .^ 2 + num.Y .^ 2) - (R + speed*num.max_iterations*num.τ) * ones(num.n, num.n);
            e = similar(ua)
            for II in idx.inside
                e[II] = abs(fwd.u[II] - ua[II])
            end
            ERR[:, i, j] .= normf(e, idx.inside, num.Δ)
            push!(u_data, fwd.u)
            push!(ua_data, sqrt.(num.X .^ 2 + num.Y .^ 2) - R * ones(num.n, num.n))
            push!(e_data, e)
            push!(grid, num.H)
        end
    end

    return u_data, ua_data, grid, ERR, e_data
end

u_data, ua_data, grid, ERR, e_data = conv_IIOE_CFL(N_array,CFL_array, ii = 8, verbose = true, speed = -1., R = 0.8)


fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)
f = Figure()
for i in (3,6,9)
    c = CFL_array[i÷3]
    ax = Axis(f[1, i], title = @sprintf "N = 64, CFL = %d" c)
    colsize!(f.layout, i, Aspect(1, 1.0))
    hidedecorations!(ax)
    ct = heatmap!(f[1,i], grid[i][2:end-1], grid[i][2:end-1], (e_data[i][2:end-1,2:end-1])')
    contour!(f[1,i], grid[i][2:end-1], grid[i][2:end-1], (ua_data[i][2:end-1,2:end-1])', levels = 0:0, linewidth = 3, color =:red)
    contour!(f[1,i], grid[i][2:end-1], grid[i][2:end-1], (u_data[i][2:end-1,2:end-1])', levels = 0:0, linewidth = 3, color =:white)
end
Colorbar(f[1, 10], flipaxis = true, label = "Normalized error", ticks = 0:1)
resize_to_layout!(f)

f2 = Figure()
fontsize_theme = Theme(fontsize = 15)
set_theme!(fontsize_theme)

colsize!(f2.layout, 1, Aspect(1, 1.0))
ax = Axis(f2[1,1], xscale = log2, yscale = log10,
xticks = N_array, xlabel = "Points per dimension", ylabel = "Error")

resize_to_layout!(f2)


for j in 1:length(CFL_array)
    for i in (2) #type of norm (1, 2, inf)
        c = CFL_array[j]
        n = N_array[j]
        a, b = fit_order(N_array, ERR[i, :, j])
        f2 = scatter!(N_array, ERR[i, :, j], markersize = 10, label = @sprintf "CFL = %d, order %.2f" c b)
        f2 = lines!(N_array, a*(1 ./ N_array).^b)
        f2 = lines!([n, n, n],ERR[i,j,:] , linestyle =:dash, color=:black)
    end
end


axislegend(position = :lb)

f2 = current_figure()


f = Figure(resolution = (1200, 1200))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
for i in (3,6,9)
    j = i÷3
    c = CFL_array[i÷3]
    ax = Axis(f[1, j], title = @sprintf "N = 64, CFL = %d" c)
    colsize!(f.layout, j, Aspect(1, 1.0))
    hidedecorations!(ax)
    ct = heatmap!(f[1,j], grid[i][2:end-1], grid[i][2:end-1], (e_data[i][2:end-1,2:end-1])')
    contour!(f[1,j], grid[i][2:end-1], grid[i][2:end-1], (ua_data[i][2:end-1,2:end-1])', levels = 0:0, linewidth = 4, color =:red)
    contour!(f[1,j], grid[i][2:end-1], grid[i][2:end-1], (u_data[i][2:end-1,2:end-1])', levels = 0:0, linewidth = 4, color =:white)
end
Colorbar(f[1, 4], flipaxis = true, label = "Normalized error", ticks = 0:1)


ax = Axis(f[2:3,1:4], xscale = log2, yscale = log10,
xticks = N_array, xlabel = "Points per dimension", ylabel = "Error")

scatter!(f[2:3,1:4], N_array, ERR[2, :, 1], label = "CFL = 1", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, ERR[2, :, 1])
lines!(f[2:3,1:4], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

scatter!(f[2:3,1:4], N_array, ERR[2, :, 2], label = "CFL = 4", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, ERR[2, :, 2])
lines!(f[2:3,1:4], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)

scatter!(f[2:3,1:4], N_array, ERR[2, :, 3], label = "CFL = 16", markersize = 25, color =:black, marker=:diamond)
a, b = fit_order(N_array, ERR[2, :, 3])
lines!(f[2:3,1:4], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dot, color =:black, label = @sprintf "order %.2f" b)

lines!(f[2:3,1:4], [16, 16, 16], ERR[2, 1, :], linewidth = 3, color =:blue)
lines!(f[2:3,1:4], [32, 32, 32], ERR[2, 2, :], linewidth = 3, color =:blue)
lines!(f[2:3,1:4], [64, 64, 64], ERR[2, 3, :], linewidth = 3, color =:blue)

axislegend(position = :lb)
resize_to_layout!(f)
rm("fit.log")
f = current_figure()

#Makie.save("./figures/paper_figures/L2_CFL1_4_16_n16_32_64_retracting.png", f)
