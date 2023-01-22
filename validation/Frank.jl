using Revise
using Flower

function init_franck!(grid, temp, R, T_inf, h)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind

    @inbounds for II in all_indices
        s = sqrt(x[II]^2 + y[II]^2)
        if s >= (R - h)
            temp[II] = T_inf*(1-(expint(0.25*s^2))/(expint(0.25*R^2)))
        end
    end
end

function init_franck!(grid, temp, R, T_inf, h, t)
    @unpack x, y, nx, ny, ind = grid
    @unpack all_indices = ind

    @inbounds for II in all_indices
        s = sqrt(x[II]^2 + y[II]^2)/√t
        if s >= (R - h)
            temp[II] = T_inf*(1-(expint(0.25*s^2))/(expint(0.25*R^2)))
        end
    end
end

N_array = 2 .^[5, 6, 7]
plt = ("partial", "full", "all")
plt2 = ("1", "2", "∞")

function conv_Frank(arr;
    case = "Sphere",
    CFL = 0.5,
    L0 = 16.,
    R = 1.56,
    TEND = 1.0,
    T_inf = -0.5,
    solid = false,
    liquid = true,
    verbose = true)

    N = length(arr)
    ERR = zeros(3,3,N,2)
    T_data = Vector{Matrix{Float64}}(undef, 0)
    u_data1 = Vector{Matrix{Float64}}(undef, 0)
    u_data2 = Vector{Matrix{Float64}}(undef, 0)
    u_data3 = Vector{Matrix{Float64}}(undef, 0)
    grid = Vector{Vector{Float64}}(undef, 0)
    radius_data = Vector{Vector{Float64}}(undef, 0)
    time_data = Vector{Vector{Float64}}(undef, 0)
    analytical_radius = [R*sqrt(1 + i) for i = 0:0.001:TEND];
    analytical_temperature_data = Vector{Matrix{Float64}}(undef, 0)

    for i in 1:N
        x = LinRange(-L0/2, L0/2, arr[i]+1)
        y = LinRange(-L0/2, L0/2, arr[i]+1)

        num = Numerical(case = case,
            T_inf = T_inf,
            L0 = L0,
            x = x,
            y = y,
            CFL = CFL,
            TEND = TEND,
            R = R,
            ϵ = 0.00
            );
        gp, gu, gv = init_meshes(num)
        opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
        @time MIXED, SOLID, LIQUID, radius = run_forward(num, gp, gu, gv,
            opS, opL, phS, phL, fwd,
            stefan = true,
            advection = true,
            heat = true,
            heat_solid_phase = false,
            heat_liquid_phase = true,
            verbose = verbose,
            show_every = 1,
            save_radius = true,
            Vmean = true
        )

        analytical_temperature = similar(phL.T)
        analytical_temperature .= 0.
        init_franck!(gp, analytical_temperature, R, T_inf, 0, 2)
        push!(u_data1, fwd.usave[1,:,:])
        push!(u_data2, fwd.usave[end÷2,:,:])
        push!(u_data3, fwd.usave[end,:,:])

        push!(T_data, phL.T)
        push!(analytical_temperature_data, analytical_temperature)


        push!(grid, gp.x[1,:])
        push!(radius_data, radius)
        time = [1 + i*num.τ for i in 0:num.max_iterations]
        push!(time_data, time)

        e = analytical_temperature - phL.T

        if solid
            ERR[1, :, i, 1] .= normf(e, MIXED, gp.geoS.cap[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(e, SOLID, gp.geoS.cap[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(e, vcat(SOLID, MIXED), gp.geoS.cap[:,:,5], num.Δ)
        elseif liquid
            ERR[1, :, i, 1] .= normf(e, MIXED, gp.geoL.cap[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(e, LIQUID, gp.geoL.cap[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(e, vcat(LIQUID, MIXED), gp.geoL.cap[:,:,5], num.Δ)
        end

        for i = 1:3
            for j = 1:3
                ERR[j, i, :, 2] .= Richardson_extrapolation(ERR[j, i, :, 1], 2)
            end
        end
    end
    return ERR, grid, analytical_radius, analytical_temperature_data, T_data, radius_data, time_data, u_data1, u_data2, u_data3
end

ERR, grid, analytical_radius, analytical_temperature_data, T_data, radius_data, time_data, u_data1, u_data2, u_data3 = conv_Frank(N_array, verbose = true, TEND= 1.0, CFL = 0.1);

N = length(N_array)

err_radius = zeros(N)

for i in 1:N
    err_radius[i] = abs(analytical_radius[end] - radius_data[i][end])
end


err_map = similar(analytical_temperature_data[3])

for II in CartesianIndices(err_map)
    a = abs(analytical_temperature_data[3][II] - T_data[3][II])
    err_map[II] = ifelse(a == 0., NaN, a)
end


my_colors = ("green", "blue", "red")













f = Figure(resolution = (1200, 1200))
# f = Figure()
fontsize_theme = Theme(fontsize = 40)
set_theme!(fontsize_theme)

# ax = Axis(f[1, 1], xlabel = "Dimensionless time", ylabel = "Dimensionless radius")
# colsize!(f.layout, 1, Aspect(1, 2))
# resize_to_layout!(f)
# xlims!(1, 2.02)
# ylims!(1.56, 2.22)

# ax.xticks = [1.0 + 0.1*i for i = 0:10]

# for i in 1:N
#     scatter!(f[1,1],time_data[i][end:-2*i:2], radius_data[i][end:-2*i:2], color=(my_colors[i], 1.0), markersize = 10, label = "N = $(N_array[i])")
# end
# lines!(f[1,1],1:1/(length(analytical_radius)-1):2, analytical_radius, label = "Analytical solution", linewidth = 5,color=:black)



struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

ax = Axis(f[1:2,1], xscale = log2, yscale = log10,
xticks = N_array, yticks = LogTicks(IntegerTicks()), xlabel = "Points per dimension", ylabel = "Error")
colsize!(f.layout, 1, Aspect(1, 2))
XX = ERR[1, 2, :, 1]

scatter!(f[1:2,1], N_array, XX, label = "Mixed cells, L-2 norm", markersize = 35, color =:black, marker=:rect)
a, b = fit_order(N_array, XX)
lines!(f[1:2,1], N_array, a*(1 ./ N_array).^b, linewidth = 5, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

YY = ERR[2, 2, :, 1]

scatter!(f[1:2,1], N_array, YY, label = "Full cells, L-2 norm", markersize = 35, color =:black, marker=:utriangle)
a, b = fit_order(N_array, YY)
lines!(f[1:2,1], N_array, a*(1 ./ N_array).^b, linewidth = 5, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)

ZZ = ERR[3, 2, :, 1]

scatter!(f[1:2,1], N_array, ZZ, label = "All cells, L-2 norm", markersize = 35, color =:black, marker=:circle)
a, b = fit_order(N_array, ZZ)
lines!(f[1:2,1], N_array, a*(1 ./ N_array).^b, linewidth = 5, linestyle =:dot, color =:black, label = @sprintf "order %.2f" b)

axislegend(position = :lb)


ax3 = Axis(f[1:2,2], xlabel = "x")
colsize!(f.layout, 2, Aspect(1, 2))
xlims!(-2.5, 2.5)
ylims!(-2.5, 2.5)

ax3.xticks = [1.560, 2.206]
hideydecorations!(ax3)
hm = heatmap!(f[1:2,2], grid[N], grid[N], err_map)#, colormap=Reverse(:greys))
resize_to_layout!(f)
Colorbar(f[1:2, 3], vertical = true, flipaxis = true,  ticks = 0:1, label = "Normalized error")#, colormap=Reverse(:greys))

for i in 1:N
    lines!(f[1:2,2],[-10000, -99999], [-10000, -99999], color=my_colors[i], linewidth = 5, label = "N = $(N_array[i])")
    contour!(f[1:2,2],grid[i], grid[i], u_data1[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1:2,2],grid[i], grid[i], u_data2[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1:2,2],grid[i], grid[i], u_data3[i], levels = 0:0, color=my_colors[i], linewidth = 5)
end
vlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)
hlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)

axislegend(position = :lb)

resize_to_layout!(f)
f = current_figure()

Makie.save("./figures/paper_figures/Frank_convergence_new.png", f)





















f = Figure(resolution = (1200, 1200))
# f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)

ax = Axis(f[1, 1], xlabel = "Dimensionless time", ylabel = "Dimensionless radius")
colsize!(f.layout, 1, Aspect(1, 2))
resize_to_layout!(f)
xlims!(1, 2.02)
ylims!(1.56, 2.22)

ax.xticks = [1.0 + 0.1*i for i = 0:10]

for i in 1:N
    scatter!(f[1,1],time_data[i][end:-2*i:2], radius_data[i][end:-2*i:2], color=(my_colors[i], 1.0), markersize = 10, label = "N = $(N_array[i])")
end
lines!(f[1,1],1:1/(length(analytical_radius)-1):2, analytical_radius, label = "Analytical solution", linewidth = 5,color=:black)

axislegend(position = :lt)

# Makie.save("./figures/paper_figures/Frank_radius.png", f)









# f = Figure(resolution = (800, 800))
f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)

ax = Axis(f[1,1], xscale = log2, yscale = log10,
xticks = N_array, yticks = LogTicks(IntegerTicks()), xlabel = "Points per dimension", ylabel = "Error")

XX = ERR[1, 2, :, 1]

scatter!(f[1,1], N_array, XX, label = "Mixed cells, L-2 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, XX)
lines!(f[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

YY = ERR[2, 2, :, 1]

scatter!(f[1,1], N_array, YY, label = "Full cells, L-2 norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, YY)
lines!(f[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)

ZZ = ERR[3, 2, :, 1]

scatter!(f[1,1], N_array, ZZ, label = "All cells, L-2 norm", markersize = 25, color =:black, marker=:circle)
a, b = fit_order(N_array, ZZ)
lines!(f[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dot, color =:black, label = @sprintf "order %.2f" b)


axislegend(position = :lb)
resize_to_layout!(f)

# Makie.save("./figures/paper_figures/Frank_error_temp.png", f)





f = Figure(resolution = (1200, 1200))
ax3 = Axis(f[1,1], xlabel = "x")

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

xlims!(-3, 3)
ylims!(-3, 3)

ax3.yticks = [1.560, 2.206]

ax3.xticks = [1.560, 2.206]
hideydecorations!(ax3)
hm = heatmap!(f[1,1], grid[N], grid[N], err_map)#, colormap=Reverse(:greys))
resize_to_layout!(f)
Colorbar(f[1, 2], vertical = true, flipaxis = true,  ticks = 0:1, label = "Normalized error")#, colormap=Reverse(:greys))

for i in 1:N
    lines!(f[1,1],[-10000, -99999], [-10000, -99999], color=my_colors[i], linewidth = 5, label = "N = $(N_array[i])")
    contour!(f[1,1],grid[i], grid[i], u_data1[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1,1],grid[i], grid[i], u_data2[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1,1],grid[i], grid[i], u_data3[i], levels = 0:0, color=my_colors[i], linewidth = 5)
end
vlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)
hlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)

axislegend(position = :lb)

resize_to_layout!(f)
f = current_figure()

# Makie.save("./figures/paper_figures/Frank_errormap_temp.png", f)