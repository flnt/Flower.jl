using Revise
using Flower

# using Revise
# using Flower

# fontsize_theme = Theme(fontsize = 30)
# set_theme!(fontsize_theme)

# ratio = 1
# L0 = 16.
# nx = 128
# ny = ratio * nx

# x = LinRange(-L0/2, L0/2, nx+1)
# y = LinRange(-ratio*L0/2, ratio*L0/2, ny+1)

# num = Numerical(case = "Sphere",
#     CFL = 0.1,
#     TEND = 1.0,
#     x = x,
#     y = y,
#     R = 1.56,
#     save_every = 1,
#     ϵ = 0.05,
#     θd = 0.0,
#     T_inf = -0.5
# )

# gp, gu, gv = init_meshes(num)
# opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

# @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
#     opS, opL, phS, phL, fwd,
#     stefan = true,
#     advection = true,
#     heat = true,
#     heat_solid_phase = false,
#     heat_liquid_phase = true,
#     verbose = true,
#     show_every = 1,
# )

# tcks = -ratio*L0/2:2:ratio*L0/2
# lim = L0 / 2

# fsolid = Figure(resolution = (1600, 1000))
# colsize!(fsolid.layout, 1, Aspect(1, 1.0))
# ax = Axis(fsolid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gp.x[1,:], gp.y[:,1], phS.T')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# # limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fsolid)

# # fsolid = current_figure()

# fliquid = Figure(resolution = (1600, 1000))
# colsize!(fliquid.layout, 1, Aspect(1, 1.0))
# ax = Axis(fliquid[1,1], aspect = 1/ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# heatmap!(gp.x[1,:], gp.y[:,1], phL.T')
# contour!(gp.x[1,:], gp.y[:,1], gp.u', levels = 0:0, color=:red, linewidrth = 3);
# contour!(gp.x[1,:], gp.y[:,1], fwd.usave[1,:,:]', levels = 0:0, color=:blue, linewidrth = 3);
# # limits!(ax, -lim, lim, -lim, lim)
# resize_to_layout!(fliquid)

# # fliquid = current_figure()

N_array = 2 .^[5, 6, 7]
plt = ("partial", "full", "all")
plt2 = ("1", "2", "∞")

function conv_Frank(x;
    case = "Sphere",
    L0 = 16.,
    R = 1.56,
    TEND = 1.0,
    T_inf = -0.5,
    solid = false,
    liquid = true,
    verbose = true)

    N = length(x)
    r = (x[end])/(x[end-1])
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
        num = Numerical(case = case,
            T_inf = T_inf,
            L0 = L0,
            n = x[i],
            CFL = 0.1,
            TEND = TEND,
            R = R
            );
        idx, idxu, idxv = set_indices(num.n);
        tmp, fwd = init_fields(num, idx, idxu, idxv);
        MIXED, SOLID, LIQUID, radius = run_forward(num, idx, idxu, idxv, tmp, fwd,
        heat = true,
        solid_phase = solid,
        liquid_phase = liquid,
        stefan = true,
        advection = true,
        verbose = verbose,
        save_radius = true
        )

        analytical_temperature = zeros(num.n,num.n)
        init_franck!(analytical_temperature, R*sqrt(1+TEND), T_inf, num.H, num.n, 0)

        push!(u_data1, fwd.usave[1,:,:])
        push!(u_data2, fwd.usave[end÷2,:,:])
        push!(u_data3, fwd.usave[end,:,:])

        push!(T_data, fwd.TL)
        push!(analytical_temperature_data, analytical_temperature)


        push!(grid, num.H)
        push!(radius_data, radius)
        time = [1 + i*num.τ for i in 0:num.max_iterations]
        push!(time_data, time)

        e = analytical_temperature - fwd.TL

        if solid
            ERR[1, :, i, 1] .= normf(e, MIXED, tmp.SOL[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(e, SOLID, tmp.SOL[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(e, vcat(SOLID, MIXED), tmp.SOL[:,:,5], num.Δ)
        elseif liquid
            ERR[1, :, i, 1] .= normf(e, MIXED, tmp.LIQ[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(e, LIQUID, tmp.LIQ[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(e, vcat(LIQUID, MIXED), tmp.LIQ[:,:,5], num.Δ)
        end

        for i = 1:3
            for j = 1:3
                ERR[j, i, :, 2] .= Richardson_extrapolation(ERR[j, i, :, 1], 2)
            end
        end
    end
    return ERR, grid, analytical_radius, analytical_temperature_data, T_data, radius_data, time_data, u_data1, u_data2, u_data3
end

ERR, grid, analytical_radius, analytical_temperature_data, T_data, radius_data, time_data, u_data1, u_data2, u_data3 = conv_Frank(N_array, verbose = false, TEND= 1.0)

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
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax = Axis(f[1, 1], xlabel = "Dimensionless time", ylabel = "Dimensionless radius")
colsize!(f.layout, 1, Aspect(1, 2))
resize_to_layout!(f)
xlims!(1, 2.02)
ylims!(1.56, 2.22)

ax.xticks = [1.0 + 0.1*i for i = 0:10]

for i in 1:N
    scatter!(f[1,1],time_data[i][end:-3*i:2], radius_data[i][end:-3*i:2], color=(my_colors[i], 1.0), markersize = 10, label = "N = $(N_array[i])")
end
lines!(f[1,1],1:1/(length(analytical_radius)-1):2, analytical_radius, label = "Analytical solution", linewidth = 5,color=:black)

axislegend(position = :lt)

ax2 = Axis(f[2, 1], xscale = log2, yscale = log2,
xticks = N_array, xlabel = "Points per dimension", ylabel = "Error")


scatter!(f[2,1], N_array, err_radius, label = "Radius error", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, err_radius)
lines!(f[2,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)


axislegend(position = :lb)


ax3 = Axis(f[1:2,2], xlabel = "x")

xlims!(-2.5, 2.5)
ylims!(-2.5, 2.5)

ax3.xticks = [1.560, 2.206]
hideydecorations!(ax3)
hm = heatmap!(f[1:2,2], grid[N], grid[N], err_map, colormap=:greys)
resize_to_layout!(f)
Colorbar(f[1:2, 3], vertical = true, flipaxis = true,  ticks = 0:1, label = "Normalized error", colormap=:greys)

for i in 1:N
    lines!(f[1:2,2],[-10000, -99999], [-10000, -99999], color=my_colors[i], linewidth = 5, label = "N = $(N_array[i])")
    contour!(f[1:2,2],grid[i], grid[i], u_data1[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1:2,2],grid[i], grid[i], u_data2[i], levels = 0:0, color=my_colors[i], linewidth = 5)
    contour!(f[1:2,2],grid[i], grid[i], u_data3[i], levels = 0:0, color=my_colors[i], linewidth = 5)
end
vlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)
hlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)

axislegend(position = :lb)
rm("fit.log")
resize_to_layout!(f)
f = current_figure()

#Makie.save("./figures/paper_figures/Frank_convergence.png", f)
