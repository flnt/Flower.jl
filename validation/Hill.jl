using Revise
using Flower

N_array = [16, 32, 64]

function conv_Hill(arr;
    case = "Planar",
    CFL = 0.5,
    L0 = 2.,
    TEND = 1.0,
    T_inf = -1.0,
    solid = false,
    liquid = true,
    verbose = true)

    N = length(arr)
    ERR = zeros(3,3,N,2)
    grid = Vector{Vector{Float64}}(undef, 0)
    radius_data = Vector{Vector{Float64}}(undef, 0)
    time_data = Vector{Vector{Float64}}(undef, 0)
    u_data1 = Vector{Matrix{Float64}}(undef, 0)
    u_data2 = Vector{Matrix{Float64}}(undef, 0)
    u_data3 = Vector{Matrix{Float64}}(undef, 0)
    analytical_temperature_data = Vector{Matrix{Float64}}(undef, 0)
    T_data = Vector{Matrix{Float64}}(undef, 0)

    for i in 1:N
        x = LinRange(0., 1.0, arr[i]+1)
        y = LinRange(0., 1.0, arr[i]+1)

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
            verbose = true,
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

        push!(radius_data, radius)

        time = [i*num.τ for i in 0:num.max_iterations]
        push!(time_data, time)
        push!(grid, gp.x[1,:])

        push!(u_data1, fwd.usave[1,:,:])
        push!(u_data2, fwd.usave[end÷2,:,:])
        push!(u_data3, fwd.usave[end,:,:])

        # e = phL.T - TL2

        push!(T_data, phL.T)
        push!(analytical_temperature_data, TLanalytical)

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
    return ERR, radius_data, time_data, grid, u_data1, u_data2, u_data3, analytical_temperature_data, T_data
end

TEND = 0.01
ERR, radius_data, time_data, grid, u_data1, u_data2, u_data3, analytical_temperature_data, T_data = conv_Hill(N_array, TEND = TEND, CFL = 0.5, verbose = false);

N = length(N_array)

my_colors = ("green", "blue", "red")

analytical_pos = LinRange(0, TEND, 100)

err_map = similar(analytical_temperature_data[3])

for II in CartesianIndices(err_map)
    a = abs(analytical_temperature_data[3][II] - T_data[3][II])
    err_map[II] = ifelse(a == 0., NaN, a)
end

# f = Figure(resolution = (1200, 1200))
fradius = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)

ax = Axis(fradius[1, 1], xlabel = "Dimensionless time", ylabel = "Dimensionless position")
colsize!(fradius.layout, 1, Aspect(1, 1.375))

# ax.xticks = [0, TEND]
# ax.yticks = [0, TEND]


for i in 1:N
    scatter!(fradius[1,1],time_data[i][end:-3*i:2], radius_data[i][end:-3*i:2], color=(my_colors[i], 0.7), markersize = 10, label = "N = $(N_array[i])")
end

lines!(fradius[1,1],analytical_pos, analytical_pos, label = "Analytical solution", linewidth = 5,color=:black)

axislegend(position = :lt)


# Makie.save("./figures/paper_figures/Hill_radius.png", fradius)











fconv = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)

ax = Axis(fconv[1,1], xscale = log2, yscale = log10,
xticks = N_array, yticks = LogTicks(IntegerTicks()), xlabel = "Points per dimension", ylabel = "Error")

# XX = ERR[1, 2, :, 2]

# scatter!(fconv[1,1], N_array, XX, label = "Mixed cells, L-2 norm", markersize = 25, color =:black, marker=:rect)
# a, b = fit_order(N_array, XX)
# lines!(fconv[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

ZZ = ERR[3, 2, :, 2]

scatter!(fconv[1,1], N_array, ZZ, label = "All cells, L-2 norm", markersize = 25, color =:black, marker=:circle)
a, b = fit_order(N_array, ZZ)
lines!(fconv[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dot, color =:black, label = @sprintf "order %.2f" b)

YY = ERR[3, 3, :, 1]

scatter!(fconv[1,1], N_array, YY, label = "All cells, L-inf norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, YY)
lines!(fconv[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)


axislegend(position = :lb)
resize_to_layout!(fconv)


Makie.save("./figures/paper_figures/Hill_conv_temp.png", fconv)



ferrmap = Figure(resolution = (1200, 1200))
ax3 = Axis(ferrmap[1,1], xlabel = "x")

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

# xlims!(-3, 3)
# ylims!(-3, 3)

# ax3.yticks = [1.560, 2.206]

# ax3.xticks = [1.560, 2.206]
hideydecorations!(ax3)
hm = heatmap!(ferrmap[1,1], grid[N], grid[N], err_map)#, colormap=Reverse(:greys))
resize_to_layout!(ferrmap)
Colorbar(ferrmap[1, 2], vertical = true, flipaxis = true,  ticks = 0:1, label = "Normalized error")#, colormap=Reverse(:greys))

for i in 1:N
    # lines!(ferrmap[1,1],[-10000, -99999], [-10000, -99999], color=my_colors[i], linewidth = 5, label = "N = $(N_array[i])")
    contour!(ferrmap[1,1],grid[i], grid[i], u_data1[i], levels = 0:0, color=my_colors[i], linewidth = 1)
    contour!(ferrmap[1,1],grid[i], grid[i], u_data2[i], levels = 0:0, color=my_colors[i], linewidth = 1)
    contour!(ferrmap[1,1],grid[i], grid[i], u_data3[i], levels = 0:0, color=my_colors[i], linewidth = 1)
end
# vlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)
# hlines!(ax3, [1.560, 2.206], color = :black, linestyle =:dash, linewidth = 3)

# axislegend(position = :lb)

resize_to_layout!(ferrmap)


# Makie.save("./figures/paper_figures/Frank_errormap_temp.png", ferrmap)