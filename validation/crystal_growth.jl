using Revise
using Flower

eps_array = [0.0004, 0.0006, 0.0008, 0.001]

function surface_tension_effect(x;
    CFL = 0.5,
    n = 200,
    L0 = 4.,
    TEND = 0.5,
    T_inf = -0.5
    )

    N = length(x)

    length_data = Vector{Vector{Float64}}(undef, 0)
    grid = Vector{Vector{Float64}}(undef, 0)
    time_data = Vector{Vector{Float64}}(undef, 0)
    κ_max_data = Vector{Vector{Float64}}(undef, 0)
    κ_min_data = Vector{Vector{Float64}}(undef, 0)
    T_data = Vector{Array{Float64,2}}(undef, 0)
    u_data = Vector{Array{Float64,3}}(undef, 0)

    for i in 1:N
        num = Numerical(T_inf = T_inf,
            ϵ_κ = x[i],
            case = "Crystal",
            L0 = 4.,
            n = n,
            CFL = CFL,
            TEND = TEND
            )

        idx, idxu, idxv = set_indices(num.n)
        tmp, fwd = init_fields(num, idx, idxu, idxv)
        fwd.TL[:,:] .= num.T_inf
        fwd.u = @. (sqrt((num.X)^2 + num.Y^2)*(0.1 + 0.02*cos(4*atan((num.X)/num.Y))) - 0.01);

        MIXED = run_forward(num, idx, idxu, idxv, tmp, fwd,
            stefan = true,
            advection = true,
            heat = true,
            liquid_phase = true,
            solid_phase = true,
            show_every = 100,
            verbose = true,
            save_length = true
            )

        fwd.lengthsave ./= fwd.lengthsave[1]

        κ_min = zeros(num.max_iterations+1)
        κ_max = zeros(num.max_iterations+1)

        for j in 1:num.max_iterations+1
            κ_min[j] = findmin(fwd.κsave[j,:,:])[1]
            κ_max[j] = findmax(fwd.κsave[j,:,:])[1]
        end

        if i == 1
            push!(grid, num.H)
            time = [i*num.τ for i in 0:num.max_iterations]
            push!(time_data, time)
        end

        push!(length_data, fwd.lengthsave)
        push!(κ_max_data, κ_max)
        push!(κ_min_data, κ_min)
        push!(T_data, fwd.TL)
        push!(u_data, fwd.usave)
    end
    return length_data, κ_max_data, κ_min_data, u_data, T_data, time_data, grid
end

length_data, κ_max_data, κ_min_data, u_data, T_data, time_data, grid = surface_tension_effect(eps_array, TEND = 0.6)

step = length(time_data[1])÷20

step2 = step÷5

my_colors = ("blue", "orange", "green", "purple")

f = Figure()
fontsize_theme = Theme(fontsize = 15)
set_theme!(fontsize_theme)

ax = Axis(f[1, 1], xlabel = "Dimensionless time", ylabel = "Normalized interface length",yscale = log10)
colsize!(f.layout, 1, Aspect(1, 1.375))

ax.xticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

for i in 1:4
    scatter!(f[1,1],time_data[1][1:step2:end], length_data[i][1:step2:end], label = "ϵ_κ = $(eps_array[i])", linewidth = 3, color = my_colors[i])
end

axislegend(position = :lt)

resize_to_layout!(f)
f = current_figure()


#Makie.save("./figures/paper_figures/crystal_interface_vs_time_logscale.png", f)

f2 = Figure(resolution = (4000, 4000))
fontsize_theme = Theme(fontsize = 80)
set_theme!(fontsize_theme)

for i in 1:4
    j = i
    ix = 1
    if i > 2
        ix = 2
        j = i-2
    end
    ax = Axis(f2[ix, j], title = @sprintf "ϵ_κ = %g" eps_array[i])#$(N_array[i])
    colsize!(f2.layout, j, Aspect(1, 1.0))
    hidedecorations!(ax)
    heatmap!(f2[ix, j], grid[1], grid[1], T_data[i]', colormap = :ice)
    for ii in 1:step:length(time_data[1])
        contour!(f2[ix, j], grid[1], grid[1], u_data[i][ii,:,:]', levels = 0:0, color=:black, linewidth = 3);
    end
end
f2 = current_figure()

#Makie.save("./figures/paper_figures/position_eps0_0002_to_eps0_0008_2.png", f2)
