using Revise
using Flower

N_array = [50, 100, 150, 200]

function grid_effect(x;
    CFL = 0.5,
    ϵ_κ = 0.0004,
    ϵ_V = 0.002,
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
            ϵ_κ = ϵ_κ,
            ϵ_V = ϵ_V,
            case = "Crystal",
            L0 = 4.,
            n = x[i],
            CFL = CFL,
            TEND = TEND
            )

        idx = set_indices(num.n)
        tmp, fwd = init_fields(num, idx)
        fwd.TL[:,:] .= num.T_inf
        fwd.u = @. (sqrt((num.X)^2 + num.Y^2)*(0.1 + 0.02*cos(4*atan((num.X)/num.Y))) - 0.01);

        MIXED = run_forward(num, idx, tmp, fwd,
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


        push!(grid, num.H)
        time = [i*num.τ for i in 0:num.max_iterations]
        push!(time_data, time)


        push!(length_data, fwd.lengthsave)
        push!(κ_max_data, κ_max)
        push!(κ_min_data, κ_min)
        push!(T_data, fwd.TL)
        push!(u_data, fwd.usave)
    end
    return length_data, κ_max_data, κ_min_data, u_data, T_data, time_data, grid
end

length_data, κ_max_data, κ_min_data, u_data, T_data, time_data, grid = grid_effect(N_array,
CFL = 0.5,
ϵ_κ = 0.0004,
ϵ_V = 0.000,
L0 = 4.,
TEND = 0.6,
T_inf = -0.5)

steps = Vector{Int64}(undef, 0)
for i in 1:length(N_array)
    push!(steps, length(time_data[i])÷20)
end

f = Figure(resolution = (4000, 4000))
fontsize_theme = Theme(fontsize = 80)
set_theme!(fontsize_theme)

for i in 1:4
    j = i
    ix = 1
    if i > 2
        ix = 2
        j = i-2
    end
    ax = Axis(f[ix, j], title = @sprintf "N = %d" N_array[i])
    colsize!(f.layout, j, Aspect(1, 1.0))
    hidedecorations!(ax)
    heatmap!(f[ix, j], grid[i], grid[i], T_data[i]', colormap = :ice)
    for ii in 1:steps[i]:length(time_data[i])
        contour!(f[ix, j], grid[i], grid[i], u_data[i][ii,:,:]', levels = 0:0, color=:black, linewidth = 3);
    end
end
f = current_figure()

#Makie.save("./figures/paper_figures/grid_effect_eps0_0004_2.png", f)
