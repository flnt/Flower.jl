using Revise
using Flower

N_array = 2 .^[5, 6, 7]
plt = ("partial", "full", "all")
plt2 = ("1", "2", "∞")

function conv_cutcell_CN(x;
    case = "Sphere",
    L0 = 2.,
    Θd = 1.0,
    ii = 4,
    R = 0.75,
    solid = true,
    liquid = false,
    verbose = true)

    N = length(x)
    r = (x[end]+1)/(x[end-1]+1)
    ERR = zeros(3,3,N,2)
    T_data = Vector{Matrix{Float64}}(undef, 0)
    u_data = Vector{Matrix{Float64}}(undef, 0)
    grid = Vector{Vector{Float64}}(undef, 0)

    for i in 1:N
        num = Numerical(case = case,
            L0 = L0,
            θd = Θd,
            n = x[i] + 1,
            max_iterations = ii*4^(i),
            R = R + 1/x[1]);
        idx, idxu, idxv = set_indices(num.n);
        tmp, fwd = init_fields(num, idx, idxu, idxv);
        MIXED, SOLID, LIQUID = run_forward(num, idx, idxu, idxv, tmp, fwd,
        heat = true,
        solid_phase = solid,
        liquid_phase = liquid,
        verbose = verbose
        )
        @show (num.max_iterations*num.τ)
        @show (num.τ)
        if solid
            ERR[1, :, i, 1] .= normf(fwd.TS, MIXED, tmp.SOL[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(fwd.TS, SOLID, tmp.SOL[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(fwd.TS, vcat(SOLID, MIXED), tmp.SOL[:,:,5], num.Δ)
            push!(T_data, fwd.TS)
        elseif liquid
            ERR[1, :, i, 1] .= normf(fwd.TL, MIXED, tmp.LIQ[:,:,5], num.Δ)
            ERR[2, :, i, 1] .= normf(fwd.TL, LIQUID, tmp.LIQ[:,:,5], num.Δ)
            ERR[3, :, i, 1] .= normf(fwd.TL, vcat(LIQUID, MIXED), tmp.LIQ[:,:,5], num.Δ)
            push!(T_data, fwd.TL)
        end
        push!(u_data, fwd.u)
        push!(grid, num.H)
    end
    try
        @assert N >= 3
    catch
        @warn "Need 3 points or more for Richardson extrapolation"
        return T_data, u_data, grid, ERR
    end
    for i = 1:3
        for j = 1:3
            ERR[j, i, :, 2] .= Richardson_extrapolation(ERR[j, i, :, 1], r)
        end
    end
    return T_data, u_data, grid, ERR
end

T_data, u_data, grid, ERR = conv_cutcell_CN(N_array,
case = "Sphere",
R = 0.75,
ii = 4,
verbose = false)


f = Figure(resolution = (1200, 1200))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

for i in 1:length(N_array)
    ax = Axis(f[1, i], title = "N = $(N_array[i])")
    colsize!(f.layout, i, Aspect(1, 1.0))
    hidedecorations!(ax)
    hm = heatmap!(f[1,i], grid[i], grid[i], T_data[i]', colormap = Reverse(:ice))
    contour!(f[1,i], grid[i],grid[i],u_data[i]', levels = 0:0, linewidth = 4, color =:red)
end

Colorbar(f[1, 4], vertical = true, flipaxis = true, label = "Temperature", ticks = 0:1, colormap = Reverse(:ice))
ax = Axis(f[2:3,1:4], xscale = log2, yscale = log10,
xticks = N_array, xlabel = "Points per dimension", ylabel = "Error")

scatter!(f[2:3,1:4], N_array, ERR[1, 2, :, 2], label = "Partial cells, L-2 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, ERR[1, 2, :, 2])
lines!(f[2:3,1:4], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

scatter!(f[2:3,1:4], N_array, ERR[2, 2, :, 2], label = "Full cells, L-2 norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, ERR[2, 2, :, 2])
lines!(f[2:3,1:4], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)

axislegend(position = :lb)

resize_to_layout!(f)

f = current_figure()
