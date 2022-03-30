using Revise
using Flower

NB_array = [0, 6, 12]

function test_NB(NB_array; n = 128, ϵ_κ = 0.005, R = 0.5, case = "Crystal")
    N = length(NB_array)
    u_data = Vector{Matrix{Float64}}(undef, 0)
    V_data = Vector{Matrix{Float64}}(undef, 0)
    for i in 1:N
        num = Numerical(ϵ_κ = ϵ_κ,
            R = R,
            case = case,
            n = n,
            max_iterations = 1,
            NB = NB_array[i],
            A = 0.2,
            N = 6
            )

        idx = set_indices(num.n)
        tmp, fwd = init_fields(num, idx)
        MIXED = run_forward(num, idx, tmp, fwd,
            stefan = true)
        for II in CartesianIndices(fwd.V)
            if fwd.V[II] == 0. fwd.V[II] = NaN end
        end
        push!(u_data, fwd.u)
        push!(V_data, fwd.V)
    end
    return u_data, V_data
end

u_data, V_data = test_NB(NB_array, n = 256, R = 0.5)

f = Figure()
fontsize_theme = Theme(fontsize = 35)
set_theme!(fontsize_theme)

for i in 1:length(NB_array)
    ax = Axis(f[1, i], title = "NB = $(NB_array[i])")
    colsize!(f.layout, i, Aspect(1, 1.0))
    hidedecorations!(ax)
    heatmap!(f[1,i], V_data[i]')
    if i != 1
        contour!(f[1,i], u_data[i]', levels = 0:0, linewidth = 2, color =(:red, 0.5))
    end
end
Colorbar(f[1, 4],label ="Normalized V", limits = (-1, 1))

resize_to_layout!(f)
f = current_figure()
