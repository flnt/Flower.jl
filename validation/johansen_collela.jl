using Revise
using Flower

N_array = 2 .^[5, 6, 7]

plt2 = ("1", "2", "∞")

function conv_cutcell_johansen_colella(x;
    case = "Sphere",
    L0 = 4.,
    T_inf = -1.0,
    R = 0.5)

    N = length(x)
    r = (x[end]+1)/(x[end-1]+1)
    ERR = zeros(3,N,2)
    V_data = Vector{Vector{Float64}}(undef, 0)
    angle_data = Vector{Vector{Float64}}(undef, 0)
    T_data = Vector{Matrix{Float64}}(undef, 0)

    for i in 1:N
        num = Numerical(case = case,
            T_inf = T_inf,
            L0 = L0,
            n = x[i],
            CFL = 0.1,
            max_iterations = 1,
            R = R);

        idx = set_indices(num.n);
        tmp, fwd = init_fields(num, idx);
        fwd.TL .= fwd.u
        MIXED, SOLID, LIQUID = run_forward(num, idx, tmp, fwd,
        stefan = true,
        verbose = false
        )
        ERR[:, i, 1] .= normf(fwd.V, MIXED, num.Δ)
        angle = similar(fwd.V)
        for II in MIXED
            angle[II] = tmp.liq_projection[II].angle
        end
        push!(V_data, fwd.V[MIXED])
        push!(angle_data, angle[MIXED])
        push!(T_data, fwd.TL)
    end
    for i = 1:3
        ERR[i, :, 2] .= Richardson_extrapolation(ERR[i, :, 1], r)
    end
    return V_data, angle_data, T_data, ERR
end

V_data, angle_data, T_data, ERR = conv_cutcell_johansen_colella(N_array,
case = "Sphere",
R = 0.5096,
T_inf = -1.0)

f = Figure()
ax = Axis(f[1, 1], xticks = -180:90:180)

for i = 1:length(N_array)
    scatter!(f[1,1], angle_data[i]*180/pi, V_data[i], markersize = 5)
    #lines!(f[1,1], angle_data[i]*180/pi, ERR[1,i,1]*ones(length(V_data[i])), markersize = 5)
end

f2 = Figure(resolution = (800, 800))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
ax = Axis(f2[1,1], xscale = log2, yscale = log10,
xticks = N_array, xlabel = "Points per dimension", ylabel = "Error")
#ylims!(2^-8, 2^-4)
colsize!(f2.layout, 1, Aspect(1, 4/3))
#=for i in (2,3) #type of norm (1, 2, inf)
    scatter!(f2,N_array, ERR[i, :, 2], label = "L-$(plt2[i]) norm", markersize = 10)
    a, b = fit_order(N_array, ERR[i, :, 2])
    lines!(f2,N_array, a*(1 ./ N_array).^b, label = @sprintf "order %.2f" b)
end=#
scatter!(f2[1,1], N_array, ERR[1, :, 2], label = "L-1 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, ERR[1, :, 2])
lines!(f2[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

#=scatter!(f2[1,1], N_array, ERR[2, :, 2], label = "L-2 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, ERR[2, :, 2])
lines!(f2[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 2, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)
=#
scatter!(f2[1,1], N_array, ERR[3, :, 2], label = "L-inf norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, ERR[3, :, 2])
lines!(f2[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)
axislegend(position = :lb)
resize_to_layout!(f2)
rm("fit.log")
f2 = current_figure()

#Makie.save("./figures/paper_figures/L2_Linf_VMIXED_n32_64_128_johansen_colella.png", f2)
