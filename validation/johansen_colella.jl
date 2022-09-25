using Revise
using Flower

N_array = 2 .^[5, 6, 7]

plt2 = ("1", "2", "∞")

function conv_cutcell_johansen_colella(arr;
    case = "Sphere",
    L0 = 4.,
    T_inf = -1.0,
    R = 0.5)

    N = length(arr)
    ERR = zeros(3,N,2)
    V_data = Vector{Vector{Float64}}(undef, 0)
    angle_data = Vector{Vector{Float64}}(undef, 0)
    T_data = Vector{Matrix{Float64}}(undef, 0)
    STD = []
    for i in 1:N

        X = LinRange(-L0/2, L0/2, arr[i])
        Y = LinRange(-L0/2, L0/2, arr[i])

        num = Numerical(case = case,
            T_inf = T_inf,
            L0 = L0,
            x = X,
            y = Y,
            θd = 0.0,
            CFL = 0.1,
            max_iterations = 1,
            R = R,
            ϵ = 0.00);

        gp, gu, gv = init_meshes(num)
        opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)
        
        phL.T .= gp.u
        # phS.T .= gp.u
        @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
            opS, opL, phS, phL, fwd,
            stefan = true,
            verbose = true,
            show_every = 1
        )
        
        ERR[:, i, 1] .= normf(gp.V, MIXED, num.Δ)
        angle = similar(gp.V)
        for II in MIXED
            angle[II] = gp.geoL.projection[II].angle
        end
        push!(V_data, gp.V[MIXED])
        push!(angle_data, angle[MIXED])
        push!(T_data, phL.T)
        push!(STD, std(gp.V[MIXED] .- mean(gp.V[MIXED])))
        for i = 1:N
            ERR[i, :, 2] .= Richardson_extrapolation(ERR[i, :, 1], 2)
        end
    end
    return V_data, angle_data, T_data, ERR, STD
end

V_data, angle_data, T_data, ERR, STD = conv_cutcell_johansen_colella(N_array,
case = "Sphere",
R = 0.5096,
T_inf = -1.0)

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

f = Figure()
ax = Axis(f[1, 1], xticks = -180:90:180)

for i = 1:length(N_array)
    scatter!(f[1,1], angle_data[i]*180/pi, V_data[i], markersize = 5)
end

f2 = Figure(resolution = (800, 800))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)
ax = Axis(f2[1,1], xscale = log2, yscale = log10,
xticks = N_array, yticks = LogTicks(IntegerTicks()), xlabel = "Points per dimension", ylabel = "Error")
colsize!(f2.layout, 1, Aspect(1, 4/3))

scatter!(f2[1,1], N_array, ERR[2, :, 2], label = "L-1 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array, ERR[2, :, 2])
lines!(f2[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

scatter!(f2[1,1], N_array, ERR[3, :, 2], label = "L-inf norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array, ERR[3, :, 2])
lines!(f2[1,1], N_array, a*(1 ./ N_array).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)
axislegend(position = :lb)
resize_to_layout!(f2)

f2

#Makie.save("./figures/paper_figures/L2_Linf_VMIXED_n32_64_128_johansen_colella.png", f2)
