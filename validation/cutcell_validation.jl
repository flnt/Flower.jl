using Revise
using Flower

N_array = 2 .^[4, 5, 6, 7, 8]
plt = ("partial", "full", "all")
plt2 = ("1", "2", "∞")

function conv_cutcell_CN(arr;
    case = "Sphere",
    L0 = 2.,
    Θd = 1.0,
    ϵ_κ = 0.0,
    ii = 4,
    R = 0.75,
    solid = true,
    liquid = false,
    verbose = true)

    N = length(arr)
    ERR = zeros(3,3,N)
    T_data = Vector{Matrix{Float64}}(undef, 0)
    e_data = Vector{Matrix{Float64}}(undef, 0)
    u_data = Vector{Matrix{Float64}}(undef, 0)
    grid = Vector{Vector{Float64}}(undef, 0)

    X = LinRange(-L0/2, L0/2, arr[N]+1)
    Y = LinRange(-L0/2, L0/2, arr[N]+1)
    num = Numerical(case = case,
        CFL = 0.5,
        L0 = L0,
        θd = Θd,
        ϵ_κ = ϵ_κ,
        save_every = 1,
        x = X,
        y = Y,
        max_iterations = ii*4^(N),
        R = R,
        ϵ = 0.0,  
        A = -0.2,
        N = 6);

    gp, gu, gv = init_meshes(num)
    opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

    # if case == "Planar"
    #     gp.u .+= num.Δ/4
    # end
    @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
        opS, opL, phS, phL, fwd,
        heat = true,
        heat_solid_phase = solid,
        heat_liquid_phase = liquid,
        verbose = verbose,
        show_every = 1,
    )

    if solid 
        base = phS.T
        cap = gp.geoS.cap[:,:,5]
    elseif liquid
        base = phL.T
        cap = gp.geoL.cap[:,:,5]
    end

    for i in 1:N-1
        X = LinRange(-L0/2, L0/2, arr[i]+1)
        Y = LinRange(-L0/2, L0/2, arr[i]+1)
        num = Numerical(case = case,
            CFL = 0.5,
            L0 = L0,
            θd = Θd,
            ϵ_κ = ϵ_κ,
            save_every = 1,
            x = X,
            y = Y,
            max_iterations = ii*4^(i),
            R = R,
            ϵ = 0.0,
            A = -0.2,
            N = 6);

        gp, gu, gv = init_meshes(num)
        opS, opL, phS, phL, fwd = init_fields(num, gp, gu, gv)

        @time MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv,
            opS, opL, phS, phL, fwd,
            heat = true,
            heat_solid_phase = solid,
            heat_liquid_phase = liquid,
            verbose = verbose,
            show_every = 1,
        )
        
        @show (num.max_iterations*num.τ)
        @show (num.τ)
        if solid
            interpolated_field = interpolate_field(base, phS.T, cap)
            e = interpolated_field - phS.T.* gp.geoS.cap[:,:,5]
            ERR[1, :, i] .= normf(e, MIXED, gp.geoS.cap[:,:,5], num.Δ)
            ERR[2, :, i] .= normf(e, SOLID, gp.geoS.cap[:,:,5], num.Δ)
            ERR[3, :, i] .= normf(e, vcat(SOLID, MIXED), gp.geoS.cap[:,:,5], num.Δ)
            push!(T_data, phS.T)
            push!(e_data, abs.(e))
        elseif liquid
            interpolated_field = interpolate_field(base, phL.T, cap)
            e = interpolated_field - phL.T .* gp.geoL.cap[:,:,5]
            ERR[1, :, i] .= normf(e, MIXED, gp.geoL.cap[:,:,5], num.Δ)
            ERR[2, :, i] .= normf(e, LIQUID, gp.geoL.cap[:,:,5], num.Δ)
            ERR[3, :, i] .= normf(e, vcat(LIQUID, MIXED), gp.geoL.cap[:,:,5], num.Δ)
            push!(T_data, phL.T)
            push!(e_data, abs.(e))
        end
        push!(u_data, gp.u)
        push!(grid, gp.x[1,:])
    end
    return T_data, u_data, e_data, grid, ERR
end

case = "Square"
R = 0.8
ii = 1
Θd = 1.0
ϵ_κ = 0.0
liquid = 0
solid = 1
verbose = true

T_data, u_data, e_data, grid, ERR = conv_cutcell_CN(N_array,
case = case,
R = R,
ii = ii,
Θd = Θd,
ϵ_κ = ϵ_κ,
liquid = Bool(liquid),
solid = Bool(solid),
verbose = verbose)

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

f = Figure(resolution = (1200, 1200))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

for i in 1:length(N_array)-1
    ax = Axis(f[1, i], title = "N = $(N_array[i])")
    colsize!(f.layout, i, Aspect(1, 1.0))
    hidedecorations!(ax)
    if case == "Crystal"
        hm = heatmap!(f[1,i], grid[i], grid[i], T_data[i]', colormap = Reverse(:redsblues))
    else
        hm = heatmap!(f[1,i], grid[i], grid[i], T_data[i]', colormap = Reverse(:ice))
    end
    contour!(f[1,i], grid[i], grid[i], u_data[i]', levels = 0:0, linewidth = 2, color =:red)
    ax = Axis(f[2, i])
    hidedecorations!(ax)
    hm2 = heatmap!(f[2,i], grid[i], grid[i], e_data[i]', colorrange = (0, max(e_data[i]...)))#
end

if case == "Crystal"
    Colorbar(f[1, length(N_array)], vertical = true, flipaxis = true, label = "Temperature", ticks = [-1, 1], colorrange = (-1, 1), colormap = Reverse(:redsblues))
else
    Colorbar(f[1, length(N_array)], vertical = true, flipaxis = true, label = "Temperature", ticks = [0, 1], colorrange = (0, 1), colormap = Reverse(:ice))
end
Colorbar(f[2, length(N_array)], vertical = true, flipaxis = true, label = "Normalized error", ticks = 0:1)

ax = Axis(f[3:length(N_array),1:length(N_array)], xscale = log2, yscale = log10,
xticks = N_array, yticks = LogTicks(IntegerTicks()), xlabel = "Points per dimension", ylabel = "Error")

XX = ERR[1, 2, 1:end-1]

scatter!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], XX, label = "Mixed cells, L-2 norm", markersize = 25, color =:black, marker=:rect)
a, b = fit_order(N_array[1:end-1], XX)
lines!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], a*(1 ./ N_array[1:end-1]).^b, linewidth = 3, linestyle =:dash, color =:black, label = @sprintf "order %.2f" b)

YY = ERR[2, 2, 1:end-1]

scatter!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], YY, label = "Full cells, L-2 norm", markersize = 25, color =:black, marker=:utriangle)
a, b = fit_order(N_array[1:end-1], YY)
lines!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], a*(1 ./ N_array[1:end-1]).^b, linewidth = 3, linestyle =:dashdot, color =:black, label = @sprintf "order %.2f" b)

ZZ = ERR[3, 2, 1:end-1]

scatter!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], ZZ, label = "All cells, L-2 norm", markersize = 25, color =:black, marker=:circle)
a, b = fit_order(N_array[1:end-1], ZZ)
lines!(f[3:length(N_array),1:length(N_array)], N_array[1:end-1], a*(1 ./ N_array[1:end-1]).^b, linewidth = 3, linestyle =:dot, color =:black, label = @sprintf "order %.2f" b)


axislegend(position = :lb)

resize_to_layout!(f)

f = current_figure()

Makie.save("./figures/validation/cutcell_$(case)_R$(R)_T$(Θd)_k$(ϵ_κ)_ite$(ii)_L$(liquid)_S$(solid).png", f)
