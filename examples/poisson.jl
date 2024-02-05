using Revise
using Flower
using IncompleteLU

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 40)
set_theme!(fontsize_theme)

prefix = "/Users/alex/Documents/PhD/papers/cutcell_NS_stationary/CPC/figures/"

function f(x, y)
    return cos(π^2 * x * y)*sin(π^2 * x * y)
end

function ∇fx(x, y)
    return π^2 * y * (cos(π^2 * x * y)^2 - sin(π^2 * x * y)^2)
end

function ∇fy(x, y)
    return π^2 * x * (cos(π^2 * x * y)^2 - sin(π^2 * x * y)^2)
end

function Δf(x, y)
    return -4 * (π^2)^2 *sin(π^2 * x*y)*cos(π^2 * x*y)*(x^2+y^2)
end

mutable struct conv_regression
    yreg::Array{Float64, 1}
    coef0::Float64
    coef1::Float64
end

function regression(x, y, x_reg)
    coeffs = fit(log.(x), log.(y), 1).coeffs
    y_reg = exp.(coeffs[2].*log.(x_reg) .+ coeffs[1])
    result = conv_regression(y_reg, -coeffs[1], -coeffs[2])

    return result
end

function dirichlet_bcs!(gp, D)
    @unpack x, y, dx, dy, mid_point, ind, α = gp

    @inbounds @threads for II in ind.inside
        x_bc = mid_point[II].x * dx[II] + x[II]
        y_bc = mid_point[II].y * dy[II] + y[II]
        D[II] = f(x_bc, y_bc)
    end
end

function neumann_bcs!(gp, N)
    @unpack x, y, dx, dy, mid_point, ind, α = gp

    @inbounds @threads for II in ind.inside
        x_bc = mid_point[II].x * dx[II] + x[II]
        y_bc = mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        N[II] = Nx * cos(α[II]+π) + Ny * sin(α[II]+π)
    end

    replace!(N, NaN=>0.0)

    return nothing
end

function robin_bcs!(gp, R)
    @unpack x, y, dx, dy, mid_point, ind, α = gp

    @inbounds @threads for II in ind.inside
        x_bc = mid_point[II].x * dx[II] + x[II]
        y_bc = mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        R[II] = Nx * cos(α[II]+π) + Ny * sin(α[II]+π) + f(x_bc, y_bc)
    end

    replace!(R, NaN=>0.0)

    return nothing
end

ϵ = 0.01
bc = neu
if is_dirichlet(bc)
    bc_str = "dir"
elseif is_neumann(bc)
    bc_str = "neu"
elseif is_robin(bc)
    bc_str = "rob"
end

L0 = 2.5

n_cases = 4
l1 = zeros(n_cases)
l2 = zeros(n_cases)
loo = zeros(n_cases)
l1_mixed = zeros(n_cases)
l2_mixed = zeros(n_cases)
loo_mixed = zeros(n_cases)
l1_full = zeros(n_cases)
l2_full = zeros(n_cases)
loo_full = zeros(n_cases)

st_case = 6

npts = 2 .^ [st_case:(st_case+n_cases-1)...]

i = 1
n = 32
# for (i,n) in enumerate(npts)
    x = LinRange(-L0 / 2.0, L0 / 2.0, n + 1)
    y = LinRange(-L0 / 2.0, L0 / 2.0, n + 1)

    num = Numerical(
        case = "Cylinder",
        x = x,
        y = y,
        CFL = 1.0,
        max_iterations = 0,
        R = 1.0,
        ϵ = ϵ,
    )

    gp, gu, gv = init_meshes(num)
    op, phS, phL, fwd, fwdS, fwdL = init_fields(num, gp, gu, gv)
    gp.u .*= -1.0

    MIXED, SOLID, LIQUID = run_forward(num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL)
    BC = Boundaries(left = per, bottom = per, right = per, top = per)

    θd = zeros(gp)
    if is_dirichlet(bc)
        dirichlet_bcs!(gp, θd)
    elseif is_neumann(bc)
        neumann_bcs!(gp, θd)
    elseif is_robin(bc)
        robin_bcs!(gp, θd)
    end

    update_ls_data(num, gp, gu, gv, gp.u, gp.κ, true, true, false)
    laps = set_matrices!(gp, gp.geoL, gu, gu.geoL, gv, gv.geoL, op.opC_pL, op.opC_uL, op.opC_vL, true, true)
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps
    A, rhs, data_A = set_poisson(bc, num, gp, θd, op.opC_pL, op.opC_uL, op.opC_vL, Lp, bc_Lp, bc_Lp_b, BC)

    b = Δf.(
        gp.x .+ getproperty.(gp.geoL.centroid, :x) .* gp.dx,
        gp.y .+ getproperty.(gp.geoL.centroid, :y) .* gp.dy
    )
    veci(rhs,gp,1) .+= op.opC_pL.M * vec(b)

    res = zeros(size(rhs))

    @time @inbounds @threads for i in 1:A.m
        @inbounds A[i,i] += 1e-10
    end
    @time res .= A \ rhs

    T = reshape(vec1(res, gp), gp)
    D = reshape(vec2(res, gp), gp)
    Tana = f.(
        gp.x .+ getproperty.(gp.geoL.centroid, :x) .* gp.dx,
        gp.y .+ getproperty.(gp.geoL.centroid, :y) .* gp.dy
    )

    for II in gp.ind.all_indices
        if gp.geoL.cap[II,5] < 1e-12
            Tana[II] = 0.0
            T[II] = 0.0
        end
    end

    vlim = maximum(abs.(Tana))
    err = abs.(Tana .- T)

    LIQUID = gp.ind.all_indices[gp.geoL.cap[:,:,5] .> (1-1e-16)]
    MIXED = gp.ind.all_indices[gp.geoL.cap[:,:,5] .<= (1-1e-16) .&& gp.geoL.cap[:,:,5] .> 1e-16]

    println("$(length(MIXED) * 100 / (length(LIQUID) + length(MIXED)))% of mixed cells")

    norm_all = normf(err, vcat(LIQUID, MIXED), gp.geoL.cap[:,:,5], num.Δ)
    println("ALL: $norm_all")
    norm_mixed = normf(err, MIXED, gp.geoL.cap[:,:,5], num.Δ)
    println("MIXED: $norm_mixed")
    norm_full = normf(err, LIQUID, gp.geoL.cap[:,:,5], num.Δ)
    println("FULL: $norm_full")

    l1[i] = norm_all[1]
    l2[i] = norm_all[2]
    loo[i] = norm_all[3]

    l1_mixed[i] = norm_mixed[1]
    l2_mixed[i] = norm_mixed[2]
    loo_mixed[i] = norm_mixed[3]

    l1_full[i] = norm_full[1]
    l2_full[i] = norm_full[2]
    loo_full[i] = norm_full[3]

    for II in gp.ind.all_indices
        if abs.(Tana[II]) < 1e-16
            Tana[II] = NaN
            T[II] = NaN
            err[II] = NaN
        end
        if abs.(D[II]) < 1e-16
            D[II] = NaN
        end
    end
# end

# x_reg = 2^st_case:2^(st_case+n_cases-1)
# y_tcks = [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

# conv_l1 = regression(npts, l1, x_reg)
# conv_l2 = regression(npts, l2, x_reg)
# conv_loo = regression(npts, loo, x_reg)

# conv_l1_mixed = regression(npts, l1_mixed, x_reg)
# conv_l2_mixed = regression(npts, l2_mixed, x_reg)
# conv_loo_mixed = regression(npts, loo_mixed, x_reg)

# conv_l1_full = regression(npts, l1_full, x_reg)
# conv_l2_full = regression(npts, l2_full, x_reg)
# conv_loo_full = regression(npts, loo_full, x_reg)

# fl1 = Figure(resolution = (1600, 1000))
# ax = Axis(fl1[1,1], aspect = 0.5, xscale=log10, yscale=log10,
#             xlabel="pts", ylabel=L"$L _ 1$ Error", xticks = npts
# )
# colsize!(fl1.layout, 1, Aspect(1, 0.5))
# fl1 = lines!(x_reg, conv_l1_mixed.yreg, label="order $(@sprintf("%.2f", conv_l1_mixed.coef1))", linewidth=6.0)
# fl1 = scatter!(npts, l1_mixed, label="Partial cells", markersize=40)
# fl1 = lines!(x_reg, conv_l1.yreg, label="order $(@sprintf("%.2f", conv_l1.coef1))", linewidth=6.0)
# fl1 = scatter!(npts, l1, marker=:rect, label="All cells", markersize=40)
# fl1 = lines!(x_reg, conv_l1_full.yreg, label="order $(@sprintf("%.2f", conv_l1_full.coef1))", linewidth=6.0)
# fl1 = scatter!(npts, l1_full, marker=:diamond, label="Full cells", markersize=40)
# axislegend(position = :lb)
# fl1 = current_figure()
# resize_to_layout!(fl1)
# Makie.save(prefix*"conv_l1_"*bc_str*"_eps$(ϵ).pdf", fl1)

# fl2 = Figure(figure_padding=(0, 50, 50, 50), resolution = (1600, 1000))
# ax = Axis(fl2[1,1], aspect = 0.5, xscale=log10, yscale=log10,
#             xlabel="pts", ylabel=L"$L _ 2$ error", xticks = npts
# )
# colsize!(fl2.layout, 1, Aspect(1, 0.5))
# lines!(x_reg, conv_l2_mixed.yreg, label="order $(@sprintf("%.2f", conv_l2_mixed.coef1))", linewidth=6.0)
# scatter!(npts, l2_mixed, label="Partial cells", markersize=40)
# lines!(x_reg, conv_l2.yreg, label="order $(@sprintf("%.2f", conv_l2.coef1))", linewidth=6.0)
# scatter!(npts, l2, marker=:rect, label="All cells", markersize=40)
# lines!(x_reg, conv_l2_full.yreg, label="order $(@sprintf("%.2f", conv_l2_full.coef1))", linewidth=6.0)
# scatter!(npts, l2_full, marker=:diamond, label="Full cells", markersize=40)
# # fl2[1,1] = Legend(fl2, ax, framevisible = false)
# axislegend(position = :lb)
# fl2 = current_figure()
# resize_to_layout!(fl2)
# Makie.save(prefix*"conv_l2_"*bc_str*"_eps$(ϵ).pdf", fl2)

# floo = Figure(figure_padding=(0, 50, 50, 50), resolution = (1600, 1000))
# ax = Axis(floo[1,1], aspect = 0.5, xscale=log10, yscale=log10,
#             xlabel="pts", ylabel=L"$L _ \infty$ error", xticks = npts
# )
# colsize!(floo.layout, 1, Aspect(1, 0.5))
# floo = lines!(x_reg, conv_loo_mixed.yreg, label="order $(@sprintf("%.2f", conv_loo_mixed.coef1))", linewidth=6.0)
# floo = scatter!(npts, loo_mixed, label="Partial cells", markersize=40)
# floo = lines!(x_reg, conv_loo.yreg, label="order $(@sprintf("%.2f", conv_loo.coef1))", linewidth=6.0)
# floo = scatter!(npts, loo, marker=:rect, label="All cells", markersize=40)
# floo = lines!(x_reg, conv_loo_full.yreg, label="order $(@sprintf("%.2f", conv_loo_full.coef1))", linewidth=6.0)
# floo = scatter!(npts, loo_full, marker=:diamond, label="Full cells", markersize=40)
# axislegend(position = :lb)
# floo = current_figure()
# resize_to_layout!(floo)
# Makie.save(prefix*"conv_loo_"*bc_str*"_eps$(ϵ).pdf", floo)

tcks = -num.L0:0.5:num.L0

fa = Figure(figure_padding=(0, 50, 50, 50), resolution = (1600, 1000))
colsize!(fa.layout, 1, Aspect(1, 1.0))
ax = Axis(fa[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], Tana', colorrange=(-vlim, vlim))
cbar = fa[1,2] = Colorbar(fa, hmap, labelpadding=0)
resize_to_layout!(fa)

ft = Figure(figure_padding=(0, 50, 50, 50), resolution = (1600, 1000))
colsize!(ft.layout, 1, Aspect(1, 1.0))
ax = Axis(ft[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], T', colorrange=(-vlim, vlim))
cbar = ft[1,2] = Colorbar(ft, hmap, labelpadding=0)
resize_to_layout!(ft)

fe = Figure(figure_padding=(0, 50, 50, 50), resolution = (1600, 1100))
colsize!(fe.layout, 1, Aspect(1, 1.0))
ax = Axis(fe[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], log10.(err)', colorrange=(-5.5, -1.0))
cbar = fe[1,2] = Colorbar(fe, hmap, labelpadding=0;
    tickformat=custom_formatter,
    minorticksvisible=true,
    minorticks=LogMinorTicks(),
)
resize_to_layout!(fe)

# Makie.save(prefix*"f_ana.pdf", fa)
# Makie.save(prefix*"f_cut.pdf", ft)
# Makie.save(prefix*"f_err.pdf", fe)

nothing

