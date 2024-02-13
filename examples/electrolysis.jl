using Revise
using Flower

#From poisson.jl

fontsize_theme = Theme(fonts=(;regular="CMU Serif"), fontsize = 60)
set_theme!(fontsize_theme)

prefix="/local/home/pr277828/flower/test/"


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
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        D[II] = f(x_bc, y_bc)
    end
end

function neumann_bcs!(gp, N)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        N[II] = Nx * cos(LS[1].α[II]+π) + Ny * sin(LS[1].α[II]+π)
    end

    replace!(N, NaN=>0.0)

    return nothing
end

function robin_bcs!(gp, R)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        R[II] = Nx * cos(LS[1].α[II]+π) + Ny * sin(LS[1].α[II]+π) + f(x_bc, y_bc)
    end

    replace!(R, NaN=>0.0)

    return nothing
end

ϵ = 0.001
bc = rob
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
n = 512
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
    gp.LS[1].u .*= -1.0

    run_forward(num, gp, gu, gv, op, phS, phL, fwd, fwdS, fwdL)
    BC = Boundaries(left = per, bottom = per, right = per, top = per)
    BC_int = [bc]

    θd = zeros(gp)
    if is_dirichlet(bc)
        dirichlet_bcs!(gp, θd)
    elseif is_neumann(bc)
        neumann_bcs!(gp, θd)
    elseif is_robin(bc)
        robin_bcs!(gp, θd)
    end

    update_all_ls_data(num, gp, gu, gv, BC_int, true, true, false)
    laps = set_matrices!(num, gp, [gp.LS[1].geoL], gu, [gu.LS[1].geoL], gv, [gv.LS[1].geoL], op.opC_pL, op.opC_uL, op.opC_vL, true, true)
    Lp, bc_Lp, bc_Lp_b, Lu, bc_Lu, bc_Lu_b, Lv, bc_Lv, bc_Lv_b = laps
    a0_p = [θd]
    A, rhs = set_poisson(BC_int, num, gp, a0_p, op.opC_pL, op.opC_uL, op.opC_vL, 1, Lp, bc_Lp, bc_Lp_b, BC, true)

    b = Δf.(
        gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx,
        gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy
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
        gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx,
        gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy
    )

    for II in gp.ind.all_indices
        if gp.LS[1].geoL.cap[II,5] < 1e-12
            Tana[II] = 0.0
            T[II] = 0.0
        end
    end

    vlim = maximum(abs.(Tana))
    err = abs.(Tana .- T)

    LIQUID = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .> (1-1e-16)]
    MIXED = gp.ind.all_indices[gp.LS[1].geoL.cap[:,:,5] .<= (1-1e-16) .&& gp.LS[1].geoL.cap[:,:,5] .> 1e-16]

    println("$(length(MIXED) * 100 / (length(LIQUID) + length(MIXED)))% of mixed cells")

    norm_all = normf(err, vcat(LIQUID, MIXED), gp.LS[1].geoL.cap[:,:,5], num.Δ)
    println("ALL: $norm_all")
    norm_mixed = normf(err, MIXED, gp.LS[1].geoL.cap[:,:,5], num.Δ)
    println("MIXED: $norm_mixed")
    norm_full = normf(err, LIQUID, gp.LS[1].geoL.cap[:,:,5], num.Δ)
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


tcks = -num.L0:0.5:num.L0

fa = Figure(figure_padding=(0, 50, 50, 50), size = (1600, 1000))
ax = Axis(fa[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], Tana', colorrange=(-vlim, vlim))
cbar = fa[1,2] = Colorbar(fa, hmap, labelpadding=0)
colsize!(fa.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fa.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fa)

ft = Figure(figure_padding=(0, 50, 50, 50), size = (1000, 1000))
ax = Axis(ft[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], T', colorrange=(-vlim, vlim))
cbar = ft[1,2] = Colorbar(ft, hmap, labelpadding=0)
colsize!(ft.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(ft.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(ft)

fd = Figure(figure_padding=(0, 50, 50, 50), size = (1600, 1000))
ax = Axis(fd[1,1], aspect = 1, xlabel=L"x", ylabel=L"y", xticks = tcks, yticks = tcks,
    xgridvisible=false, ygridvisible=false)
hmap = heatmap!(gp.x[1,:], gp.y[:,1], D', colorrange=(-vlim, vlim))
cbar = fd[1,2] = Colorbar(fd, hmap, labelpadding=0)
colsize!(fd.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fd.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fd)

x_cen = getproperty.(gp.LS[1].mid_point[gp.LS[1].MIXED], :x)
y_cen = getproperty.(gp.LS[1].mid_point[gp.LS[1].MIXED], :y)
x_mix = gp.x[gp.LS[1].MIXED] .+ x_cen .* gp.dx[gp.LS[1].MIXED]
y_mix = gp.y[gp.LS[1].MIXED] .+ y_cen .* gp.dy[gp.LS[1].MIXED]
θ = atan.(y_mix, x_mix) .* 180 ./ pi

perm = sortperm(θ)
θp = θ[perm]
Dp = D[gp.LS[1].MIXED][perm]

fd2 = Figure(size = (1600, 1000))
ax = Axis(fd2[1,1], aspect=1, xlabel=L"\theta", ylabel=L"p ^ \gamma", xticks = -180:45:180)
lines!(θp[.!isnan.(Dp)], Dp[.!isnan.(Dp)], linewidth = 3)
limits!(ax, -180, -135, -0.6, 0.6)
colsize!(fd2.layout, 1, widths(ax.scene.viewport[])[1])
rowsize!(fd2.layout, 1, widths(ax.scene.viewport[])[2])
resize_to_layout!(fd2)

Makie.save(prefix*"f_ana.pdf", fa)
Makie.save(prefix*"f_bulk.pdf", ft)
Makie.save(prefix*"f_bound.pdf", fd2)
#Makie.save(prefix*"f_err.pdf", fe)

nothing

