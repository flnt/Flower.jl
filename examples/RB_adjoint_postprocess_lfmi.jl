using JLD2
using Flower
using Interpolations

Ratio = 4
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

JLD2.@load "/home/tf/Documents/RB_opt/opt_data_lfmi_new.jld2" num gp gu gv fwd_des opt_p opt_S opt_L opt_u opt_uu opt_RB res

store = zeros(length(res.trace), 2)
for i in axes(store,1)
    store[i, 1] = res.trace[i].iteration
    store[i, 2] = res.trace[i].value
end

f = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f[1,1], yscale = log10, xlabel = "Iteration", ylabel = L" J / J_0")

lines!(f[1,1], store[:,1], store[:,2]./store[1,2], color =:black, linewidth = 3)
scatter!(f[1,1], store[:,1], store[:,2]./store[1,2], markersize = 10, color =:black, marker=:rect)

f = current_figure()
Makie.save("./figures/paper_figures/RB_opt_costalot.png", f)

store_p = zeros(length(opt_p), 10)
for i in axes(store,1)
    store_p[i, 1] = opt_p[i][1]
    store_p[i, 2] = opt_p[i][2]
    store_p[i, 3] = opt_p[i][3]
    store_p[i, 4] = opt_p[i][4]
    store_p[i, 5] = opt_p[i][5]
    store_p[i, 6] = opt_p[i][6]
    store_p[i, 7] = opt_p[i][7]
    store_p[i, 8] = opt_p[i][8]
    store_p[i, 9] = opt_p[i][9]
    store_p[i, 10] = opt_p[i][10]
end

f2 = Figure()
fontsize_theme = Theme(fontsize = 20)
set_theme!(fontsize_theme)
ax = Axis(f2[1,1], xlabel = "Iteration", ylabel = "Actuator coefficients")

lines!(f2[1,1], store[:,1],  store_p[:, 1], linewidth = 3, label = L"a_1")
scatter!(f2[1,1], store[:,1], store_p[:, 1], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 2], linewidth = 3, label = L"a_2")
scatter!(f2[1,1], store[:,1], store_p[:, 2], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 3], linewidth = 3, label = L"a_3")
scatter!(f2[1,1], store[:,1], store_p[:, 3], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 4], linewidth = 3, label = L"a_4")
scatter!(f2[1,1], store[:,1], store_p[:, 5], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 6], linewidth = 3, label = L"a_5")
scatter!(f2[1,1], store[:,1], store_p[:, 6], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 7], linewidth = 3, label = L"a_6")
scatter!(f2[1,1], store[:,1], store_p[:, 7], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 8], linewidth = 3, label = L"a_7")
scatter!(f2[1,1], store[:,1], store_p[:, 8], markersize = 10, marker=:rect)
lines!(f2[1,1], store[:,1],  store_p[:, 9], linewidth = 3, label = L"a_8")
scatter!(f2[1,1], store[:,1], store_p[:, 10], markersize = 10, marker=:rect)
axislegend(position = :rc)
f2 = current_figure()
Makie.save("./figures/paper_figures/RB_opt_costalot_coeff.png", f2)