using JLD2
using Flower

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_2.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+03.jld2"

num_Ra1000 = num
gp_Ra1000 = gp
gu_Ra1000 = gu
gv_Ra1000 = gv
phS_Ra1000 = phS
phL_Ra1000 = phL
fwd_Ra1000 = fwd
av_height1_Ra1000 = av_height1
eRa_Ra1000 = eRa

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_2.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+04.jld2"

num_Ra10000 = num
gp_Ra10000 = gp
gu_Ra10000 = gu
gv_Ra10000 = gv
phS_Ra10000 = phS
phL_Ra10000 = phL
fwd_Ra10000 = fwd
av_height1_Ra10000 = av_height1
eRa_Ra10000 = eRa

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_2.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+05.jld2"

num_Ra100000 = num
gp_Ra100000 = gp
gu_Ra100000 = gu
gv_Ra100000 = gv
phS_Ra100000 = phS
phL_Ra100000 = phL
fwd_Ra100000 = fwd
av_height1_Ra100000 = av_height1
eRa_Ra100000 = eRa

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_4.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+05.jld2"

num_Ra100000_long = num
gp_Ra100000_long = gp
gu_Ra100000_long = gu
gv_Ra100000_long = gv
phS_Ra100000_long = phS
phL_Ra100000_long = phL
fwd_Ra100000_long = fwd
av_height1_Ra100000_long = av_height1
eRa_Ra100000_long = eRa

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_2.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_2.0e+05.jld2"

num_Ra200000 = num
gp_Ra200000 = gp
gu_Ra200000 = gu
gv_Ra200000 = gv
phS_Ra200000 = phS
phL_Ra200000 = phL
fwd_Ra200000 = fwd
av_height1_Ra200000 = av_height1
eRa_Ra200000 = eRa

tcks = -ratio*L0/2:2:ratio*L0/2
lim = L0 / 2
fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, ratio))
for i = 1:size(fwd_Ra100000_long.TLsave,1)
    ax = Axis(fp[1,1], aspect = ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
    hidedecorations!(ax)
    heatmap!(gp_Ra100000_long.y[:,1], gp_Ra100000_long.x[1,:], fwd_Ra100000_long.Tsave[i,:,:], colormap = :thermal)#, colorrange=(num.θd, max(fwd_Ra100000_long.Tsave[5,:,:]...)))
    contour!(gp_Ra100000_long.y[:,1], gp_Ra100000_long.x[1,:], fwd_Ra100000_long.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
    resize_to_layout!(fp)
    Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", 100000))-i$(i).png", fp)
end

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)


fh = Figure(resolution = (1600, 1000))
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Dimensionless height")
lines!(fwd_Ra1000.time, av_height1_Ra1000, linewidth = 3, label = "Ra = 10³")
lines!(fwd_Ra10000.time, av_height1_Ra10000, linewidth = 3, label = "Ra = 10⁴")
lines!(fwd_Ra100000_long.time, av_height1_Ra100000, linewidth = 3, label = "Ra = 10⁵")
lines!(fwd_Ra200000.time, av_height1_Ra200000, linewidth = 3, label = "Ra = 2×10⁵")
resize_to_layout!(fh)
axislegend(position = :rb)
fu = current_figure()
# Makie.save("/home/tf/All/Thesis/Figures/Chapter5/height_Ra_TM0.png", fu)

feRa = Figure(resolution = (1600, 1000))
colsize!(feRa.layout, 1, Aspect(1, 1.0))
ax = Axis(feRa[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Effective Rayleigh number",
    yscale = log10, yticks = LogTicks(IntegerTicks()))
# ax.yticks = [1e-6]
lines!(fwd_Ra1000.time[2:end], eRa_Ra1000[2:end], linewidth = 3, label = "Ra = 10³")
lines!(fwd_Ra10000.time[2:end], eRa_Ra10000[2:end], linewidth = 3, label = "Ra = 10⁴")
lines!(fwd_Ra100000_long.time[2:end], eRa_Ra100000[2:end], linewidth = 3, label = "Ra = 10⁵")
lines!(fwd_Ra200000.time[2:end], eRa_Ra200000[2:end], linewidth = 3, label = "Ra = 2×10⁵")
hlines!(ax, [1707.76], color = :black, linestyle =:dash, linewidth = 3, label = "Raₑ = 1707.76")
resize_to_layout!(feRa)
axislegend(position = :rb)
feRa = current_figure()
# Makie.save("/home/tf/All/Thesis/Figures/Chapter5/eRa_TM0.png", feRa)

