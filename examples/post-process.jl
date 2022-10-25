using JLD2
using Flower

Ratio = 8
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

# Ra = 10^3 iter 2e4
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




# Ra 10^4 iter 2e4
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




# Ra 10^4 iter 4e4
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_4.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+04.jld2"

height = zeros(size(fwd.usave)[1:2])

for i in 1:size(fwd.usave,1)
    for j in 1:size(fwd.usave,2)
        itp = LinearInterpolation(reverse(fwd.usave[i,j,2:end-1]), reverse(gp.x[1,2:end-1]))
        height[i,j] = itp(0.0)
    end
end

av_height1 = mean(height, dims=2)[:,1] .+ 0.5 .- 0.05

eRa = av_height1.^3 .* 1e4 .* (1.0 .- num.θd)

num_Ra10000_long = num
gp_Ra10000_long = gp
gu_Ra10000_long = gu
gv_Ra10000_long = gv
phS_Ra10000_long = phS
phL_Ra10000_long = phL
fwd_Ra10000_long = fwd
av_height1_Ra10000_long = radius
eRa_Ra10000_long = eRa




# Ra 5x10^4 iter 4e4
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_64-ny_512-ratio_8-maxiter_4.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_5.0e+04.jld2"

num_Ra50000_long = num
gp_Ra50000_long = gp
gu_Ra50000_long = gu
gv_Ra50000_long = gv
phS_Ra50000_long = phS
phL_Ra50000_long = phL
fwd_Ra50000_long = fwd
av_height1_Ra50000_long = av_height1
eRa_Ra50000_long = eRa





# Ra 10^5 iter 2e4
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




# Ra 10^5 iter 4e4
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




# Ra 2x10^5 iter 2e4
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


# Ra 5x10^5 iter 4e4
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_80-ny_640-ratio_8-maxiter_4.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_5.0e+05.jld2"

num_Ra500000 = num
gp_Ra500000 = gp
gu_Ra500000 = gu
gv_Ra500000 = gv
phS_Ra500000 = phS
phL_Ra500000 = phL
fwd_Ra500000 = fwd
av_height1_Ra500000 = av_height1



height = zeros(size(fwd.usave)[1:2])

for i in 1:size(fwd.usave,1)
    for j in 1:size(fwd.usave,2)
        itp = LinearInterpolation(reverse(fwd.usave[i,j,2:end-1]), reverse(gp.x[1,2:end-1]))
        height[i,j] = itp(0.0)
    end
end

av_height1_Ra500000 = mean(height, dims=2)[:,1] .+ 0.5 .- 0.05

eRa = av_height1_Ra500000.^3 .* 1e4 .* (1.0 .- num.θd)

eRa_Ra500000 = eRa




#Figures
struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)


# fh = Figure(resolution = (1600, 1000))
fh = Figure()
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Dimensionless height")
lines!(fwd_Ra1000.time, av_height1_Ra1000, linewidth = 3, label = "Ra = 10³")
lines!(fwd_Ra10000_long.time, av_height1_Ra10000_long.+0.4499000000000003, linewidth = 3, label = "Ra = 10⁴")
# lines!(fwd_Ra50000_long.time, av_height1_Ra50000_long, linewidth = 3, label = "Ra = 5×10⁴")
lines!(fwd_Ra100000_long.time, av_height1_Ra100000_long, linewidth = 3, label = "Ra = 10⁵")
lines!(fwd_Ra500000.time, av_height1_Ra500000, linewidth = 3, label = "Ra = 5×10⁵")
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()
# Makie.save("/home/tf/All/Thesis/Figures/Chapter5/height_Ra_TM0.png", fu)

# feRa = Figure(resolution = (1600, 1000))
feRa = Figure()
colsize!(feRa.layout, 1, Aspect(1, 1.0))
ax = Axis(feRa[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Effective Rayleigh number",
    yscale = log10, yticks = LogTicks(IntegerTicks()))
# ax.yticks = [1e-6]
lines!(fwd_Ra1000.time[4:end], eRa_Ra1000[4:end]./4, linewidth = 3, label = "Ra = 10³")
lines!(fwd_Ra10000_long.time[4:end], eRa_Ra10000_long[4:end]./4, linewidth = 3, label = "Ra = 10⁴")
# lines!(fwd_Ra50000_long.time[4:end], eRa_Ra50000_long[4:end]./4, linewidth = 3, label = "Ra = 5×10⁴")
lines!(fwd_Ra100000_long.time[4:end], eRa_Ra100000_long[4:end]./4, linewidth = 3, label = "Ra = 10⁵")
lines!(fwd_Ra500000.time[4:end], eRa_Ra500000[4:end]./4, linewidth = 3, label = "Ra = 5×10⁵")
hlines!(ax, [1707.76], color = :black, linestyle =:dash, linewidth = 3, label = "Raₑ = 1707.76")
resize_to_layout!(feRa)
axislegend(position = :rb)
feRa = current_figure()
# Makie.save("/home/tf/All/Thesis/Figures/Chapter5/eRa_TM0.png", feRa)
















#Figures iterations
fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, Ratio))
for i = 1:size(fwd_Ra10000_long.TLsave,1)
    ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
    hidedecorations!(ax)
    heatmap!(gp_Ra10000_long.y[:,1], gp_Ra10000_long.x[1,:], fwd_Ra10000_long.Tsave[i,:,:], colormap = :thermal)
    contour!(gp_Ra10000_long.y[:,1], gp_Ra10000_long.x[1,:], fwd_Ra10000_long.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
    resize_to_layout!(fp)
    Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", 10000))-i$(i).png", fp)
end


fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, Ratio))
for i = 1:size(fwd_Ra100000_long.TLsave,1)
    ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
    hidedecorations!(ax)
    heatmap!(gp_Ra100000_long.y[:,1], gp_Ra100000_long.x[1,:], fwd_Ra100000_long.Tsave[i,:,:], colormap = :thermal)
    contour!(gp_Ra100000_long.y[:,1], gp_Ra100000_long.x[1,:], fwd_Ra100000_long.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
    resize_to_layout!(fp)
    Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", 100000))-i$(i).png", fp)
end


fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, Ratio))
for i = 1:size(fwd_Ra50000_long.TLsave,1)
    ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
    hidedecorations!(ax)
    heatmap!(gp_Ra50000_long.y[:,1], gp_Ra50000_long.x[1,:], fwd_Ra50000_long.Tsave[i,:,:], colormap = :thermal)
    contour!(gp_Ra50000_long.y[:,1], gp_Ra50000_long.x[1,:], fwd_Ra50000_long.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
    resize_to_layout!(fp)
    Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", 50000))-i$(i).png", fp)
end