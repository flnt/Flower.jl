using JLD2
using Flower
using Interpolations

Ratio = 8
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

#-------------------------------------------------#
#EXTRACT DATA
#-------------------------------------------------#

#RAYLEIGH 10^3
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/rayleigh-data/St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+03.jld2"

num_St0_1_Ra3 = num
gp_St0_1_Ra3 = gp
gu_St0_1_Ra3 = gu
gv_St0_1_Ra3 = gv
phS_St0_1_Ra3 = phS
phL_St0_1_Ra3 = phL
fwd_St0_1_Ra3 = fwd
Ra3 = 10^3

#RAYLEIGH 10^4
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/rayleigh-data/St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+04.jld2"

num_St0_1_Ra4 = num
gp_St0_1_Ra4 = gp
gu_St0_1_Ra4 = gu
gv_St0_1_Ra4 = gv
phS_St0_1_Ra4 = phS
phL_St0_1_Ra4 = phL
fwd_St0_1_Ra4 = fwd
Ra4 = 10^4

#RAYLEIGH 10^5
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/rayleigh-data/St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+05.jld2"

num_St0_1_Ra5 = num
gp_St0_1_Ra5 = gp
gu_St0_1_Ra5 = gu
gv_St0_1_Ra5 = gv
phS_St0_1_Ra5 = phS
phL_St0_1_Ra5 = phL
fwd_St0_1_Ra5 = fwd
Ra5 = 10^5


#RAYLEIGH 10^6
JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/rayleigh-data/St_0.4-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_2.5-Ra_1.0e+06.jld2"

num_St0_4_Ra6 = num
gp_St0_4_Ra6 = gp
gu_St0_4_Ra6 = gu
gv_St0_4_Ra6 = gv
phS_St0_4_Ra6 = phS
phL_St0_4_Ra6 = phL
fwd_St0_4_Ra6 = fwd
Ra6 = 10^6



hRa3 = zeros(size(fwd_St0_1_Ra3.usave)[1:2])
hRa4 = zeros(size(fwd_St0_1_Ra4.usave)[1:2])
hRa5 = zeros(size(fwd_St0_1_Ra5.usave)[1:2])
hRa6 = zeros(size(fwd_St0_4_Ra6.usave)[1:2])

X = (hRa3, hRa4, hRa5, hRa6)
Y = (fwd_St0_1_Ra3, fwd_St0_1_Ra4, fwd_St0_1_Ra5, fwd_St0_4_Ra6)
Z = (gp_St0_1_Ra3, gp_St0_1_Ra4, gp_St0_1_Ra5, gp_St0_4_Ra6)

for (x,y,z) in zip(X,Y,Z)
    for i in 1:size(y.usave,1)
        for j in 1:size(y.usave,2)
            if sum(y.usave[i,j,2:end-1]) != 0.0
                itp = LinearInterpolation(reverse(y.usave[i,j,2:end-1]), reverse(z.x[1,2:end-1]))
                x[i,j] = itp(0.0)
            end
        end
    end
end

avhRa3 = mean(hRa3, dims=2)[:,1] # .+ 0.5 .- 0.05
avhRa4 = mean(hRa4, dims=2)[:,1] # .+ 0.5 .- 0.05
avhRa5 = mean(hRa5, dims=2)[:,1] # .+ 0.5 .- 0.05
avhRa6 = mean(hRa6, dims=2)[:,1] # .+ 0.5 .- 0.05

iRa3 = findmax(avhRa3)[2]
iRa4 = findmax(avhRa4)[2]
iRa5 = findmax(avhRa5)[2]
iRa6 = findmax(avhRa6)[2]

avhRa3_f = avhRa3[1:iRa3] .+ 0.5 .- 0.05
avhRa4_f = avhRa4[1:iRa4] .+ 0.5 .- 0.05
avhRa5_f = avhRa5[1:iRa5] .+ 0.5 .- 0.05
avhRa6_f = avhRa6[1:iRa6] .+ 0.5 .- 0.05


eRa3 = (1/Ratio)*avhRa3_f.^3 .* 1e3 .* (1.0 .- num_St0_1_Ra3.θd)
eRa4 = (1/Ratio)*avhRa4_f.^3 .* 1e4 .* (1.0 .- num_St0_1_Ra4.θd)
eRa5 = (1/Ratio)*avhRa5_f.^3 .* 1e5 .* (1.0 .- num_St0_1_Ra5.θd)
eRa6 = (1/Ratio)*avhRa6_f.^3 .* 1e6 .* (1.0 .- num_St0_4_Ra6.θd)

tmp_Ra3 = findall(x -> x <= 1707.76, eRa3)
tmp_Ra4 = findall(x -> x <= 1707.76, eRa4)
tmp_Ra5 = findall(x -> x <= 1707.76, eRa5)
tmp_Ra6 = findall(x -> x <= 1707.76, eRa6)

iRae_Ra3 = length(tmp_Ra3)
iRae_Ra4 = length(tmp_Ra4)
iRae_Ra5 = length(tmp_Ra5)
iRae_Ra6 = length(tmp_Ra6)

tRae_Ra3 = fwd_St0_1_Ra3.time[iRae_Ra3]
tRae_Ra4 = fwd_St0_1_Ra4.time[iRae_Ra4]
tRae_Ra5 = fwd_St0_1_Ra5.time[iRae_Ra5]
tRae_Ra6 = fwd_St0_4_Ra6.time[iRae_Ra6]

val_Rae_Ra3 = eRa3[iRae_Ra3]
val_Rae_Ra4 = eRa4[iRae_Ra4]
val_Rae_Ra5 = eRa5[iRae_Ra5]
val_Rae_Ra6 = eRa6[iRae_Ra6]

avh_Rae_Ra3 = avhRa3_f[iRae_Ra3]
avh_Rae_Ra4 = avhRa4_f[iRae_Ra4]
avh_Rae_Ra5 = avhRa5_f[iRae_Ra5]
avh_Rae_Ra6 = avhRa6_f[iRae_Ra6]



struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fh = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Normalized height")
lines!(fwd_St0_1_Ra3.time[1:iRa3], avhRa3_f, linewidth = 7, label = "Ra = 10³")
lines!(fwd_St0_1_Ra4.time[1:iRa4], avhRa4_f, linewidth = 7, label = "Ra = 10⁴")
lines!(fwd_St0_1_Ra5.time[1:iRa5], avhRa5_f, linewidth = 7, label = "Ra = 10⁵")
lines!(fwd_St0_4_Ra6.time[1:iRa6], avhRa6_f, linewidth = 7, label = "Ra = 10⁶")
scatter!(tRae_Ra5, avh_Rae_Ra5, color =:black, markersize = 25)
scatter!(tRae_Ra6, avh_Rae_Ra6, color =:black, markersize = 25, label = "Raₑ = 1707.76")
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()


feRa = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(feRa.layout, 1, Aspect(1, 1.0))
ax = Axis(feRa[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Raₑ",
    yscale = log10, yticks = LogTicks(IntegerTicks()))
lines!(fwd_St0_1_Ra3.time[2:iRa3], eRa3[2:end], linewidth = 7, label = "Ra = 10³")
lines!(fwd_St0_1_Ra4.time[2:iRa4], eRa4[2:end], linewidth = 7, label = "Ra = 10⁴")
lines!(fwd_St0_1_Ra5.time[2:iRa5], eRa5[2:end], linewidth = 7, label = "Ra = 10⁵")
lines!(fwd_St0_4_Ra6.time[4:iRa6], eRa6[4:end], linewidth = 7, label = "Ra = 10⁶")
hlines!(ax, [1707.76], color = :black, linestyle =:dash, linewidth = 7, label = "Raₑ = 1707.76")
resize_to_layout!(feRa)
axislegend(position = :rb)
feRa = current_figure()


time = Makie.Observable(1)
u = @lift(fwd_St0_1_Ra4.usave[$time,:,:])
T = @lift(fwd_St0_1_Ra4.Tsave[$time,:,:])
fp = Figure(resolution = (2560, 320))
# fp = Figure()
colsize!(fp.layout, 1, Aspect(1, Ratio))
ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
hidedecorations!(ax)
heatmap!(fp[1,1],gp_St0_1_Ra4.y[:,1], gp_St0_1_Ra4.x[1,:], T, colormap = :thermal)
contour!(fp[1,1],gp_St0_1_Ra4.y[:,1], gp_St0_1_Ra4.x[1,:], u, levels = 0:0, color=:black, linewidth = 2);
resize_to_layout!(fp)
trim!(fp.layout)
# contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
framerate = 48
timestamps = range(2, iRa4, step=1)
record(fp, "St0_1_Ra4.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end

# St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+03
# St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+04
# St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+05
# St_0.1-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_10.0-Ra_1.0e+06

# St_0.2-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_5.0-Ra_1.0e+03
# St_0.2-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_5.0-Ra_1.0e+04
# St_0.2-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_5.0-Ra_1.0e+05
# St_0.2-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_5.0-Ra_1.0e+06

# St_0.4-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_2.5-Ra_1.0e+03
# St_0.4-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_2.5-Ra_1.0e+04
# St_0.4-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_2.5-Ra_1.0e+05
# St_0.4-nx_80-ny_640-ratio_8-maxiter_1.0e+05-TM_0.0-T1_1.0-T2_0.0-λ_2.5-Ra_1.0e+06

xs = LinRange(0, 10, 100)
ys = LinRange(0, 15, 100)
zs = [cos(x) * sin(y) for x in xs, y in ys]

contour(gp_St0_1_Ra5.y[:,1], gp_St0_1_Ra5.x[1,:], fwd_St0_1_Ra5.usave[700,:,:], levels = 0:0, color=:black, linewidth = 2)

using GLMakie