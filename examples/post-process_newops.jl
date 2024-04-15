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

# #RAYLEIGH 10^6
# JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_128_ny_32_ratio_4_maxiter_3.0e+03_TM_0.0_T1_0.7_T2_-0.3_St_0.1_Ra_1.0e+06.jld2"

# gp_Ra5, gu_Ra5, gv_Ra5 = init_meshes(num)
# num_Ra5 = num
# fwd_Ra5 = fwd
# Ra5 = 10^6
# tmax_Ra5 = findmax(fwd_Ra5.t)[2]
# hRa5 = zeros(tmax_Ra5)
# nx_Ra5 = size(num_Ra5.x,1) - 1
# ny_Ra5 = size(num_Ra5.y,1) - 1

# for t in 1:length(hRa5)
#     h = 0.
#     update_ls_data(num_Ra5, gp_Ra5, gu_Ra5, gv_Ra5, 1, fwd_Ra5.u[1,t,:,:], fwd_Ra5.κ[1,t,:,:], [Stefan()], Stefan(), true, false, true, false)
#     for jj in gp_Ra5.LS[1].MIXED
#         h += gp_Ra5.y[jj] + num_Ra5.Δ * gp_Ra5.LS[1].mid_point[jj].y 
#     end
#     hRa5[t] =   h / length(gp_Ra5.LS[1].MIXED) + 1 / 2
# end

#RAYLEIGH 10^5
JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_512_ny_64_ratio_8_maxiter_1.0e+04_TM_0.0_T1_0.7_T2_-0.3_St_10.0_Ra_1.0e+05.jld2"

gp_Ra5, gu_Ra5, gv_Ra5 = init_meshes(num)
num_Ra5 = num
fwd_Ra5 = fwd
Ra5 = 10^5
tmax_Ra5 = findmax(fwd_Ra5.t)[2]
hRa5 = zeros(tmax_Ra5)
nx_Ra5 = size(num_Ra5.x,1) - 1
ny_Ra5 = size(num_Ra5.y,1) - 1

# for t in 1:length(hRa5)
#     h = 0.
#     update_ls_data(num_Ra5, gp_Ra5, gu_Ra5, gv_Ra5, 1, fwd_Ra5.u[1,t,:,:], fwd_Ra5.κ[1,t,:,:], [Stefan()], Stefan(), true, false, true, false)
#     idx = intersect(gp_Ra5.LS[1].MIXED, gp_Ra5.ind.all_indices[:,1:512])
#     for jj in idx #gp_Ra5.LS[1].MIXED
#         # h += gp_Ra5.y[jj] + num_Ra5.Δ * gp_Ra5.LS[1].mid_point[jj].y
#         tmp = gp_Ra5.y[jj] + num_Ra5.Δ * gp_Ra5.LS[1].mid_point[jj].y + 0.5 
#         if h < tmp
#             h = tmp
#         end
#     end
#     hRa5[t] =   h 
#     # hRa5[t] =   h / length(idx) + 1 / 2 #h / length(gp_Ra5.LS[1].MIXED) + 1 / 2
# end

#RAYLEIGH 2.10^5
JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_512_ny_64_ratio_8_maxiter_1.0e+04_TM_0.0_T1_0.7_T2_-0.3_St_10.0_Ra_2.0e+05.jld2"

gp_Ra25, gu_Ra25, gv_Ra25 = init_meshes(num)
num_Ra25 = num
fwd_Ra25 = fwd
Ra25 = 2*10^5
tmax_Ra25 = findmax(fwd_Ra25.t)[2]
hRa25 = zeros(tmax_Ra25)
nx_Ra25 = size(num_Ra25.x,1) - 1
ny_Ra25 = size(num_Ra25.y,1) - 1

# for t in 1:length(hRa25)
#     h = 0.
#     update_ls_data(num_Ra25, gp_Ra25, gu_Ra25, gv_Ra25, 1, fwd_Ra25.u[1,t,:,:], fwd_Ra25.κ[1,t,:,:], [Stefan()], Stefan(), true, false, true, false)
#     idx = intersect(gp_Ra25.LS[1].MIXED, gp_Ra25.ind.all_indices[:,1:512])
#     for jj in idx #gp_Ra25.LS[1].MIXED
#         # h += gp_Ra25.y[jj] + num_Ra25.Δ * gp_Ra25.LS[1].mid_point[jj].y 
#         tmp = gp_Ra25.y[jj] + num_Ra25.Δ * gp_Ra25.LS[1].mid_point[jj].y + 0.5 
#         if h < tmp
#             h = tmp
#         end
#     end
#     hRa25[t] =   h 
#     # hRa25[t] =   h / length(idx) + 1 / 2 #h / length(gp_Ra25.LS[1].MIXED) + 1 / 2
# end

#RAYLEIGH 4.10^5
JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_512_ny_64_ratio_8_maxiter_1.0e+04_TM_0.0_T1_0.7_T2_-0.3_St_10.0_Ra_4.0e+05.jld2"

gp_Ra45, gu_Ra45, gv_Ra45 = init_meshes(num)
num_Ra45 = num
fwd_Ra45 = fwd
Ra45 = 4*10^5
tmax_Ra45 = findmax(fwd_Ra45.t)[2]
hRa45 = zeros(tmax_Ra45)
nx_Ra45 = size(num_Ra45.x,1) - 1
ny_Ra45 = size(num_Ra45.y,1) - 1

# for t in 1:length(hRa45)
#     h = 0.
#     update_ls_data(num_Ra45, gp_Ra45, gu_Ra45, gv_Ra45, 1, fwd_Ra45.u[1,t,:,:], fwd_Ra45.κ[1,t,:,:], [Stefan()], Stefan(), true, false, true, false)
#     idx = intersect(gp_Ra45.LS[1].MIXED, gp_Ra45.ind.all_indices[:,1:512])
#     for jj in idx #gp_Ra45.LS[1].MIXED
#         # h += gp_Ra45.y[jj] + num_Ra45.Δ * gp_Ra45.LS[1].mid_point[jj].y 
#         tmp = gp_Ra45.y[jj] + num_Ra45.Δ * gp_Ra45.LS[1].mid_point[jj].y + 0.5 
#         if h < tmp
#             h = tmp
#         end
#     end
#     hRa45[t] =   h 
#     # hRa45[t] =   h / length(idx) + 1 / 2 #h / length(gp_Ra45.LS[1].MIXED) + 1 / 2
# end

#RAYLEIGH 5.10^4
JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_512_ny_64_ratio_8_maxiter_1.0e+04_TM_0.0_T1_0.7_T2_-0.3_St_10.0_Ra_5.0e+04.jld2"

gp_Ra54, gu_Ra54, gv_Ra54 = init_meshes(num)
num_Ra54 = num
fwd_Ra54 = fwd
Ra54 = 5*10^4
tmax_Ra54 = findmax(fwd_Ra54.t)[2]
hRa54 = zeros(tmax_Ra54)
nx_Ra54 = size(num_Ra54.x,1) - 1
ny_Ra54 = size(num_Ra54.y,1) - 1

#RAYLEIGH 5.10^3
JLD2.@load "/home/tf/Documents/Flower_figures/newops_nx_512_ny_64_ratio_8_maxiter_1.0e+04_TM_0.0_T1_0.7_T2_-0.3_St_10.0_Ra_5.0e+03.jld2"

gp_Ra53, gu_Ra53, gv_Ra53 = init_meshes(num)
num_Ra53 = num
fwd_Ra53 = fwd
Ra53 = 4*10^5
tmax_Ra53 = findmax(fwd_Ra53.t)[2]
hRa53 = zeros(tmax_Ra53)
nx_Ra53 = size(num_Ra53.x,1) - 1
ny_Ra53 = size(num_Ra53.y,1) - 1

#-------------------------------------------------#
#POST PROCESS DATA
#-------------------------------------------------#

X = (hRa5, hRa5, hRa25, hRa45, hRa54, hRa53)
Y = (fwd_Ra5, fwd_Ra5, fwd_Ra25, fwd_Ra45, fwd_Ra54, fwd_Ra53)
Z = (gp_Ra5, gp_Ra5, gp_Ra25, gp_Ra45, gp_Ra54, gp_Ra53)
NX = (nx_Ra5, nx_Ra5, nx_Ra25, nx_Ra45, nx_Ra54, nx_Ra53)
NY = (ny_Ra5, ny_Ra5, ny_Ra25, ny_Ra45,  ny_Ra54, ny_Ra53)

for (x,y,z,nx,ny) in zip(X,Y,Z,NX,NY)
    for t in 1:length(x)
        sum = 0
        for i in 1:nx
            for j in 1:(ny-1)
                diff = y.u[1,t,j,i]*y.u[1,t,j+1,i]
                if diff < 0
                    pos = z.y[j+1,1] + y.u[1,t,j+1,i]
                    sum += pos 
                end
            end
        end
        x[t] = (sum / nx + 0.5)
    end
end



# X = (hRa5, hRa25, hRa45)
# Y = (fwd_Ra5, fwd_Ra25, fwd_Ra45)
# Z = (gp_Ra5, gp_Ra25, gp_Ra45)
# NX = (nx_Ra5, nx_Ra25, nx_Ra45)
# NY = (ny_Ra5, ny_Ra25, ny_Ra45)

# for (x,y,z) in zip(X,Y,Z)
#     for t in 1:length(x)
#         for j in 1:size(y.usave,2)
#             if sum(y.usave[t,j,2:end-1]) != 0.0
#                 itp = LinearInterpolation(reverse(y.usave[i,j,2:end-1]), reverse(z.x[1,2:end-1]))
#                 x[i,j] = itp(0.0)
#             end
#         end
#     end
# end
a = 0.3

eRa5 = Ra5*a*(hRa5).^3
eRa25 = Ra25*a*(hRa25).^3
eRa45 = Ra45*a*(hRa45).^3
eRa54 = Ra54*a*(hRa25).^3
eRa53 = Ra53*a*(hRa45).^3

tmp_Ra5 = findall(x -> x <= 1707.76, eRa5)
tmp_Ra25 = findall(x -> x <= 1707.76, eRa25)
tmp_Ra45 = findall(x -> x <= 1707.76, eRa45)
tmp_Ra54 = findall(x -> x <= 1707.76, eRa54)
tmp_Ra53 = findall(x -> x <= 1707.76, eRa53)

iRae_Ra5 = length(tmp_Ra5)
iRae_Ra25 = length(tmp_Ra25)
iRae_Ra45 = length(tmp_Ra45)
iRae_Ra54 = length(tmp_Ra54)
iRae_Ra53 = length(tmp_Ra53)

aeRa5 = (eRa5[tmp_Ra5[end]+1] - eRa5[tmp_Ra5[end]])/(fwd_Ra5.t[tmp_Ra5[end]+1] - fwd_Ra5.t[tmp_Ra5[end]])
teRa5 = (1707.76 - eRa5[tmp_Ra5[end]])/aeRa5 + fwd_Ra5.t[tmp_Ra5[end]]
ahRa5 = (hRa5[tmp_Ra5[end]+1] - hRa5[tmp_Ra5[end]])/(fwd_Ra5.t[tmp_Ra5[end]+1] - fwd_Ra5.t[tmp_Ra5[end]])
p_hRa5 = ahRa5 * teRa5 + hRa5[tmp_Ra5[end]] - ahRa5 * fwd_Ra5.t[tmp_Ra5[end]]

aeRa25 = (eRa25[tmp_Ra25[end]+1] - eRa25[tmp_Ra25[end]])/(fwd_Ra25.t[tmp_Ra25[end]+1] - fwd_Ra25.t[tmp_Ra25[end]])
teRa25 = (1707.76 - eRa25[tmp_Ra25[end]])/aeRa25 + fwd_Ra25.t[tmp_Ra25[end]]
ahRa25 = (hRa25[tmp_Ra25[end]+1] - hRa25[tmp_Ra25[end]])/(fwd_Ra25.t[tmp_Ra25[end]+1] - fwd_Ra25.t[tmp_Ra25[end]])
p_hRa25 = ahRa25 * teRa25 + hRa25[tmp_Ra25[end]] - ahRa25 * fwd_Ra25.t[tmp_Ra25[end]]

aeRa45 = (eRa45[tmp_Ra45[end]+1] - eRa45[tmp_Ra45[end]])/(fwd_Ra45.t[tmp_Ra45[end]+1] - fwd_Ra45.t[tmp_Ra45[end]])
teRa45 = (1707.76 - eRa45[tmp_Ra45[end]])/aeRa45 + fwd_Ra45.t[tmp_Ra45[end]]
ahRa45 = (hRa45[tmp_Ra45[end]+1] - hRa45[tmp_Ra45[end]])/(fwd_Ra45.t[tmp_Ra45[end]+1] - fwd_Ra45.t[tmp_Ra45[end]])
p_hRa45 = ahRa45 * teRa45 + hRa45[tmp_Ra45[end]] - ahRa45 * fwd_Ra45.t[tmp_Ra45[end]]

aeRa54 = (eRa54[tmp_Ra54[end]+1] - eRa54[tmp_Ra54[end]])/(fwd_Ra54.t[tmp_Ra54[end]+1] - fwd_Ra54.t[tmp_Ra54[end]])
teRa54 = (1707.76 - eRa54[tmp_Ra54[end]])/aeRa54 + fwd_Ra54.t[tmp_Ra54[end]]
ahRa54 = (hRa54[tmp_Ra54[end]+1] - hRa54[tmp_Ra54[end]])/(fwd_Ra54.t[tmp_Ra54[end]+1] - fwd_Ra54.t[tmp_Ra54[end]])
p_hRa54 = ahRa54 * teRa54 + hRa54[tmp_Ra54[end]] - ahRa54 * fwd_Ra54.t[tmp_Ra45[end]]

# aeRa53 = (eRa53[tmp_Ra53[end]+1] - eRa53[tmp_Ra53[end]])/(fwd_Ra53.t[tmp_Ra53[end]+1] - fwd_Ra53.t[tmp_Ra53[end]])
# teRa53 = (1707.76 - eRa53[tmp_Ra53[end]])/aeRa53 + fwd_Ra53.t[tmp_Ra53[end]]
# ahRa53 = (hRa53[tmp_Ra53[end]+1] - hRa53[tmp_Ra53[end]])/(fwd_Ra53.t[tmp_Ra53[end]+1] - fwd_Ra53.t[tmp_Ra53[end]])
# p_hRa53 = ahRa53 * teRa53 + hRa53[tmp_Ra53[end]] - ahRa53 * fwd_Ra53.t[tmp_Ra53[end]]


# tRae_Ra3 = 0.5*(fwd_Ra3.time[iRae_Ra3] + fwd_Ra3.time[iRae_Ra3+1])
# tRae_Ra4 = 0.5*(fwd_Ra4.time[iRae_Ra4] + fwd_Ra4.time[iRae_Ra4+1])
# tRae_Ra5 = 0.5*(fwd_Ra5.time[iRae_Ra5] + fwd_Ra5.time[iRae_Ra5+1])
# tRae_Ra5 = 0.5*(fwd_Ra5.time[iRae_Ra5] + fwd_Ra5.time[iRae_Ra5+1])

# iRae_Ra3 = length(tmp_Ra3)
# iRae_Ra4 = length(tmp_Ra4)
# val_Rae_Ra3 = 0.5*(eRa3[iRae_Ra3] + eRa3[iRae_Ra3+1])
# val_Rae_Ra4 = 0.5*(eRa4[iRae_Ra4] + eRa4[iRae_Ra4+1])
# val_Rae_Ra5 = 0.5*(eRa5[iRae_Ra5] + eRa5[iRae_Ra5+1])
# val_Rae_Ra5 = 0.5*(eRa5[iRae_Ra5] + eRa5[iRae_Ra5+1])

# avh_Rae_Ra3 = 0.5*(avhRa3_f[iRae_Ra3] + avhRa3_f[iRae_Ra3+1])
# avh_Rae_Ra4 = 0.5*(avhRa4_f[iRae_Ra4] + avhRa4_f[iRae_Ra4+1])
# avh_Rae_Ra5 = 0.5*(avhRa5_f[iRae_Ra5] + avhRa5_f[iRae_Ra5+1])
# avh_Rae_Ra5 = 0.5*(avhRa5_f[iRae_Ra5] + avhRa5_f[iRae_Ra5+1])

#-------------------------------------------------#
#PLOT DATA
#-------------------------------------------------#

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fh = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Time", ylabel = "Average height")
ax.yticks = (0.1:0.2:0.9)
ax.xticks = (0:0.1:1.1)
# xlims!(ax, 0, 1)
ylims!(ax, 0, 1)
lines!(fwd_Ra53.t, hRa53, linewidth = 9, label = "Ra = 5 × 10³")
lines!(fwd_Ra54.t, hRa54, linewidth = 9, label = "Ra = 5 × 10⁴")
lines!(fwd_Ra5.t, hRa5, linewidth = 9, label = "Ra = 1 × 10⁵")
lines!(fwd_Ra25.t, hRa25, linewidth = 9, label = "Ra = 2 × 10⁵")
lines!(fwd_Ra45.t, hRa45, linewidth = 9, label = "Ra = 4 × 10⁵")
# lines!(fwd_Ra4.time[1:iRa4], avhRa4_f, linewidth = 7, label = "Ra = 10⁴")
# lines!(fwd_Ra5.time[1:iRa5], avhRa5_f, linewidth = 7, label = "Ra = 10⁵")
# lines!(fwd_Ra5.time[1:iRa5], avhRa5_f, linewidth = 7, label = "Ra = 10⁶")
# scatter!(St*teRa5, p_hRa5, color =:black, markersize = 15, label = "Raₑ = 1707.76")
# scatter!(St*teRa25, p_hRa25, color =:black, markersize = 15)
# scatter!(St*teRa45, p_hRa45, color =:black, markersize = 15)
# scatter!(tRae_Ra5, avh_Rae_Ra5, color =:black, markersize = 25, label = "Raₑ = 1707.76")
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()

Makie.save("./figures/avheight.png", fh)


feRa = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(feRa.layout, 1, Aspect(1, 1.0))
ax = Axis(feRa[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Raₑ",
    yscale = log10, yticks = LogTicks(IntegerTicks()))
lines!(St*fwd_Ra5.t, eRa5 .+ 0.001, linewidth = 7, label = "Ra = 10⁵")
lines!(St*fwd_Ra25.t, eRa25 .+ 0.001, linewidth = 7, label = "Ra = 2 × 10⁵")
lines!(St*fwd_Ra45.t, eRa45 .+ 0.001, linewidth = 7, label = "Ra = 4 × 10⁵")
lines!(St*fwd_Ra54.t, eRa54 .+ 0.001, linewidth = 7, label = "Ra = 5 × 10⁴")
lines!(St*fwd_Ra53.t, eRa53 .+ 0.001, linewidth = 7, label = "Ra = 5 × 10³")
hlines!(ax, [1707.76], color = :black, linestyle =:dash, linewidth = 7, label = "Raₑ = 1707.76")
resize_to_layout!(feRa)
axislegend(position = :rb)
feRa = current_figure()

# Makie.save("./figures/eRa.png", feRa)

fp = Figure(resolution = (1600, 1000))
colsize!(fp.layout, 1, Aspect(1, Ratio))
# for i = 1:length(fwd_Ra45.t)
    ax = Axis(fp[1,1], aspect = Ratio)  # customized as you see fit
    hidedecorations!(ax)
    heatmap!(transpose(fwd_Ra45.u[1,100,:,:]), colormap = :lightrainbow, clim = (-0.3, 0.7))
    # contour!(gp_Ra5.y[:,1], gp_Ra5.x[1,:], fwd_Ra5.usave[i,:,:], levels = 0:0, color=:red, linewidth = 2);
    resize_to_layout!(fp)
    # Makie.save("/home/Documents/Ra_$(@sprintf("%.1e", Ra45))-i$(i).png", fp)
# end


# # fp = Figure(resolution = (1600, 1000))
# # colsize!(fp.layout, 1, Aspect(1, Ratio))
# # for i = 1:iRa5
# #     ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# #     hidedecorations!(ax)
# #     heatmap!(gp_Ra5.y[:,1], gp_Ra5.x[1,:], fwd_Ra5.Tsave[i,:,:], colormap = :thermal)
# #     contour!(gp_Ra5.y[:,1], gp_Ra5.x[1,:], fwd_Ra5.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
# #     resize_to_layout!(fp)
# #     Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", Ra5))_λ2-i$(i).png", fp)
# # end


# #-------------------------------------------------#
# #VIDEOS
# #-------------------------------------------------#



# time = Makie.Observable(1)
# u = @lift(fwd_Ra5.usave[$time,:,:])
# T = @lift(fwd_Ra5.Tsave[$time,:,:])
# # fp = Figure(resolution = (1600, 200))
# fp = Figure()
# colsize!(fp.layout, 1, Aspect(1, Ratio))
# ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# hidedecorations!(ax)
# heatmap!(fp[1,1],gp_Ra5.y[:,1], gp_Ra5.x[1,:], T, colormap = :thermal)
# contour!(fp[1,1],gp_Ra5.y[:,1], gp_Ra5.x[1,:], u, levels = 0:0, color=:black, linewidth = 2);
# resize_to_layout!(fp)
# trim!(fp.layout)
# # contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
# framerate = 15
# timestamps = range(2, iRa5, step=1)
# record(fp, "test.mp4", timestamps;
#         framerate = framerate) do t
#     time[] = t
# end