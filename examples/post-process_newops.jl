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

#-------------------------------------------------#
#POST PROCESS DATA
#-------------------------------------------------#

X = (hRa5, hRa25)
Y = (fwd_Ra5, fwd_Ra25)
Z = (gp_Ra5, gp_Ra25)
NX = (nx_Ra5, nx_Ra25)
NY = (ny_Ra5, ny_Ra25)

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
        x[t] = sum / nx
    end
end


eRa5 = Ra5*(hRa5 .+ 0.5).^3
eRa25 = Ra25*(hRa25 .+ 0.5).^3

tmp_Ra5 = findall(x -> x <= 1707.76, eRa5)
tmp_Ra25 = findall(x -> x <= 1707.76, eRa25)

iRae_Ra5 = length(tmp_Ra5)
iRae_Ra25 = length(tmp_Ra25)

# tRae_Ra3 = 0.5*(fwd_Ra3.time[iRae_Ra3] + fwd_Ra3.time[iRae_Ra3+1])
# tRae_Ra4 = 0.5*(fwd_Ra4.time[iRae_Ra4] + fwd_Ra4.time[iRae_Ra4+1])
# tRae_Ra5 = 0.5*(fwd_Ra5.time[iRae_Ra5] + fwd_Ra5.time[iRae_Ra5+1])
# tRae_Ra6 = 0.5*(fwd_Ra6.time[iRae_Ra6] + fwd_Ra6.time[iRae_Ra6+1])

# iRae_Ra3 = length(tmp_Ra3)
# iRae_Ra4 = length(tmp_Ra4)
# val_Rae_Ra3 = 0.5*(eRa3[iRae_Ra3] + eRa3[iRae_Ra3+1])
# val_Rae_Ra4 = 0.5*(eRa4[iRae_Ra4] + eRa4[iRae_Ra4+1])
# val_Rae_Ra5 = 0.5*(eRa5[iRae_Ra5] + eRa5[iRae_Ra5+1])
# val_Rae_Ra6 = 0.5*(eRa6[iRae_Ra6] + eRa6[iRae_Ra6+1])

# avh_Rae_Ra3 = 0.5*(avhRa3_f[iRae_Ra3] + avhRa3_f[iRae_Ra3+1])
# avh_Rae_Ra4 = 0.5*(avhRa4_f[iRae_Ra4] + avhRa4_f[iRae_Ra4+1])
# avh_Rae_Ra5 = 0.5*(avhRa5_f[iRae_Ra5] + avhRa5_f[iRae_Ra5+1])
# avh_Rae_Ra6 = 0.5*(avhRa6_f[iRae_Ra6] + avhRa6_f[iRae_Ra6+1])

#-------------------------------------------------#
#PLOT DATA
#-------------------------------------------------#

# struct IntegerTicks end

# Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fh = Figure(resolution = (1600, 1600))
# fontsize_theme = Theme(fontsize = 50)
# set_theme!(fontsize_theme)
# colsize!(fh.layout, 1, Aspect(1, 1.0))
# ax = Axis(fh[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Normalized height")
# lines!(fwd_Ra3.time[1:iRa3], avhRa3_f, linewidth = 7, label = "Ra = 10³")
# lines!(fwd_Ra4.time[1:iRa4], avhRa4_f, linewidth = 7, label = "Ra = 10⁴")
# lines!(fwd_Ra5.time[1:iRa5], avhRa5_f, linewidth = 7, label = "Ra = 10⁵")
# lines!(fwd_Ra6.time[1:iRa6], avhRa6_f, linewidth = 7, label = "Ra = 10⁶")
# scatter!(tRae_Ra5, avh_Rae_Ra5, color =:black, markersize = 25)
# scatter!(tRae_Ra6, avh_Rae_Ra6, color =:black, markersize = 25, label = "Raₑ = 1707.76")
# resize_to_layout!(fh)
# axislegend(position = :rb)
# fh = current_figure()

# # Makie.save("./figures/avheight.png", fh)


# feRa = Figure(resolution = (1600, 1600))
# fontsize_theme = Theme(fontsize = 50)
# set_theme!(fontsize_theme)
# colsize!(feRa.layout, 1, Aspect(1, 1.0))
# ax = Axis(feRa[1,1], aspect = 1, xlabel = "Dimensionless time", ylabel = "Raₑ",
#     yscale = log10, yticks = LogTicks(IntegerTicks()))
# lines!(fwd_Ra3.time[2:iRa3], eRa3[2:end], linewidth = 7, label = "Ra = 10³")
# lines!(fwd_Ra4.time[2:iRa4], eRa4[2:end], linewidth = 7, label = "Ra = 10⁴")
# lines!(fwd_Ra5.time[2:iRa5], eRa5[2:end], linewidth = 7, label = "Ra = 10⁵")
# lines!(fwd_Ra6.time[4:iRa6], eRa6[4:end], linewidth = 7, label = "Ra = 10⁶")
# hlines!(ax, [1707.76], color = :black, linestyle =:dash, linewidth = 7, label = "Raₑ = 1707.76")
# resize_to_layout!(feRa)
# axislegend(position = :rb)
# feRa = current_figure()

# Makie.save("./figures/eRa.png", feRa)

# # fp = Figure(resolution = (1600, 1000))
# # colsize!(fp.layout, 1, Aspect(1, Ratio))
# # for i = 1:iRa6
# #     ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# #     hidedecorations!(ax)
# #     heatmap!(gp_Ra6.y[:,1], gp_Ra6.x[1,:], fwd_Ra6.Tsave[i,:,:], colormap = :thermal)
# #     contour!(gp_Ra6.y[:,1], gp_Ra6.x[1,:], fwd_Ra6.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
# #     resize_to_layout!(fp)
# #     Makie.save("/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/Figures_iterations/Ra_$(@sprintf("%.1e", Ra6))-i$(i).png", fp)
# # end


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
# u = @lift(fwd_Ra6.usave[$time,:,:])
# T = @lift(fwd_Ra6.Tsave[$time,:,:])
# # fp = Figure(resolution = (1600, 200))
# fp = Figure()
# colsize!(fp.layout, 1, Aspect(1, Ratio))
# ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# hidedecorations!(ax)
# heatmap!(fp[1,1],gp_Ra6.y[:,1], gp_Ra6.x[1,:], T, colormap = :thermal)
# contour!(fp[1,1],gp_Ra6.y[:,1], gp_Ra6.x[1,:], u, levels = 0:0, color=:black, linewidth = 2);
# resize_to_layout!(fp)
# trim!(fp.layout)
# # contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
# framerate = 15
# timestamps = range(2, iRa6, step=1)
# record(fp, "test.mp4", timestamps;
#         framerate = framerate) do t
#     time[] = t
# end