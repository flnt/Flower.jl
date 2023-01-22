using JLD2
using Flower
using Interpolations

Ratio = 1
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

St = 0.2
TM = 0.0
LH = 5.0
Ra = 1.0e+06

JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/RB_opt_data/opt_St_0.2-nx_64-ny_64-ratio_1-maxiter_3.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_5.0-Ra_1.0e+06.jld2"

# JLD2.@load "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/julia-data/rayleigh-nx_120-ny_960-ratio_8-maxiter_6.0e+04-TM_0.0-T1_1.0-T2_0.0-λ_2.0-Ra_1.0e+05.jld2"

h = zeros(size(fwd.usave)[1:2]);

for i in 1:size(fwd.usave,1)
    for j in 1:size(fwd.usave,2)
        if sum(fwd.usave[i,j,2:end-1]) != 0.0
            itp = LinearInterpolation(reverse(fwd.usave[i,j,2:end-1]), reverse(gp.x[1,2:end-1]))
            h[i,j] = itp(0.0)
        end
    end
end

avh = mean(h, dims=2)[:,1];

II =  length(avh); #findmax(avh)[2];

avh_f = avh[1:II] .+ 0.5 .- 0.05;

eRa = (1/Ratio)*avh_f.^3 .* Ra .* (1.0 .- num.θd);

tmp = findall(x -> x <= 1707.76, eRa);

iRae = length(tmp)

# fp = Figure(resolution = (267, 200))
# colsize!(fp.layout, 1, Aspect(1, Ratio))
# for i = 2:II
#     ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
#     hidedecorations!(ax)
#     heatmap!(gp.y[:,1], gp.x[1,:], fwd.Tsave[i,:,:], colormap = :thermal)
#     contour!(gp.y[:,1], gp.x[1,:], fwd.usave[i,:,:], levels = 0:0, color=:black, linewidth = 2);
#     resize_to_layout!(fp)
#     Makie.save("/home/tf/All/Soutenance/images/TM$(TM)_St$(St)_LH$(LH)_Ra_$(@sprintf("%.1e", Ra))-i$(i).png", fp)
# end

time = Makie.Observable(1)
u = @lift(fwd.usave[$time,:,:])
T = @lift(fwd.Tsave[$time,:,:])
# fp = Figure(resolution = (2560, 320))
fp = Figure(resolution = (400, 300))
colsize!(fp.layout, 1, Aspect(1, Ratio))
ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
hidedecorations!(ax)
hidespines!(ax)
heatmap!(fp[1,1],gp.y[:,1], gp.x[1,:], T, colormap = :thermal)
contour!(fp[1,1],gp.y[:,1], gp.x[1,:], u, levels = 0:0, color=:red, linewidth = 3);
resize_to_layout!(fp)
trim!(fp.layout)
# contour!(num.H, num.H, u, levels = 0:0, color=:black, linewidth = 2);
framerate = 48
timestamps = range(2, II, step=1)
record(fp, "/media/tf/be46b01f-eb0c-4298-a160-d10e4e87b3b9/RB_opt_data/opt_TM$(TM)_St$(St)_LH$(LH)_Ra_$(@sprintf("%.1e", Ra)).mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end

# fp = Figure(resolution = (400, 300))
# colsize!(fp.layout, 1, Aspect(1, Ratio))
# ax = Axis(fp[1,1], aspect = Ratio, xticks = tcks, yticks = tcks)  # customized as you see fit
# hidedecorations!(ax)
# hidespines!(ax)
# resize_to_layout!(fp)

using GLMakie

f = Figure()
ax = Axis3(f[1,1], aspect = (8, 2, 1)) 
# Makie.hidedecorations!(ax)
# ax.xticks = 1000:1001
# ax.yticks = 1000:1001
# ax.zticks = 1000:1001
ax.xlabel = "X"
ax.ylabel = "Time"
ax.zlabel = "Y"
for i = 10:II
contour3d!(f[1,1],gp.y[:,1],gp.x[1,:],  fwd.usave[i,:,:].+avh_f[i],transformation=(:xz, 0),
    transparency = true, levels = avh_f[i]:avh_f[i], linewidth = 10,
    colormap=Reverse(:thermal), fill = true, colorrange = (0, avh_f[II]) );
end
resize_to_layout!(f)
f