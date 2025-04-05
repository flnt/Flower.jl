using JLD2
using Flower
using Interpolations

Ratio = 4
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

#JLD2.@save "/home/tf/Documents/RB_opt/opt_data.jld2" 
# num gp gu gv fwd_des opt_p opt_S opt_L opt_u opt_uu opt_RB res


# JLD2.@load "/home/tf/Documents/RB_opt/RB_nx_128_ny_32_ratio_4_tend_1.0e+00_TM_0.0_T1_0.7_T2_-0.3_St_1.0_Ra_5.0e+03.jld2"

# temp5103 = temp
# ls5103 = ls
# RB5103 = RB

struct IntegerTicks end

Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fh = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Time", ylabel = "Average height")
ax.yticks = (0.1:0.2:1.0)
ax.xticks = (0:0.2:1.1)
xlims!(ax, 0, 1.025)
ylims!(ax, 0, 1.0)
lines!(fwd.RB[1,1:end-1280], fwd.RB[2,1:end-1280], linewidth = 9, label = "Desired solution", color=:red)
lines!(RB8104[1,1:end-1000], RB8104[2,1:end-1000], linewidth = 9, label = "Iteration 0", color=:purple)
lines!(RB4104[1,1:end-2], RB4104[2,1:end-2], linewidth = 9, label = "Ra = 4 × 10⁴", color=:blue)
lines!(RB1104[1,1:end-2], RB1104[2,1:end-2], linewidth = 9, label = "Ra = 1 × 10⁴", color=:green)
# lines!(RB5103[1,:], RB5103[2,:], linewidth = 9, label = "Ra = 5 × 10³")
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()

Makie.save("/home/tf/Documents/RB_opt/RB_opt_avheight.png", fh)

fh = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], yscale = log10, aspect = 1, xlabel = "Time", ylabel = "Effective Rayleigh number")
# ax.yticks = (0.:1000:1*0.25*findmax(RB1105[3,:])[1]+500)
ax.xticks = (0:0.2:1.1)
xlims!(ax, 0, 1.025)
ylims!(ax, 1, 1*0.25*findmax(RB1105[3,:])[1]+500)
lines!(RB1105[1,1:end-1280], 1*RB1105[3,1:end-1280]/4, linewidth = 9, label = "Ra = 1 × 10⁵", color=:red)
lines!(RB8104[1,1:end-1000], 1*RB8104[3,1:end-1000]/4, linewidth = 9, label = "Ra = 8 × 10⁴", color=:purple)
lines!(RB4104[1,:], 1*RB4104[3,:]/4, linewidth = 9, label = "Ra = 4 × 10⁴", color=:blue)
lines!(RB1104[1,:], 1*RB1104[3,:]/4, linewidth = 9, label = "Ra = 1 × 10⁴", color=:green)
# lines!(RB5103[1,:], 0.7*RB5103[3,:]/4, linewidth = 9, label = "Ra = 5 × 10³")
lines!(RB1104[1,:], 1707.7 .+ 0.0*RB1104[3,:]/4, linewidth = 9, label = "Ra = 1707.7", color=:black)
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()

Makie.save("/home/tf/Documents/RB_opt/RB_opt_effectiveRa.png", fh)


fh = Figure(resolution = (1600, 1600))
fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)
colsize!(fh.layout, 1, Aspect(1, 1.0))
ax = Axis(fh[1,1], aspect = 1, xlabel = "Time", ylabel = "Effective Rayleigh number")
# ax.yticks = (0.:1000:1*0.25*findmax(RB1105[3,:])[1]+500)
ax.xticks = (0:0.2:1.1)
xlims!(ax, 0, 1.025)
ylims!(ax, 1, 1*0.25*findmax(RB1105[3,:])[1]+500)
lines!(RB1105[1,1:end-1280], 1*RB1105[3,1:end-1280]/4, linewidth = 9, label = "Ra = 1 × 10⁵", color=:red)
lines!(RB8104[1,1:end-1000], 1*RB8104[3,1:end-1000]/4, linewidth = 9, label = "Ra = 8 × 10⁴", color=:purple)
lines!(RB4104[1,1:end-1], 1*RB4104[3,1:end-1]/4, linewidth = 9, label = "Ra = 4 × 10⁴", color=:blue)
lines!(RB1104[1,1:end-1], 1*RB1104[3,1:end-1]/4, linewidth = 9, label = "Ra = 1 × 10⁴", color=:green)
# lines!(RB5103[1,:], 0.7*RB5103[3,:]/4, linewidth = 9, label = "Ra = 5 × 10³")
lines!(RB1104[1,:], 1707.7 .+ 0.0*RB1104[3,:]/4, linewidth = 9, label = "Ra = 1707.7", color=:black)
resize_to_layout!(fh)
axislegend(position = :rt)
fh = current_figure()

Makie.save("/home/tf/Documents/RB_opt/RB_opt_effectiveRa_linear.png", fh)
for i = 1:11
    # Set adjustable temperature limits (change these as needed)
    temp_min = -0.3    # e.g. minimum temperature
    temp_max = 0.7 # e.g. maximum temperature
    
    F1 = Figure(size = (1600, 400))
    
    # Create the axis without xlims, ylims, ticklabelsize, or labelfontsize keywords.
    ax = Axis(F1[1, 1],
        aspect = ratio,
        xlabel = "x", 
        ylabel = "y",
        xticks = -2:0.5:2,
        yticks = [-0.5, 0, 0.5]
    )
    
    # Set the axis limits after creation.
    xlims!(ax, (-2, 2))
    ylims!(ax, (-0.5, 0.5))
    
    # Increase the font sizes by setting the properties directly.
    ax.xticklabelsize[] = 30         # size of tick labels on the x axis
    ax.yticklabelsize[] = 30         # size of tick labels on the y axis
    ax.xlabelsize[]     = 34         # size of the x-axis label
    ax.ylabelsize[]     = 34         # size of the y-axis label
    
    # Plot the filled contour (using clims to control the color range)
    hm = contourf!(ax, gp.x[1, :], gp.y[:, 1],
        # min.(abs.(abs.(opt_S[1]' + opt_L[1]') - abs.(fwdS_des.T[end,:,:]' + fwdL_des.T[end,:,:]')), 1.99), 
        temp1105[i,:,:]', 
        colormap = :dense,
        levels = 30,
        clims = (-0.31, 0.71)
    )
    # hm.colorrange[] = (temp_min, temp_max)
    hm.clims[] = (temp_min, temp_max)

    # Add a contour line (for example, at level 0) in red.
    contour!(ax, gp.x[1, :], gp.y[:, 1], ls1105[i,:,:]',
        levels = [0],
        color = :red,
        linewidth = 5
    )
    
    
    # Add a colorbar to the figure next to the axis.
    Colorbar(F1[1, 2];
        colormap = :dense,
        label = "Temperature",
        labelsize = 34,
        ticklabelsize = 30,
        width = 40,
        height = 300,
        limits = (-0.31, 0.71)
    )
    
    # resize_to_layout!(F1)
    F1  # display the figure
    Makie.save("/home/tf/Documents/RB_opt/RB_opt_Ra1105_i$(i-1).png", F1)
end