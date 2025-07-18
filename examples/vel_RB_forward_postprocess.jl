using JLD2
using Flower
using Interpolations

Ratio = 4
L0 = 1.
tcks = -Ratio*L0/2:2:Ratio*L0/2
lim = L0 / 2

JLD2.@load "/home/tf/Documents/RB_opt/vel_RB_nx_128_ny_32_ratio_4_tend_1.0e+00_TM_0.0_T1_0.7_T2_-0.3_St_1.0_Ra_1.0e+05.jld2"

temp1105 = temp
ls1105 = ls
RB1105 = RB
uvel1105 = uvel
vvel1105 = vvel
num1105 = num

JLD2.@load "/home/tf/Documents/RB_opt/RB_nx_128_ny_32_ratio_4_tend_1.0e+00_TM_0.0_T1_0.7_T2_-0.3_St_1.0_Ra_8.0e+04.jld2"

temp8104 = temp
ls8104 = ls
RB8104 = RB

JLD2.@load "/home/tf/Documents/RB_opt/RB_nx_128_ny_32_ratio_4_tend_1.0e+00_TM_0.0_T1_0.7_T2_-0.3_St_1.0_Ra_4.0e+04.jld2"

temp4104 = temp
ls4104 = ls
RB4104 = RB

JLD2.@load "/home/tf/Documents/RB_opt/RB_nx_128_ny_32_ratio_4_tend_1.0e+00_TM_0.0_T1_0.7_T2_-0.3_St_1.0_Ra_1.0e+04.jld2"

temp1104 = temp
ls1104 = ls
RB1104 = RB

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
lines!(RB1105[1,1:end-1280], RB1105[2,1:end-1280], linewidth = 9, label = "Ra = 1 × 10⁵", color=:red)
lines!(RB8104[1,1:end-1000], RB8104[2,1:end-1000], linewidth = 9, label = "Ra = 8 × 10⁴", color=:purple)
lines!(RB4104[1,1:end-2], RB4104[2,1:end-2], linewidth = 9, label = "Ra = 4 × 10⁴", color=:blue)
lines!(RB1104[1,1:end-2], RB1104[2,1:end-2], linewidth = 9, label = "Ra = 1 × 10⁴", color=:green)
# lines!(RB5103[1,:], RB5103[2,:], linewidth = 9, label = "Ra = 5 × 10³")
resize_to_layout!(fh)
axislegend(position = :rb)
fh = current_figure()

Makie.save("/home/tf/Documents/RB_opt/vel_RB_opt_avheight.png", fh)

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

Makie.save("/home/tf/Documents/RB_opt/vel_RB_opt_effectiveRa.png", fh)


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

Makie.save("/home/tf/Documents/RB_opt/vel_RB_opt_effectiveRa_linear.png", fh)


for i = 2:11
    vort, mag = compute_vorticity_staggered(uvel1105[i,:,:], vvel1105[i,:,:], num.Δ, num.Δ)
    # Set adjustable temperature limits (change these as needed)
    temp_min = -0.3    # e.g. minimum temperature
    temp_max = 0.7 # e.g. maximum temperature
    
    F1 = Figure(size = (3200, 400))
    
    # Create the axis without xlims, ylims, ticklabelsize, or labelfontsize keywords.
    ax = Axis(F1[1, 1],
        aspect = Ratio,
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
        # vort[:,:]'.*(ls1105[i,:,:]' .> 0) .- 0.5, 
        colormap = :dense,
        levels = 30,
        clims = (temp_min, temp_max)
    )
    # hm.colorrange[] = (temp_min, temp_max)
    hm.clims[] = (temp_min, temp_max)


    # contour!(ax, gp.x[1, :], gp.y[:, 1], vort[:,:]'.*(ls1105[i,:,:]' .> 0),
    #     levels = 20,
    #     color = :black,
    #     linewidth = 1.5
    # )
    
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
        limits = (temp_min, temp_max)
    )
    



    temp_min = minimum(vort)    # e.g. minimum temperature
    temp_max = maximum(vort)# e.g. maximum temperature
    print(temp_min, temp_max)
    ax2 = Axis(F1[1, 3],
        aspect = Ratio,
        xlabel = "x", 
        ylabel = "y",
        xticks = -2:0.5:2,
        yticks = [-0.5, 0, 0.5]
    )
    
    # Set the ax2is limits after creation.
    xlims!(ax2, (-2, 2))
    ylims!(ax2, (-0.5, 0.5))
    
    # Increase the font sizes by setting the properties directly.
    ax2.xticklabelsize[] = 30         # size of tick labels on the x ax2is
    ax2.yticklabelsize[] = 30         # size of tick labels on the y ax2is
    ax2.xlabelsize[]     = 34         # size of the x-ax2is label
    ax2.ylabelsize[]     = 34         # size of the y-axis label
    
    # Plot the filled contour (using clims to control the color range)
    hm = contourf!(ax2, gp.x[1, :], gp.y[:, 1],
        # min.(abs.(abs.(opt_S[1]' + opt_L[1]') - abs.(fwdS_des.T[end,:,:]' + fwdL_des.T[end,:,:]')), 1.99), 
        # temp1105[i,:,:]', 
        vort[:,:]'.*(ls1105[i,:,:]' .> 0), 
        colormap = :viridis,
        levels = 10,
        clims = (temp_min, temp_max)
    )
    # hm.colorrange[] = (temp_min, temp_max)
    hm.clims[] = (temp_min, temp_max)


    # contour!(ax, gp.x[1, :], gp.y[:, 1], vort[:,:]'.*(ls1105[i,:,:]' .> 0),
    #     levels = 20,
    #     color = :black,
    #     linewidth = 1.5
    # )
    
    # Add a contour line (for example, at level 0) in red.
    contour!(ax2, gp.x[1, :], gp.y[:, 1], ls1105[i,:,:]',
        levels = [0],
        color = :red,
        linewidth = 5
    )

    # Add a colorbar to the figure next to the axis.
    Colorbar(F1[1, 4];
        colormap = :viridis,
        label = "Vorticity",
        labelsize = 34,
        ticklabelsize = 30,
        width = 40,
        height = 300,
        limits = (temp_min*num.Δ, temp_max*num.Δ)
    )
    # resize_to_layout!(F1)
    F1  # display the figure
    Makie.save("/home/tf/Documents/RB_opt/vort_RB_opt_Ra1105_i$(i-1).png", F1)
end








for i = 2:11
    vort, mag = compute_vorticity_staggered(uvel1105[i,:,:], vvel1105[i,:,:], num.Δ, num.Δ)
    # Set adjustable temperature limits (change these as needed)

    temp_min = minimum(vort)    # e.g. minimum temperature
    temp_max = maximum(vort)# e.g. maximum temperature
    
    F111 = Figure(size = (1600, 400))
    
    # Create the axis without xlims, ylims, ticklabelsize, or labelfontsize keywords.
    ax = Axis(F111[1, 1],
        aspect = Ratio,
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
        # temp1105[i,:,:]', 
        vort[:,:]'.*(ls1105[i,:,:]' .> 0), 
        colormap = :viridis,
        levels = 30,
        clims = (temp_min, temp_max)
    )
    # hm.colorrange[] = (temp_min, temp_max)
    hm.clims[] = (temp_min, temp_max)


    # contour!(ax, gp.x[1, :], gp.y[:, 1], vort[:,:]'.*(ls1105[i,:,:]' .> 0),
    #     levels = 20,
    #     color = :black,
    #     linewidth = 1.5
    # )
    
    # Add a contour line (for example, at level 0) in red.
    contour!(ax, gp.x[1, :], gp.y[:, 1], ls1105[i,:,:]',
        levels = [0],
        color = :red,
        linewidth = 5
    )

    # Add a colorbar to the figure next to the axis.
    Colorbar(F111[1, 2];
        colormap = :viridis,
        label = "Vorticity",
        labelsize = 34,
        ticklabelsize = 30,
        width = 40,
        height = 300,
        limits = (temp_min*num.Δ, temp_max*num.Δ)
    )
    



    
    # resize_to_layout!(F111)
    F111  # display the figure
    Makie.save("/home/tf/Documents/RB_opt/vort_RB_opt_Ra1105_i$(i-1).png", F111)
end









"""
    compute_vorticity_staggered(uvel, vvel, dx, dy)

Compute the 2D vorticity field on a C‑grid:
  • uvel  ∈ ℝ^(nx,   ny+1)
  • vvel  ∈ ℝ^(nx+1, ny)
Returns ω ∈ ℝ^(nx, ny) at cell centers.
Assumes periodic in x, no‑slip (zero) at j=1,ny+1 for uvel,
and zero vvel at j=1,ny boundaries if needed.
"""
function compute_vorticity_staggered(uvel::AbstractMatrix, vvel::AbstractMatrix, dx::Real, dy::Real)
    nx, ny1 = size(uvel)
    nx1, ny   = size(vvel)
    @assert nx1 == nx+1 && ny1 == ny+1

    ω = zeros(eltype(uvel), nx, ny)
    mag = zeros(eltype(uvel), nx, ny)

    # periodic wrap in i-direction
    periodic_i(i) = i < 1  ? nx+1 : i > nx+1 ? 1 : i

    for j in 1:ny, i in 1:nx
        # ∂v/∂x at center i: using vvel[i+1,j] and vvel[i,j]
        im_v = periodic_i(i)
        ip_v = periodic_i(i+1)
        dv_dx = (vvel[ip_v, j] - vvel[im_v, j]) / dx
        vv = (vvel[ip_v, j] - vvel[im_v, j]) / 2.
        # ∂u/∂y at center j: using uvel[i,j+1] and uvel[i,j]
        du_dy = (uvel[i, j+1] - uvel[i, j]) / dy
        uu = (uvel[i, j+1] - uvel[i, j])/ 2.
        w_temp = dv_dx - du_dy
        # if w_temp < 0.
        #     w_temp = max(w_temp, -1000)
        # else
        #     w_temp = min(w_temp, 1000)
        # end
        ω[i, j] = w_temp #dv_dx - du_dy
        mag[i, j] = sqrt(uu*uu + vv*vv)
    end

    return ω, mag
end
