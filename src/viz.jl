struct LogMinorTicks end
	
function MakieLayout.get_minor_tickvalues(
        ::LogMinorTicks, scale, tickvalues, vmin, vmax
)
    vals = Float64[]
    for (lo, hi) in zip(
            @view(tickvalues[1:end-1]),
            @view(tickvalues[2:end])
        )
        interval = hi-lo
        steps = log10.(LinRange(10^lo, 10^hi, 11))
        append!(vals, steps[2:end-1])
    end
    vals
end

custom_formatter(values) = map(
    v -> "10" * Makie.UnicodeFun.to_superscript("$(round(v, digits = 1))"),
    values
)

function plot_grid(grid; linewidth=0.5, limitsx=false, limitsy=false, hide=false, stepx=1, stepy=1)
    x = grid.x_nodes
    y = grid.y_nodes
    if isa(limitsx, Tuple{Float64,Float64}) || isa(limitsx, Tuple{Int,Int})
        lx = limitsx        
    else
        lx = (min(x...), max(x...))
    end
    if isa(limitsy, Tuple{Float64,Float64}) || isa(limitsy, Tuple{Int,Int})
        ly = limitsy
    else
        ly = (min(y...), max(y...))
    end

    fig = Figure(resolution = (1300, 1000))
    ax = Axis(fig[1,1], aspect=DataAspect(), xlabel="x", ylabel="y",
            xtickalign=0,  ytickalign=0)
    # ax.xticks = [-0.5, 0.0, 0.5]
    # ax.yticks = [-0.5, 0.0, 0.5]
    ax.xticks = collect(-15:5:30)
    ax.yticks = collect(-15:5:15)
    if hide
        hidedecorations!(ax, ticks=false, ticklabels=false)
    end
    # for i = 1:stepx:grid.nx+1
    #     lines!(ones(grid.ny+1).*x[i], y, linewidth=linewidth, color=:black)
    # end
    # for i = 1:stepy:grid.ny+1
    #     lines!(x, ones(grid.nx+1).*y[i], linewidth=linewidth, color=:black)
    # end
    x = collect(-15:1:30)
    y = collect(-15:1:15)
    lines!(x, -15*ones(size(x)), color=:black, linewidth = 5)
    lines!(x, 15*ones(size(x)), color=:black, linewidth = 5)
    lines!(-15*ones(size(y)), y, color=:black, linewidth = 5)
    lines!(30*ones(size(y)), y, color=:black, linewidth = 5)

    x = collect(-0.75:0.1:0.75)
    y = collect(-0.75:0.1:0.75)
    lines!(x, -0.75*ones(size(x)), color=:black, linewidth = 3)
    lines!(x, 0.75*ones(size(x)), color=:black, linewidth = 3)
    lines!(-0.75*ones(size(y)), y, color=:black, linewidth = 3)
    lines!(0.75*ones(size(y)), y, color=:black, linewidth = 3)

    contour!(grid.x[1,:], grid.y[:,1], grid.u', levels=[0.0], color=:red, linewidth=2.0)
    limits!(ax, lx[1], 33, ly[1], ly[2])
    # resize_to_layout!(fig)

    return fig
end

function make_video(num, fwd, grid, field="u";
                    title_prefix=field, title_suffix="", xlabel="x", ylabel="y", colormap=:viridis,
                    minv=0.0, maxv=0.0, limitsx=false, limitsy=false, framerate=24, step=1, step0=1, stepf=size(fwd.usave,1))
    x = grid.x[1,:]
    y = grid.y[:,1]
    plot_hmap = true
    if field == "T"
        z = fwd.Tsave[step0:stepf,:,:]
        u = fwd.usave[step0:stepf,:,:]
        colormap = Reverse(:ice)
    elseif field == "u"
        z = fwd.Uxsave[step0:stepf,:,:]
        u = fwd.uusave[step0:stepf,:,:]
    elseif field == "v"
        z = fwd.Uysave[step0:stepf,:,:]
        u = fwd.uvsave[step0:stepf,:,:]
    elseif field == "ucorr"
        z = fwd.Uxcorrsave[step0:stepf,:,:]
        u = fwd.uusave[step0:stepf,:,:]
    elseif field == "vcorr"
        z = fwd.Uycorrsave[step0:stepf,:,:]
        u = fwd.uvsave[step0:stepf,:,:]
    elseif field == "p"
        z = fwd.psave[step0:stepf,:,:].*num.τ
        u = fwd.usave[step0:stepf,:,:]
    elseif field == "ϕ"
        z = fwd.ϕsave[step0:stepf,:,:].*num.τ
        u = fwd.usave[step0:stepf,:,:]
    elseif field == "κ"
        z = fwd.κsave[step0:stepf,:,:]
        u = fwd.usave[step0:stepf,:,:]
    else
        plot_hmap = false
        z = fwd.usave[step0:stepf,:,:]
        u = fwd.usave[step0:stepf,:,:]
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64,Float64}) || isa(limitsx, Tuple{Int,Int})
        lx = limitsx
    else
        lx = (min(x...), max(x...))
    end
    if isa(limitsy, Tuple{Float64,Float64}) || isa(limitsy, Tuple{Int,Int})
        ly = limitsy
    else
        ly = (min(y...), max(y...))
    end

    obs = Observable{Int32}(1)
    iterator = range(0, size(z, 1)-1, step=step)

    fontsize_theme = Theme(fontsize = 30)
    set_theme!(fontsize_theme)

    fig = Figure(resolution = (1600, 1000))
    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    ax  = Axis(fig[1,1], aspect=DataAspect(), xlabel=xlabel, ylabel=ylabel,
            title=field, xtickalign=0,  ytickalign=0)
    if plot_hmap
        if !var_colorrange
            hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap=colormap, colorrange=(minv, maxv))
        else
            hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap=colormap)
        end
    end
    contour!(x, y, @lift(u[$obs,:,:]'), levels=[0.0], color=:red, linewidth=2.0)
    if plot_hmap
        cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding=0)
    end
    limits!(ax, lx[1], lx[2], ly[1], ly[2])
    colgap!(fig.layout, 10)
    rowgap!(fig.layout, 10)
    resize_to_layout!(fig)

    vid = record(fig, title_prefix*field*"_field"*title_suffix*".mp4", iterator; framerate = framerate) do it
        obs[] = it+1
    end

    return nothing
end