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
    if hide
        hidedecorations!(ax, ticks=false, ticklabels=false)
    end
    for i = 1:stepx:grid.nx+1
        lines!(ones(grid.ny+1).*x[i], y, linewidth=linewidth, color=:black)
    end
    for i = 1:stepy:grid.ny+1
        lines!(x, ones(grid.nx+1).*y[i], linewidth=linewidth, color=:black)
    end

    contour!(grid.x[1,:], grid.y[:,1], grid.u', levels=[0.0], color=:red, linewidth=2.0)
    limits!(ax, lx[1], lx[2], ly[1], ly[2])
    # resize_to_layout!(fig)

    return fig
end

function make_video(grid, field_u, field=nothing;
                    title_prefix="video", title_suffix="", xlabel="x", ylabel="y", colormap=:viridis,
                    minv=0.0, maxv=0.0, limitsx=false, limitsy=false, var=1,
                    framerate=24, step=1, step0=1, stepf=size(field_u,1))
    x = grid.x[1,:]
    y = grid.y[:,1]

    u = field_u[step0:stepf,:,:]
    plot_hmap = true
    if isnothing(field)
        plot_hmap = false
    else
        if length(size(field)) == 2
            if var==1
                z = reshape(field[step0:stepf,1:grid.ny*grid.nx], (stepf-step0+1,grid.ny,grid.nx))
            else
                z = reshape(field[step0:stepf,grid.ny*grid.nx+1:end], (stepf-step0+1,grid.ny,grid.nx))
            end
        else
            z = field[step0:stepf,:,:]
        end
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64,Float64}) || isa(limitsx, Tuple{Int,Int})
        lx = limitsx
    else
        lx = (min(x...)-grid.dx[1,1]/2, max(x...)+grid.dx[1,end]/2)
    end
    if isa(limitsy, Tuple{Float64,Float64}) || isa(limitsy, Tuple{Int,Int})
        ly = limitsy
    else
        ly = (min(y...)-grid.dy[1,1]/2, max(y...)+grid.dy[end,1]/2)
    end

    obs = Observable{Int32}(1)
    iterator = range(0, size(u, 1)-1, step=step)

    fontsize_theme = Theme(fontsize = 30)
    set_theme!(fontsize_theme)

    fig = Figure(resolution = (1600, 1000))
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    ax  = Axis(fig[1,1], aspect=DataAspect(), xlabel=xlabel, ylabel=ylabel,
            xtickalign=0,  ytickalign=0)
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

    vid = record(fig, title_prefix*title_suffix*".mp4", iterator; framerate = framerate) do it
        obs[] = it+1
    end

    return nothing
end