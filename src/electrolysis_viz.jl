


"""
From plot_grid
"""
function plot_grid_fig!(fig,ax,
    num, grid;
    linewidth = 0.5,
    limitsx = false,
    limitsy = false,
    hide = false,
    skipx = 0,
    skipy = 0,
    xscale =1.0,
    yscale = 1.0,
    )

    x = grid.x_nodes ./xscale
    y = grid.y_nodes ./yscale
  
    for i = 1:skipx+1:grid.nx+1
        lines!(ones(grid.ny+1) .* x[i], y, linewidth = linewidth, color = :black)
    end
    for i = 1:skipy+1:grid.ny+1
        lines!(x, ones(grid.nx+1) .* y[i], linewidth = linewidth, color = :black)
    end

end


"""
    make_video(
        grid, field_u, field = nothing;
        title_prefix = "video",
        title_suffix = "",
        xlabel = L"x",
        ylabel = L"y",
        colormap = :viridis,
        minv = 0.0,
        maxv = 0.0,
        limitsx = false,
        limitsy = false,
        var = 1,
        framerate = 24,
        step = 1,
        step0 = 1,
        stepf = size(field_u, 1)
    )

Generates and saves a video displaying a `field` and an interface `field_u`.
"""
function make_video_vec(
    num, grid, field_u, field = nothing,
     vecu = nothing, vecv=nothing;
    title_prefix = "video",
    title_suffix = "",
    xlabel = L"x",
    ylabel = L"y",
    colormap = :viridis,
    sz = (1600, 1000),
    minv = 0.0,
    maxv = 0.0,
    limitsx = false,
    limitsy = false,
    var = 1,
    framerate = 24,
    step = 1,
    step0 = 1,
    stepf = size(field_u, 2),
    xscale = 1, 
    yscale = 1, 
    xticks = nothing,
    yticks = nothing,
    scalscale = 1, 
    scalelabel = nothing,
    scalticks = nothing,
    )

    x = grid.x[1,:] ./xscale
    y = grid.y[:,1] ./yscale

    u = field_u[:,step0:stepf,:,:]
    plot_hmap = true
    if isnothing(field)
        plot_hmap = false
    else
        if length(size(field)) == 2
            if var == 1
                z = reshape(
                    field[step0:stepf, 1:grid.ny*grid.nx], 
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            else
                z = reshape(
                    field[step0:stepf, grid.ny*grid.nx+1:2*grid.ny*grid.nx],
                    (stepf-step0+1, grid.ny,grid.nx)
                )
            end
        else
            z = field[step0:stepf,:,:]
        end

        z = z ./scalscale
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64, Float64}) || isa(limitsx, Tuple{Int, Int})
        lx = limitsx
    else
        lx = (min(x...) - grid.dx[1,1] / 2, max(x...) + grid.dx[1,end] / 2)
    end
    if isa(limitsy, Tuple{Float64, Float64}) || isa(limitsy, Tuple{Int, Int})
        ly = limitsy
    else
        ly = (min(y...) - grid.dy[1,1] / 2, max(y...) + grid.dy[end,1] / 2)
    end

    obs = Observable{Int32}(1)
    iterator = range(0, size(u, 2) - 1, step=step)

    # fig = Figure(size = sz)
    fig = Figure()
    ax  = Axis(fig[1,1], 
    # aspect=DataAspect(), 
    aspect=1,
    xlabel = xlabel, ylabel = ylabel,
    xticks=xticks, yticks=yticks,
        xtickalign = 0,  ytickalign = 0)
    if plot_hmap
        print("scalticks", scalticks, "xticks", xticks, "yticks", yticks,"\n")

        if !var_colorrange
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap, 
                colorrange = (minv, maxv))
        else
            # if !isnothing(scalticks)
            #     print("scalticks", scalticks)

            #     hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap = colormap, ticks=scalticks)
            # else
            hmap = heatmap!(ax, x, y, @lift(z[$obs,:,:]'), colormap = colormap)
            # end
        end
    end
    if !plot_hmap
        contour!(x, y, @lift(u[1,$obs,:,:]'), levels = -10:grid.dx[1,1]:10, linewidth = 2.0)
    end
    for iLS in 1:num.nLS
        contour!(x, y, @lift(u[iLS,$obs,:,:]'), levels = [0.0], color = :red, linewidth = 3.0)
    end

    #TODO
    if !isnothing(vecu)
        arrows!(grid.x[1,:], grid.y[:,1], vecu[1,:,:], vecv[1,:,:], 
        arrowsize = 10, 
        lengthscale = 0.3,
        # arrowcolor = strength, 
        # linecolor = strength,
        )
    end

    # contourf!(x, y, u[2,1,:,:]', levels = 0:0, color = :red, linewidth = 3, extendlow = :gray);
    if plot_hmap
        if !isnothing(scalticks)
            print("scalticks", scalticks)
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0, ticks=scalticks)
            
            # Colorbar(fig[1, 2], hmap, 
            # label = "Reverse sequential colormap",
            # ticks=scalticks)
            # cbar.ticks = (scalticks)
            #https://stackoverflow.com/questions/73555826/modifying-label-and-tick-on-color-bar-of-plots-jl-plot

            Colorbar(fig[1, 2], hmap, ticks = -1:0.25:1)

        else
            cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding = 0)
        end

    end

    if !isnothing(scalelabel)
        # Label(fig[1, 1, Top()], halign = :right, scalelabel)
        Label(fig[1, 2, Top()], halign = :center, scalelabel)
    end

    # limits!(ax, lx[1], lx[2], ly[1], ly[2])
    # colgap!(fig.layout, 10)
    # rowgap!(fig.layout, 10)
    # colsize!(fig.layout, 1, widths(ax.scene.viewport[])[1])
    # rowsize!(fig.layout, 1, widths(ax.scene.viewport[])[2])
    # resize_to_layout!(fig)

    vid = record(
            fig, title_prefix*title_suffix*".mp4", iterator;
            framerate = framerate
        ) do it
        obs[] = it+1
    end

    return vid
end
