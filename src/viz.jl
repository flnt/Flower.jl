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
    v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
)

function make_video(num, fwd, grid, field="u";
                    title_prefix=field, title_suffix="", xlabel="x", ylabel="y", colormap=:viridis,
                    minv=0.0, maxv=0.0, limitsx=false, limitsy=false, framerate=24, step=1, step0=1)
    x = grid.x[1,:]
    y = grid.y[:,1]
    if field == "T"
        z = fwd.Tsave[step0:end,:,:]
        u = fwd.usave[step0:end,:,:]
        colormap = Reverse(:ice)
    elseif field == "u"
        z = fwd.Uxsave[step0:end,:,:]
        u = fwd.uusave[step0:end,:,:]
    elseif field == "v"
        z = fwd.Uysave[step0:end,:,:]
        u = fwd.uvsave[step0:end,:,:]
    elseif field == "p"
        z = fwd.psave[step0:end,:,:].*num.τ
        u = fwd.usave[step0:end,:,:]
    else
        @error ("Possible fields are T, u, v and p")
        return nothing
    end

    if minv == maxv == 0.0
        var_colorrange = true
    else
        var_colorrange = false
    end

    if isa(limitsx, Tuple{Float64,Float64}) || isa(limitsy, Tuple{Int,Int})
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
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    ax  = Axis(fig[1,1], aspect=1, xlabel=xlabel, ylabel=ylabel,
            title=field, xtickalign=0,  ytickalign=0)
    if !var_colorrange
        hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap=colormap, colorrange=(minv, maxv))
    else
        hmap = heatmap!(x, y, @lift(z[$obs,:,:]'), colormap=colormap)
    end
    contour!(x, y, @lift(u[$obs,:,:]'), levels=[0.0], color=:red, linewidth=2.0)
    cbar = fig[1,2] = Colorbar(fig, hmap, labelpadding=0)
    limits!(ax, lx[1], lx[2], ly[1], ly[2])
    resize_to_layout!(fig)

    vid = record(fig, title_prefix*field*"_field"*title_suffix*".mp4", iterator; framerate = framerate) do it
        obs[] = it+1
    end

    return nothing
end