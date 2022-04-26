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