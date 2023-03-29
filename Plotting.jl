using Plots

Plots.__init__() # workaround with sysimage: https://github.com/JuliaLang/PackageCompiler.jl/issues/786
try
    pgfplotsx()
catch
    println("WARNING: Since you don't have LaTeX with the PGFPlots package installed,")
    println("         the program falls back to using the GR plotting backend.")
    println("         Therefore, the plots may not look fully as intended!")
    ENV["GKSwstype"] = "100" # headless mode (https://discourse.julialang.org/t/plotting-from-a-server/74345/4)
    gr()
end
default(
    minorticks = 10,
    labelfontsize = 12, # default: 11
    legendfontsize = 11, # default: 8
    annotationfontsize = 11, # default: 14
    tickfontsize = 10, # default: 8
    #titlefontsize = 14, # default: 14
    legend_font_halign = :left,
    fontfamily = "Computer Modern",
    framestyle = :box,
)

# a hack to "rasterize" a scatter plot
# TODO: plotting scatter points as rects lowers file size?
function scatterheatmaps!(xss, yss, colors, labels, xlims, ylims; nbins=50, kwargs...)
    xgrid = range(xlims...; length=nbins+1)
    ygrid = range(ylims...; length=nbins+1)
    dx = xgrid[2] - xgrid[1]
    dy = ygrid[2] - ygrid[1]
    zgrid = zeros(Int, nbins, nbins)
    for (i, xs, ys) in zip(1:length(xss), xss, yss)
        for (x, y) in zip(xs, ys)
            ix = min(floor(Int, (x - xgrid[1]) / dx) + 1, nbins)
            iy = min(floor(Int, (y - ygrid[1]) / dy) + 1, nbins)
            zgrid[iy,ix] = i
        end
    end
    xgrid = (xgrid[1:end-1] .+ xgrid[2:end]) ./ 2
    ygrid = (ygrid[1:end-1] .+ ygrid[2:end]) ./ 2
    #zgrid = zgrid .!= 0
    colorgrad = cgrad([:white, colors...], alpha=1.0, categorical=false)
    #histogram2d!(xs, ys; bins=(range(xlims...; length=nbins), range(ylims...; length=nbins)), color = colorgrad, colorbar = :none, kwargs...)
    heatmap!(xgrid, ygrid, zgrid; alpha = 1.0, color = colorgrad, colorbar = :none, kwargs...)
    for (color, label) in zip(colors, labels)
        scatter!([xgrid[1]-1], [ygrid[1]-1]; color = color, markerstrokewidth=0, shape = :rect, label = label) # dummy plot with label
    end
    # TODO: hack a label into the legend with a ghost scatter plot?
end
