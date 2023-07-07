using CairoMakie


# https://docs.makie.org/stable/examples/plotting_functions/heatmap/

function mandelbrot(x, y)
    z = c = x + y*im
    for i in 1:30.0; abs(z) > 2 && return i; z = z^2 + c; end; 0
end

fig, ax, hm = heatmap(-2:0.1:1, -1.1:0.1:1.1, mandelbrot, colormap = Reverse(:deep))


Colorbar(fig[:, end+1], hm)

fig

f = Figure(fontsize = 20)
Axis(f[1, 1],
    title = "Square grid",
    # $\overline{G_{av}}$ does not work
    titlesize = 30,
    xlabel = L"Iteration step $k$",
    # xlabel = L"Iteration step $k$",
    xlabelsize = 30,
    ylabel = L"$GC$",
    ylabelsize = 30,
    # xtickformat = "{:.1e}"
)
