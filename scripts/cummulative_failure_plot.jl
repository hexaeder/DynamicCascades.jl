###
### Cummulative failures plot
###


"""
GPT generated Cumulative failures (node + line), narrow=solid, wide=dashed.
- Ignores type (node/line) -> just counts events.
- Removes the FIRST line failure (per your request).
"""

function cumulative_failure_plot(sol1, sol4;
        resolution=(900,600),
        lw=4,
        xlim=:auto)

    # --- collect times ---
    # narrow
    t_nodes_n = collect(sol1.failures_nodes.t)
    t_lines_n = collect(sol1.failures.t)
    !isempty(t_lines_n) && popfirst!(t_lines_n)  # remove first line failure
    t_all_n = sort(vcat(t_nodes_n, t_lines_n))

    # wide
    t_nodes_w = collect(sol4.failures_nodes.t)
    t_lines_w = collect(sol4.failures.t)
    !isempty(t_lines_w) && popfirst!(t_lines_w)  # remove first line failure
    t_all_w = sort(vcat(t_nodes_w, t_lines_w))


    # --- choose xlim automatically if requested ---
    if xlim === :auto
        max_t = maximum(vcat(t_all_n, t_all_w); init=0.0)
        xlim = (0.0, max_t * 1.02 + 1e-6)   # small padding
    end

    # --- build step series ---
    x_n = vcat(xlim[1], t_all_n)
    y_n = 0:length(t_all_n)

    x_w = vcat(xlim[1], t_all_w)
    y_w = 0:length(t_all_w)

    # --- plot ---
    fig = Figure(resolution=resolution, figure_padding=30)
    ax  = Axis(fig[1,1]; xlabel="Time (s)", ylabel="Cumulative failures")

    stairs!(ax, x_n, y_n; color=:blue, linewidth=lw, linestyle=:solid, label="narrow")
    stairs!(ax, x_w, y_w; color=:black, linewidth=lw, linestyle=:solid,  label="wide")

    xlims!(ax, xlim...)
    ylims!(ax, 0, max(length(t_all_n), length(t_all_w)) + 1)

    axislegend(ax; position=:lt, framevisible=false)

    return fig
end



fig = cumulative_failure_plot(sol1, sol4; xlim=(0, 75))
save(joinpath(MA_DIR, string("WS/illustrative_plot/cummulative_failures.svg")), fig)