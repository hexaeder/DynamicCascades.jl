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

################################# phase angle ##################################
I = 1
D = 0.5
K = 1
P0 = 1
P1 = 2
Z = 4*I*K - D^2

# w= sqrt(4*I*K-D^2) / (2*I)
# g = D / (2*I)
# β = (PO - P1) / (K*sin(atan(w / g)))
# α = atan(w / g)


function phase_nadir(I, D; K=1, P0=1, P1=2)
    abs(P1 - P0) / K * (1 + exp((-D*π) / sqrt(4*I*K - D^2)))
end

fig, ax, hm = Makie.heatmap(0:0.1:2, -0:0.1:1, phase_nadir, colormap = Reverse(:deep))

Colorbar(fig[:, end+1], hm)
fig

################################# frequency ####################################
I = 1
D = 0.5
K = 1
P0 = 1
P1 = 2
Z = 4*I*K - D^2

# w= sqrt(4*I*K-D^2) / (2*I)
# g = D / (2*I)
# β = (PO - P1) / (K*sin(atan(w / g)))
# α = atan(w / g)

function frequency_nadir(I, D; K=1, P0=1, P1=2)
    abs(P1 - P0) / K * (1 + exp((-D*π) / sqrt(4*I*K - D^2)))
end


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
