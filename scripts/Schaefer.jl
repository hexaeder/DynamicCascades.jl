using DynamicCascades
using LightGraphs
using MetaGraphs
using Statistics
using GLMakie
using GraphMakie
# using CairoMakie

DIR = "/Users/hw/MA/Forschungsbeleg/figures/"
set_theme!(theme_minimal(), fontsize=20)

orange = Makie.RGB([227, 114, 34]./255...)
gray = Makie.RGB([142, 144, 143]./255...)
cyan = Makie.RGB([0, 159, 218]./255...)
green = Makie.RGB([105, 146, 58]./255...)
blueish = Makie.RGB([124, 174, 175]./255...)
####
#### using two steadystates
####
network = import_system(:schaefer2018; γ=0.8)
sol = simulate(network;
               initial_fail=[5],
               trip_lines=false,
               tspan=(0., 1000.),
               solverargs=(;dtmax=0.1));
# inspect_solution(sol)
t1 = 0.0
t2 = sol.sol.t[end]

node_size = [40, 50, 40, 40, 50]
node_color = [blueish, orange, blueish, blueish, orange]
gpargs = gparguments(sol, t1; colortype=:abssteady)
fig, ax, p = graphplot(network; gpargs..., node_size, node_color)
ax.aspect = DataAspect()
p.edge_width[] = 2 .* p.edge_width[]
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR, "schaefer_1.png"), fig)

gpargs = gparguments(sol, t2; colortype=:abssteady)
fig, ax, p = graphplot(network; gpargs..., node_size, node_color)
ax.aspect = DataAspect()
p.edge_width[] = 2 .* p.edge_width[]
p.edge_color[][5] = Makie.RGB(1,1,1)
notify(p.edge_color)
hidedecorations!(ax); hidespines!(ax)
save(joinpath(DIR, "schaefer_2.png"), fig)


####
#### Transient without failure
####
function create_animation(sol, tmin, tmax, file; showrating=false)
    fig = Figure(resolution=(1600,600))
    fig[1,1] = nwax = Axis(fig, width=500, tellheigth=false, aspect=DataAspect())
    fig[1,2] = flowax = Axis(fig, heigth=600, xlabel="Time (s)", ylabel="absolut flow in PU")

    t = Node(0.0)
    rating = get_prop(network, Edge(1,2), :rating)
    hlines!(flowax, rating, linewidth=3, color=gray, visible=showrating)
    gpargs = gparguments(sol, t;
                         colortype=:abssteady,
                         offlinecolor=Makie.RGB(1,1,1))
    p = graphplot!(nwax, network; gpargs..., node_size=30, node_color)
    hidedecorations!(nwax); hidespines!(nwax)
    for idx in 1:ne(network)
        (ts, S) = seriesforidx(sol.load_S, idx)
        lines!(flowax, ts, S, linewidth=5)
    end
    hlines!(flowax, rating, linewidth=5, color=gray, visible=showrating)
    vlines!(flowax, t, linewidth=3, color=gray, visible=@lift($t>tmin))
    xlims!(flowax, tmin, tmax)
    ylims!(flowax, 0, 1.1)

    time, fps = 10, 30
    tspan = range(tmin, tmax, length=fps*time)
    save(joinpath(DIR, "../videos", file*".png"), fig)
    record(fig, joinpath(DIR, "../videos", file*".mp4"), tspan) do time
        t[] = time
    end
end

sol = simulate(network;
               initial_fail=[5],
               trip_lines=false,
               tspan=(0., 1000.),
               solverargs=(;dtmax=0.1));
create_animation(sol, 0, 15, "schaefer_transient"; showrating=true)

sol = simulate(network;
               initial_fail=[5],
               trip_lines=true,
               tspan=(0., 1000.),
               solverargs=(;dtmax=0.1));
create_animation(sol, 0, 7, "schaefer_transient_fail"; showrating=true)
