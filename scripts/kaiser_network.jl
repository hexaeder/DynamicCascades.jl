using Serialization
using DynamicCascades
using LightGraphs
using MetaGraphs
using ProgressMeter
using DataFrames
using Statistics
using GLMakie
set_theme!(resolution=(1600, 1200))
set_theme!(resolution=(1000, 800))

####
#### create Data
####
dir_ins = joinpath(RAWRESULTS_DIR, timestamp()*"_kaiser_insulator")
mkpath(dir_ins)
damping = 1.0
loadt = 0.1
@showprogress for seed in 1:20
    network = import_system(:kaiser2020; gen_γ=damping, load_τ=loadt, seed)
    x_static = steadystate(network);

    Threads.@threads for i in 1:ne(network)
        filename = joinpath(dir_ins, pstring(;seed=lpad(seed,2,"0"), fail=lpad(i,2,"0"), damping, loadt)*".dat")
        simulate(network;
                 x_static,
                 tspan=(0.0, 1000.),
                 initial_fail=[i],
                 trip_lines=false,
                 filename, verbose=false)
    end
end

dir_woi = joinpath(RAWRESULTS_DIR, timestamp()*"_kaiser_wo_insulator")
mkpath(dir_woi)
@showprogress for seed in 1:20
    network = import_system(:kaiser2020; gen_γ=damping, load_τ=loadt, seed)
    rem_edge!(network, 2, 22)
    rem_edge!(network, 1, 21)
    set_prop!(network, 1, 22, :X, 0.005)
    set_prop!(network, 2, 21, :X, 0.005)
    x_static = steadystate(network);

    Threads.@threads for i in 1:ne(network)
        filename = joinpath(dir_woi, pstring(;seed=lpad(seed,2,"0"), fail=lpad(i,2,"0"), damping, loadt)*".dat")
        simulate(network;
                 x_static,
                 tspan=(0.0, 1000.),
                 initial_fail=[i],
                 trip_lines=false,
                 filename, verbose=false)
    end
end

####
#### load data
####
dir_ins = "/Users/hw/MAScratch/2021-06-21_kaiser_insulator"
df_ins = DataFrame()
@showprogress for file in readdir(dir_ins)
    sol = deserialize(joinpath(dir_ins, file))

    failedge = collect(edges(sol.network))[sol.initial_fail[begin]]
    region = get_prop(sol.network, failedge, :region)
    region == 0 && continue # don't care about fails in the insulator right now

    subdf = maxload_dist_time(sol)
    subdf[!, :seed] .= get_prop(sol.network, :seed)
    subdf[!, :fail] .= sol.initial_fail[begin]

    append!(df_ins, subdf, cols=:union)
end

dir_woi = "/Users/hw/MAScratch/2021-06-21_kaiser_wo_insulator"
df_woi = DataFrame()
@showprogress for file in readdir(dir_woi)
    sol = deserialize(joinpath(dir_woi, file))

    failedge = collect(edges(sol.network))[sol.initial_fail[begin]]
    region = get_prop(sol.network, failedge, :region)
    region == 0 && continue # don't care about fails in the insulator right now

    subdf = maxload_dist_time(sol)
    subdf[!, :seed] .= get_prop(sol.network, :seed)
    subdf[!, :fail] .= sol.initial_fail[begin]

    append!(df_woi, subdf, cols=:union)
end

same_ins = df_ins[df_ins.region .== :same, :]
other_ins = df_ins[df_ins.region .== :other, :]
same_woi = df_woi[df_woi.region .== :same, :]
other_woi = df_woi[df_woi.region .== :other, :]


####
#### Plot flow density for different seeds
####
fig = Figure()
fig[1,1] = ax = Axis(fig)
for group in groupby(df_ins, :seed)
    Makie.density!(ax, group.staticA)
end
fig = Figure()
fig[1,1] = ax = Axis(fig)
for group in groupby(df_woi, :seed)
    Makie.density!(ax, group.staticA)
end
Makie.hist(same_ins.d; bins=[i+0.5 for i in minimum(same_ins.d)-1:maximum(same_ins.d)])
Makie.hist!(other_ins.d; bins=[i+0.5 for i in minimum(other_ins.d)-1:maximum(other_ins.d)])

Makie.density(same_ins.absdiff)
Makie.density!(other_ins.absdiff)

Makie.density(same_ins.reldiff)
Makie.density!(other_ins.reldiff)
xlims!(0,5)

####
#### Plot absolut difference from initial state
####
# CairoMakie.activate!()
fig = Figure()
fig[1,1] = ax = Axis(fig,
                     limits=(nothing, (0,0.1)),
                     xlabel="Distance between edges",
                     ylabel="Max overload (absolut difference to steady state in PU)",
                     xticks=1:9,
                     title="Mean overload and distribution over distance")
boundary = (0.0, 0.5)
npoints = 1000
width=1.4
col = ax.palette[:color][]
v1 = violin!(ax, same_ins.d, same_ins.absdiff; side=:left, boundary, npoints, width, show_median=false, color=(col[1], 0.3))
v2 = violin!(ax, other_ins.d, other_ins.absdiff; side=:right, boundary, npoints, width, show_median=false, color=(col[2], 0.3))
v3 = violin!(ax, other_woi.d, other_woi.absdiff; side=:right, boundary, npoints, width, show_median=false, color=(col[3], 0.3))
smean = [Point2f0(g.d[begin], mean(g.absdiff)) for g in groupby(same_ins, :d)]
omean = [Point2f0(g.d[begin], mean(g.absdiff)) for g in groupby(other_ins, :d)]
womean = [Point2f0(g.d[begin], mean(g.absdiff)) for g in groupby(other_woi, :d)]
linewidth=5
markersize=20
l1 = lines!(ax, smean; linewidth); # scatter!(ax, smean; markersize)
l2 = lines!(ax, omean; linewidth); # scatter!(ax, omean; markersize)
l3 = lines!(ax, womean; linewidth); # scatter!(ax, omean; markersize)

fig[1,1] = Legend(fig,
                  [[v1,l1], [v2,l2], [v3,l3]],
                  ["same region",
                   "other region (insulator)",
                   "other region (w/o insulator)"],
                  tellwidth=false,
                  halign=:right, valign=:top)
save("over_load_over_dist.pdf", fig)

####
#### Plot relative difference from initial state
####
smask = abs.(same_ins.staticA) .> 0.2
omask = abs.(other_ins.staticA) .> 0.2
scatter(same_ins.absdiff, same_ins.reldiff; color=same_ins.staticA)
scatter(same_ins.absdiff[smask], same_ins.reldiff[smask]; color=same_ins.staticA)
boundary = (0.0, 10)
npoints = 2000
width = 1.0
violin(same_ins.d[smask], same_ins.absdiff[smask]; side=:left, boundary, npoints, width, show_median=false)
violin!(other_ins.d[omask], other_ins.absdiff[omask]; side=:right, boundary, npoints, width, show_median=false)
ylims!(0,0.5)
scatterlines!([Point2f0(g.d[begin], mean(g.absdiff)) for g in groupby(same_ins[smask,:], :d)])
scatterlines!([Point2f0(g.d[begin], mean(g.absdiff)) for g in groupby(other_ins[omask,:], :d)])

####
#### Plot steady state difference
####
scatter(same_ins.d.-0.05, (same_ins.staticA.-same_ins.staticB))
scatter!(other_ins.d.+0.05, (other_ins.staticA.-other_ins.staticB))
npoints=1000
violin!(same_ins.d.-0.1, (same_ins.staticA.-same_ins.staticB); side=:left, width=2.0)

scatter!(same_woi.d.-0.10, (same_woi.staticA.-same_woi.staticB))
scatter!(other_woi.d.+0.05, (other_woi.staticA.-other_woi.staticB))

GLMakie.activate!()
####
#### Plot time over distance
####
scatter(same_ins.d, same_ins.t)
scatter!(other_ins.d, other_ins.t)
scatter!(other_woi.d, other_woi.t)

####
#### Simulate specific failure
####
network = import_system(:kaiser2020; gen_γ=damping, load_τ=loadt, seed)
sol = simulate(network; tspan=(0.0, 1000.),
         initial_fail=[38], trip_lines=false);

CairoMakie.activate!()
t = Observable(0.0)
fig = Figure()
fig[1,1] = ax = Axis(fig)
graphplot!(ax, sol, t; colortype=:abssteady)
hidedecorations!(ax); hidespines!(ax)
fig[1,1] = Label(fig, @lift("time = " * lpad(round($t, digits=2), 5) * " s"), halign=:left, valign=:top, tellwidth=false, tellheight=false)

framerate = 25
time = 10
t1, t2 = 0.0, 10.0
record(fig, "animation.mp4", range(t1, t2, length=framerate*time); framerate) do time
    t[] = time
end
