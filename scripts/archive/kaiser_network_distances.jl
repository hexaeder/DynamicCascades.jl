using Serialization
using DynamicCascades
using Graphs
using MetaGraphs
using ProgressMeter
using DataFrames
using Statistics
using GLMakie

DIR = "/Users/hw/MA/Forschungsbeleg/figures/"
orange = Makie.RGB([227, 114, 34]./255...)
gray = Makie.RGB([142, 144, 143]./255...)
cyan = Makie.RGB([0, 159, 218]./255...)
green = Makie.RGB([105, 146, 58]./255...)
blueish = Makie.RGB([124, 174, 175]./255...)

set_theme!(theme_minimal(), fontsize=20)
set_theme!(resolution=(800, 600))

####
#### create Data
####
dir1 = joinpath(RESULTS_DIR, timestamp()*"_kaiser_insulator")
mkpath(dir1)
damping = 1.0
loadt = 0.1
@showprogress for seed in 1:20
    network = import_system(:kaiser2020; gen_γ=damping, load_τ=loadt, seed)
    x_static = steadystate(network);

    Threads.@threads for i in 1:ne(network)
        filename = joinpath(dir1, pstring(;seed=lpad(seed,2,"0"), fail=lpad(i,2,"0"), damping, loadt)*".dat")
        simulate(network;
                 x_static,
                 tspan=(0.0, 1000.),
                 initial_fail=[i],
                 trip_lines=false,
                 filename, verbose=false)
    end
end

dir2 = joinpath(RESULTS_DIR, timestamp()*"_kaiser_wo_insulator")
mkpath(dir2)
@showprogress for seed in 1:20
    network = import_system(:kaiser2020; gen_γ=damping, load_τ=loadt, seed)
    rem_edge!(network, 2, 22)
    rem_edge!(network, 1, 21)
    set_prop!(network, 1, 22, :X, 0.005)
    set_prop!(network, 2, 21, :X, 0.005)
    x_static = steadystate(network);

    Threads.@threads for i in 1:ne(network)
        filename = joinpath(dir2, pstring(;seed=lpad(seed,2,"0"), fail=lpad(i,2,"0"), damping, loadt)*".dat")
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
dir1 = "/Users/hw/MAScratch/2021-06-21_kaiser_insulator"
df1 = DataFrame()
@showprogress for file in readdir(dir1)
    sol = deserialize(joinpath(dir1, file))

    failedge = collect(edges(sol.network))[sol.initial_fail[begin]]
    region = get_prop(sol.network, failedge, :region)
    region == 0 && continue # don't care about fails in the insulator right now

    subdf = maxload_dist_time(sol)
    subdf[!, :seed] .= get_prop(sol.network, :seed)
    subdf[!, :fail] .= sol.initial_fail[begin]

    append!(df1, subdf, cols=:union)
end

dir2 = "/Users/hw/MAScratch/2021-06-21_kaiser_wo_insulator"
df2 = DataFrame()
@showprogress for file in readdir(dir2)
    sol = deserialize(joinpath(dir2, file))

    failedge = collect(edges(sol.network))[sol.initial_fail[begin]]
    region = get_prop(sol.network, failedge, :region)
    region == 0 && continue # don't care about fails in the insulator right now

    subdf = maxload_dist_time(sol)
    subdf[!, :seed] .= get_prop(sol.network, :seed)
    subdf[!, :fail] .= sol.initial_fail[begin]

    append!(df2, subdf, cols=:union)
end

same1 = df1[df1.region .== :same, :]
other1 = df1[df1.region .== :other, :]
same2 = df2[df2.region .== :same, :]
other2 = df2[df2.region .== :other, :]

####
#### Plot flow density for different seeds
####
fig = Figure()
fig[1,1] = ax = Axis(fig)
for group in groupby(df1, :seed)
    Makie.density!(ax, group.staticA)
end
fig = Figure()
fig[1,1] = ax = Axis(fig)
for group in groupby(df2, :seed)
    Makie.density!(ax, group.staticA)
end
Makie.hist(same1.d; bins=[i+0.5 for i in minimum(same1.d)-1:maximum(same1.d)])
Makie.hist!(other1.d; bins=[i+0.5 for i in minimum(other1.d)-1:maximum(other1.d)])

Makie.density(same1.absdiff)
Makie.density!(other1.absdiff)

Makie.density(same1.reldiff)
Makie.density!(other1.reldiff)
xlims!(0,5)

####
#### Plot absolut difference from initial state
####
CairoMakie.activate!()
fig = Figure()
fig[1,1] = ax = Axis(fig,
                     limits=((0.1, 9.9), (0,0.11)),
                     xlabel="Distance between edges",
                     ylabel="Max overload (difference to steady state in PU)",
                     xticks=1:9,
                     title="Mean overload and distribution over distance")
# boundary = (0.0, 0.5)
npoints = 1000
col = ax.palette[:color][]
# v1 = violin!(ax, same1.d, same1.absdiff; side=:left, boundary, npoints, width=1.4, show_median=false, color=(col[1], 0.3),)
# v2 = violin!(ax, other1.d, other1.absdiff; side=:right, boundary, npoints, width=1.4, show_median=false, color=(col[2], 0.3),)
# v3 = violin!(ax, other2.d, other2.absdiff; side=:right, boundary, npoints, width=1.4, show_median=false, color=(col[3], 0.3),)

v1 = violin!(ax, same1.d, same1.absdiff; side=:left, npoints, width=1.4, show_median=false, color=(col[1], 0.4),)
v2 = violin!(ax, other1.d, other1.absdiff; side=:right, npoints, width=1.4, show_median=false, color=(col[2], 0.4),)
# v3 = violin!(ax, other2.d, other2.absdiff; side=:right, npoints, width=1.4, show_median=false, color=(col[3], 0.3),)

g = groupby(same1, :d; sort=true)

s_d =  [g.d[1] for g in groupby(same1, :d; sort=true)]
o_d =  [g.d[1] for g in groupby(other1, :d; sort=true)]
wo_d = [g.d[1] for g in groupby(other2, :d; sort=true)]

s_mean =  [mean(g.absdiff) for g in groupby(same1, :d; sort=true)]
o_mean =  [mean(g.absdiff) for g in groupby(other1, :d; sort=true)]
wo_mean = [mean(g.absdiff) for g in groupby(other2, :d; sort=true)]

s_median =  [median(g.absdiff) for g in groupby(same1, :d; sort=true)]
o_median =  [median(g.absdiff) for g in groupby(other1, :d; sort=true)]
wo_median = [median(g.absdiff) for g in groupby(other2, :d; sort=true)]

# s_std =  [std(g.absdiff) for g in groupby(same1, :d; sort=true)]
# o_std =  [std(g.absdiff) for g in groupby(other1, :d; sort=true)]
# wo_std = [std(g.absdiff) for g in groupby(other2, :d; sort=true)]

# band!(ax, s_d, s_mean+s_std, s_mean-s_std)
# band!(ax, o_d, o_mean+o_std, o_mean-o_std)
# band!(ax, wo_d, wo_mean+wo_std, wo_mean-wo_std)

linewidth=3
markersize=15
l1 = lines!(ax, s_d, s_mean; linewidth, label="same region");
s1 = scatter!(ax, s_d, s_mean; markersize);

l2 = lines!(ax, o_d,o_mean; linewidth, label="other region");
s2 = scatter!(ax, o_d, o_mean; markersize);
axislegend(framevisible=false)

save(joinpath(DIR, "load_over_dist_1.pdf"), fig)

l3 = lines!(ax, wo_d, wo_mean; linewidth, label="other region (w/o insulator)");
s3 = scatter!(ax, wo_d, wo_mean; markersize);

save(joinpath(DIR, "load_over_dist_2.pdf"), fig)

# fig[1,2] = Legend(fig,
#                   [[v1,l1, s1], [v2,l2], [v3,l3]],
#                   ["same region",
#                    "other region (insulator)",
#                    "other region (w/o insulator)"],
#                   tellwidth=false,
#                   halign=:right, valign=:top)
# save("over_load_over_dist.pdf", fig)

####
#### Plot absolut difference over Reistance differenc
####
# CairoMakie.activate!()
same1.dbin = (d -> round(d*2)/2).(same1.dres)
other1.dbin = (d -> round(d*2)/2).(other1.dres)
other2.dbin = (d -> round(d*2)/2).(other2.dres)

fig = Figure()
fig[1,1] = ax = Axis(fig,
                     limits=((0.0, 4.0), (0,0.11)),
                     xlabel="resistance distance between edges",
                     ylabel="Max overload (difference to steady state in PU)",
                     xticks=1:9,
                     title="Mean overload and distribution over distance")
# boundary = (0.0, 0.5)
npoints = 1000
col = ax.palette[:color][]
v1 = violin!(ax, same1.dbin, same1.absdiff; side=:left, npoints, width=0.4, show_median=false, color=(col[1], 0.4),)
v2 = violin!(ax, other1.dbin, other1.absdiff; side=:right, npoints, width=0.4, show_median=false, color=(col[2], 0.4),)

s_d =  [g.dbin[1] for g in groupby(same1, :dbin; sort=true)]
o_d =  [g.dbin[1] for g in groupby(other1, :dbin; sort=true)]
wo_d = [g.dbin[1] for g in groupby(other2, :dbin; sort=true)]

s_samplesize =  [size(g)[1] for g in groupby(same1, :dbin; sort=true)]
o_samplesize =  [size(g)[1] for g in groupby(other1, :dbin; sort=true)]
wo_samplesize = [size(g)[1] for g in groupby(other2, :dbin; sort=true)]

s_mean =  [mean(g.absdiff) for g in groupby(same1, :dbin; sort=true)]
o_mean =  [mean(g.absdiff) for g in groupby(other1, :dbin; sort=true)]
wo_mean = [mean(g.absdiff) for g in groupby(other2, :dbin; sort=true)]

s_median =  [median(g.absdiff) for g in groupby(same1, :dbin; sort=true)]
o_median =  [median(g.absdiff) for g in groupby(other1, :dbin; sort=true)]
wo_median = [median(g.absdiff) for g in groupby(other2, :dbin; sort=true)]

linewidth=3
markersize=15
l1 = lines!(ax, s_d, s_mean; linewidth, label="same region");
s1 = scatter!(ax, s_d, s_mean; markersize);

l2 = lines!(ax, o_d,o_mean; linewidth, label="other region");
s2 = scatter!(ax, o_d, o_mean; markersize);
axislegend(framevisible=false)

# save(joinpath(DIR, "load_over_dist_1.pdf"), fig)

l3 = lines!(ax, wo_d, wo_mean; linewidth, label="other region (w/o insulator)");
s3 = scatter!(ax, wo_d, wo_mean; markersize);

save(joinpath(DIR, "load_over_dist_2.pdf"), fig)












####
#### Plot relative difference from initial state
####
smask = abs.(same1.staticA) .> 0.2
omask = abs.(other1.staticA) .> 0.2
scatter(same1.absdiff, same1.reldiff; color=same1.staticA)
scatter(same1.absdiff[smask], same1.reldiff[smask]; color=same1.staticA)
boundary = (0.0, 10)
npoints = 2000
width = 1.0
violin(same1.d[smask], same1.absdiff[smask]; side=:left, boundary, npoints, width, show_median=false)
violin!(other1.d[omask], other1.absdiff[omask]; side=:right, boundary, npoints, width, show_median=false)
ylims!(0,0.5)
scatterlines!([Point2f(g.d[begin], median(g.absdiff)) for g in groupby(same1[smask,:], :d)])
scatterlines!([Point2f(g.d[begin], mean(g.absdiff)) for g in groupby(other1[omask,:], :d)])

####
#### Plot steady state difference
####
scatter(same1.d.-0.05, (same1.staticA.-same1.staticB))
scatter!(other1.d.+0.05, (other1.staticA.-other1.staticB))
npoints=1000
violin!(same1.d.-0.1, (same1.staticA.-same1.staticB); side=:left, width=2.0)

scatter!(same2.d.-0.10, (same2.staticA.-same2.staticB))
scatter!(other2.d.+0.05, (other2.staticA.-other2.staticB))

GLMakie.activate!()
####
#### Plot time over distance
####
scatter(same1.d, same1.t)
scatter!(other1.d, other1.t)
scatter!(other2.d, other2.t)

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
