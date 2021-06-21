using LightGraphs
using DataFrames

export max_flow, edge_distance_to, edge_distances

max_flow(c::SolutionContainer; load=:load_S, kwargs...) = max_flow(getproperty(c, load); kwargs...)

function max_flow(load::SavedValues; relative=true)
    staticA = abs.(load.saveval[begin])
    staticB = abs.(load.saveval[end])

    absolutes = abs.(reduce(hcat, load.saveval))
    max_and_idx = map(findmax, eachrow(absolutes))
    maxflow = [val for (val, _) in max_and_idx]

    times = [load.t[idx] for (_, idx) in max_and_idx]

    return (times, staticA, staticB, maxflow)
end

function edge_distances(c::SolutionContainer)
    @assert length(c.initial_fail) == 1 "Multiple initial fails in dataset!"
    idx = c.initial_fail[begin]
    edge_distance_to(c.network, idx) do g, v1, v2
        length(a_star(g, v1, v2))
    end
end

function edge_distance_to(f, g, idx)
    edge = collect(edges(g))[idx]
    v1, v2 = edge.src, edge.dst
    distances = [min(f(g, v1, e.src),
                     f(g, v1, e.dst),
                     f(g, v2, e.src),
                     f(g, v2, e.dst)) + 1
                 for e in edges(g)]
    distances[idx] = 0
    return distances
end

export maxload_dist_time
function maxload_dist_time(c::SolutionContainer)
    network = c.network

    failedge = collect(edges(network))[c.initial_fail[begin]]
    same = get_prop(network, failedge, :region)
    other = same == 1 ? 2 : 1

    sameidx = edgeidx_by_region(network, same)
    i = findfirst(isequal(c.initial_fail[begin]), sameidx)
    deleteat!(sameidx, i)
    otheridx = edgeidx_by_region(network, other)

    t, staticA, staticB, max = max_flow(c)
    d = edge_distances(c)

    region = vcat([:same for i in 1:length(sameidx)], [:other for i in 1:length(otheridx)])
    idx = vcat(sameidx, otheridx)
    d = d[idx]
    t = t[idx]
    staticA = staticA[idx]
    staticB = staticB[idx]
    max = max[idx]
    absdiff = max .- staticA
    reldiff = max./staticA .- 1
    return DataFrame(; idx, region, d, t, staticA, staticB, max, absdiff, reldiff)
end
