using Dates
using Printf
using Graphs
using OrdinaryDiffEq
using Makie
using Makie.GeometryBasics

export timestamp, pstring
export edgeidx_by_region, vertexidx_by_region

timestamp() = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

function pstring(n=2; kwargs...)
    ret = ""
    for (k,v) in kwargs
        ret *= String(k)
        if v isa AbstractFloat
            ret *= @sprintf("%.3f", v)
            repr(round(5.214, digits=8))
        elseif v isa AbstractString
            ret *= v
        else
            ret *= repr(v)
        end
        ret *= "_"
    end
    return replace(ret, " "=>"")[1:end-1]
end

function edgeidx_by_region(g, region)
    regions = get_prop(g, edges(g), :region)
    mask = regions .== region
    return (1:ne(g))[mask]
end

function vertexidx_by_region(g, region)
    regions = get_prop(g, 1:nv(g), :region)
    mask = regions .== region
    return (1:nv(g))[mask]
end

import MetaGraphs: set_prop!, get_prop

KEY_ITER = Union{AbstractUnitRange,Vector,AbstractEdgeIter}
"""
    set_prop!(g, keys::Iterable, prop::Symbol, vals::Iterable)

Set same property `prop` with different values `vals` for differet identifiers `keys`.
"""
function set_prop!(g, keys::KEY_ITER, prop::Symbol, vals::Vector)
    length(keys) == length(vals) || throw(ArgumentError("keys and vals needs to be of same length!"))
    for (k, val) in zip(keys, vals)
        if !ismissing(val)
            set_prop!(g, k, prop, val)
        end
    end
end

set_prop!(g, keys::KEY_ITER, p::Symbol, val) = set_prop!(g, keys, p, [val for v in keys])

"""
    get_prop(g, keys::Iterable, prop::Symbol)

Get same property `prop` with different values `vals` for differet keys.
"""
function get_prop(g, keys::KEY_ITER, prop::Symbol)
    [has_prop(g, k, prop) ? get_prop(g, k, prop) : missing for k in keys]
end

import DiffEqCallbacks: SavedValues
function (sv::SavedValues)(t)
    if t < sv.t[begin] || t > sv.t[end]
        throw(ArgumentError("t out of bounds!"))
    elseif t == sv.t[end]
        return sv.saveval[end]
    else
        i2 = findfirst(x -> x > t, sv.t)
        i1 = i2 - 1
        v1, v2 = sv.saveval[i1], sv.saveval[i2]
        t1, t2 = sv.t[i1], sv.t[i2]
        return @. v1 + (t-t1)/(t2-t1) * (v2-v1)
    end
end

export seriesforidx
function seriesforidx(sv::Union{SavedValues,ODESolution}, idx; cuttail=true, T=Float32)
    x = Vector{T}(sv.t)
    ytable = sv isa SavedValues ? sv.saveval : sv.u
    y = [T(l[idx]) for l in ytable]

    cuttail && remove_zero_tail!(x, y)

    return (x, y)
end

function remove_zero_tail!(x, y)
    @assert length(x) == length(y)
    tail = lastindex(y) + 1
    while y[tail - 1] ≈ 0
        tail -= 1
    end
    if tail <= length(y)
        idxs = tail:lastindex(y)
        deleteat!(x, idxs)
        deleteat!(y, idxs)
    end
    return (x, y)
end

export thetashape
function thetashape(θ; loc=(0,0), r=0.5, b=0.05, pin=20, pout=100)
    points = Vector{Point2f0}(undef, pin+pout)
    # inner half circle
    for (i, ϕ) in enumerate(range(π/2 + θ, 3/2*π + θ, length=pin))
        points[i] = (b*sin(ϕ), b*cos(ϕ)).+loc
    end
    # outer half circle
    α = atan(b, r)
    for (i, ϕ) in enumerate(range(2π - α + θ, α + θ, length=pout))
        points[i+pin] = (r*sin(ϕ), r*cos(ϕ)).+loc
    end
    points
end

export edgeidx
function edgeidx(g::AbstractGraph, s, d)
    s, d = min(s, d), max(s, d)
    for (i, e) in enumerate(edges(g))
        if e.src==s && e.dst==d
            return i
        end
    end
    nothing
end

export getedge
getedge(g::AbstractGraph, i) = collect(edges(g))[i]

