using NetworkDynamics
using OrdinaryDiffEq
using DiffEqCallbacks

export nd_model, get_parameters, calculate_apparant_power, plot_apparent_power, initial_fail_cb

import NetworkDynamics.network_dynamics
"""
    network_dynamics(gen_ver, slack_vert, load_vert, flow_edge, network::IEEE_Network)

Create a ND-object based on a IEEE_Network struct for the given vertice and edge types.
"""
function network_dynamics(generator_vertex, slack_vertex, load_vertex, flow_edge, network::IEEE_Network) where {G, S, L}
    vertex_array = Array{ODEVertex}([load_vertex for v in vertices(network.graph)])
    for i in 1:network.num_of_buses
        if network.bustype[i] == 2
            vertex_array[i] = generator_vertex
        elseif network.bustype[i] == 3
            vertex_array[i] = slack_vertex
        end
    end
    vertex_array = [v for v in vertex_array] # try to narrow type
    # edges = [flow_edge for i in 1:network.num_of_lines]
    nd = network_dynamics(vertex_array, flow_edge, network.graph)
end

function powerflow!(e, v_s, v_d, K, t)
    e[1] = K * sin(v_s[1] - v_d[1])
    nothing
end

function swing_equation!(dv, v, edges, (power, inertia, τ), t)
    # dθ = ω
    dv[1] = v[2]
    # dω = (P - γ ω + flowsum)/Inertia
    # dv[2] = power - 0.1 * v[2]
    dv[2] = power - τ * v[2]
    for e in edges
        dv[2] += e[1]
    end
    dv[2] = dv[2] / inertia
    nothing
end

#= Following Bergen & Hill (1981) we model loads dynamically:
D⋅dϕᵢ/dt = Pᵢ + ∑ⱼ Kᵢⱼ⋅sin(ϕⱼ-ϕᵢ)
The parameter D defines a time scale for the dynamics. For D → 0 this model approaches
the behavior of the algebraic load model (⩯ instantanious power balance). =#
function dynamic_load!(dv, v, edges, (power, inertia, τ), t)
    # dynamic_load_time = 0.1
    dv[1] = power
    for e in edges
        dv[1] += e[1]
    end
    dv[1] = dv[1] / τ
    nothing
end

function nd_model(rts96::IEEE_Network; gen_τ = 0.1, slack_τ = 0.1, load_τ = 0.1)
    powerflow_edge = StaticEdge(f! = powerflow!, dim = 1, coupling=:antisymmetric)
    swing_vertex = ODEVertex(f! = swing_equation!, dim = 2, sym=[:θ, :ω])
    load_vertex = ODEVertex(f! = dynamic_load!, dim = 1, sym=[:θ])

    nd = network_dynamics(swing_vertex, swing_vertex, load_vertex, powerflow_edge, rts96)

    # second we generate the parameter tuple for the simulation
    p = get_parameters(rts96, gen_τ, slack_τ, load_τ)

    cb = get_overload_callback(rts96)

    return (;nd, p, cb)
end

function get_parameters(rts96::IEEE_Network, gen_τ, slack_τ, load_τ)
    edge_p = copy(rts96.active_power_coupling)
    power = copy(rts96.active_power)
    inertia = copy(rts96.inertia)

    f(::Val{1}) = gen_τ
    f(::Val{2}) = slack_τ
    f(::Val{3}) = load_τ
    τ = [f(Val(bustype)) for bustype in rts96.bustype]

    vertex_p = collect(zip(power, inertia, τ))
    return (vertex_p, edge_p)
end

function calculate_apparant_power(u, p, t, nd::nd_ODE_Static, rts96::IEEE_Network)
    out = zeros(rts96.num_of_lines)
    calculate_apparant_power!(out, u, p, t, nd, rts96)
    return out
end

function calculate_apparant_power!(S, u, p, t, nd::nd_ODE_Static, rts96::IEEE_Network; gd=nd(u, p, t, GetGD))
    for (i, e) in enumerate(edges(nd.graph))
        Vs = rts96.voltage[e.src]
        Vd = rts96.voltage[e.dst]
        X = - 1.0/rts96.susceptance[i]
        θs = get_vertex(gd, e.src)[1]
        θd = get_vertex(gd, e.dst)[1]
        K = p[2][i]

        sqr = √(Vs^2 + Vd^2 - 2*abs(Vs)*abs(Vd)*cos(θd-θs))
        S[i] = abs(max(Vs, Vd)/X * sqr * sign(K))
    end
    nothing
end

function calculate_active_power(u, p, t, nd::nd_ODE_Static, rts96::IEEE_Network; gd=nd(u, p, t, GetGD))
    P = zeros(rts96.num_of_lines)
    for (i, e) in enumerate(edges(nd.graph))
        P[i] = get_edge(gd, i)[1]
        # is equivalent to
        # X = - 1.0/rts96.susceptance[i]
        # P[i] = Vs*Vd/X * sin(θs - θd) * sign(K) - get_edge(gd, i)[1]
    end
    return P
end

"""
    get_overload_callback(rts96::IEEE_Network)

This function returns a constructor for the overload callback with two kw arguments
  - `load=nothing` : provide `SavedValues` type for load values
  - `failurs=nothing` : provide `SavedValues` type where to save the failed lines

This is all a bit hacky. I am creating a SavingCallback for the `load` values. However
the SavingCallback is missing some values right befor the discontinuity. Those values
are injected inside the shutdown affect. In order for this to work the affect needs to
bump the `saveiter` counter of the other callback. Very ugly.
"""
function get_overload_callback(rts96::IEEE_Network)
    function gen_cb(;load_S=nothing, load_P=nothing, failures=nothing, verbose=true)

        save_S_fun = let rts96=rts96
            (u, t, integrator) -> calculate_apparant_power(u, integrator.p, t, integrator.f.f, rts96)
        end
        save_P_fun = let rts96=rts96
            (u, t, integrator) -> calculate_active_power(u, integrator.p, t, integrator.f.f, rts96)
        end


        scb_S = load_S === nothing ? nothing : SavingCallback(save_S_fun, load_S)
        scb_P = load_P === nothing ? nothing : SavingCallback(save_P_fun, load_P)

        condition = let _current_load=zeros(rts96.num_of_lines), rts96=rts96
            (out, u, t, integrator) -> begin
                # upcrossing through zero triggers condition
                calculate_apparant_power!(_current_load, u, integrator.p, t, integrator.f.f, rts96)
                out .= _current_load .- rts96.emergency_rating
                nothing
            end
        end

        affect! = let failures=failures, verbose=verbose, load_S=load_S, load_P=load_P, scb_S=scb_S, scb_P=scb_P
            (integrator, event_idx) -> begin
                verbose && println("Shutdown line $event_idx at t = $(integrator.t)")

                if failures !== nothing
                    push!(failures.t, integrator.t)
                    push!(failures.saveval, event_idx)
                end
                if scb_S !== nothing
                    scb_S.affect!.saveiter += 1
                    push!(load_S.t, integrator.t)
                    push!(load_S.saveval, save_S_fun(integrator.u, integrator.t, integrator))
                end
                if scb_P !== nothing
                    scb_P.affect!.saveiter += 1
                    push!(load_P.t, integrator.t)
                    push!(load_P.saveval, save_P_fun(integrator.u, integrator.t, integrator))
                end
                edge_p = integrator.p[2]
                edge_p[event_idx] = 0.0
                nothing
            end
        end

        vcb = VectorContinuousCallback(condition, affect!, rts96.num_of_lines;
                                       save_positions=(true,true),
                                       affect_neg! = nothing, #only trigger on upcrossing -> 0
                                       abstol=10eps(),reltol=0)
        return CallbackSet(vcb, scb_P, scb_S)
    end
    return gen_cb
end

function initial_fail_cb(idxs, time)
    function affect!(integrator)
        for i in idxs
            integrator.p[2][i] = 0.0
        end
    end

    PresetTimeCallback(time,affect!)
end

function plot_apparent_power(sv, line_idx; limits=nothing, kwargs...)
    data = Vector{Vector{Float64}}(undef, length(line_idx))
    for (i, line_i) in enumerate(line_idx)
        data[i] = map(x->getindex(x, line_i), sv.saveval)
    end
    labels = reshape(["line $i" for i in line_idx], 1, :)
    p = plot(sv.t, data, labels=labels; kwargs...)

    if limits isa IEEE_Network
        trange = [sv.t[begin], sv.t[end]]
        limits = limits.emergency_rating[line_idx]
        labels = reshape(["limit $i" for i in line_idx], 1, :)
        plot!(p, trange, transpose(hcat(limits, limits)), labels=labels)
    end
    return p
end
