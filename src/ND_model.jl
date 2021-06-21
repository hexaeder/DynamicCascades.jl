using Serialization
using NetworkDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq

export nd_model, get_parameters, calculate_apparant_power, initial_fail_cb
export steadystate, simulate, SolutionContainer

import NetworkDynamics.network_dynamics
"""
    network_dynamics(gen_ver, slack_vert, load_vert, flow_edge, network::MetaGraph)

Create a ND-object based on a IEEE_Network struct for the given vertice and edge types.
"""
function network_dynamics(graph::MetaGraph; load_vertex, generator_vertex, slack_vertex, flow_edge)
    vertex_array = Array{ODEVertex}([load_vertex for v in vertices(graph.graph)])
    for i in 1:nv(graph)
        bustype = get_prop(graph, i, :type)
        if bustype == :gen
            vertex_array[i] = generator_vertex
        elseif bustype == :slack
            vertex_array[i] = slack_vertex
        end
    end
    vertex_array = [v for v in vertex_array] # try to narrow type
    # edges = [flow_edge for i in 1:graph.num_of_lines]
    nd = network_dynamics(vertex_array, flow_edge, SimpleGraph(graph))
end

####
#### Models
####

function powerflow!(e, v_s, v_d, K, t)
    e[1] = - K * sin(v_s[1] - v_d[1])
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
    # TODO check inertia units and definition
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

struct SolutionContainer{G,T}
    network::G
    initial_fail::Vector{Int}
    failtime::Float64
    trip_lines::Bool
    sol::T
    load_S::SavedValues{Float64,Vector{Float64}}
    load_P::SavedValues{Float64,Vector{Float64}}
    failures::SavedValues{Float64,Int}
end

function project_theta!(nd, x0)
    θidx = idx_containing(nd, "θ")
    for (i, θ) in enumerate(θidx)
        n = (θ + π) ÷ 2π
        x0[i] = θ - n * 2π
    end
end

function steadystate(network; project=false, verbose=false)
    verbose && println("Find steady state...")
    (nd, p) = nd_model(network)
    x0 = zeros(length(nd.syms));
    x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())

    project && project_theta!(nd, x0)

    θidx = idx_containing(nd, "θ")
    ex = extrema(x_static[θidx])
    if ex[1] < -π/2 || ex[2] > π/2
        @warn "Steadystate: θ ∈ $ex, consider using project=true flag!"
    end

    dx = similar(x_static)
    nd(dx, x_static, p, 0.0)
    @assert isapprox(dx, zeros(length(dx)), atol=1e-11) "No steady state found $(maximum(abs.(dx)))"

    return x_static
end

function simulate(network;
                  verbose=true,
                  x_static=steadystate(network; verbose),
                  initial_fail=Int[],
                  failtime=1.0,
                  tspan=(0., 100.),
                  trip_lines=true,
                  filename=nothing,
                  terminate_steady_state=true)
    (nd, p, overload_cb) = nd_model(network);
    prob = ODEProblem(nd, copy(x_static), tspan, p);

    failures = SavedValues(Float64, Int);
    load_S = SavedValues(Float64, Vector{Float64});
    load_P = SavedValues(Float64, Vector{Float64});

    cbs = CallbackSet(overload_cb(;trip_lines, load_S, load_P, failures, verbose));
    if !isempty(initial_fail)
        cbs = CallbackSet(InitialFailCB(initial_fail, failtime), cbs)
    end
    if terminate_steady_state
        tmin = isempty(initial_fail) ? 0.0 : failtime
        cbs = CallbackSet(cbs, TerminateSteadyStateAfter(tmin))
    end

    verbose && println("Run simulation for trips on lines $(initial_fail)")
    sol = solve(prob, AutoTsit5(Rosenbrock23()), callback=cbs, progress=true);
    container = SolutionContainer(network,
                                  initial_fail, failtime, trip_lines,
                                  sol, load_S, load_P, failures)

    if terminate_steady_state && sol.t[end] < tspan[2]
        verbose && println("Terminated on steady state at $(sol.t[end])")
    else
        @warn "Did not reach steady state! (lines $trip_lines)"
    end

    if filename !== nothing
        return serialize(filename, container)
    else
        return container
    end
end

function nd_model(network::MetaGraph)
    @assert isapprox(sum(describe_nodes(network).P), 0, atol=1e-14) "Power sum not zero!"
    flow_edge = StaticEdge(f! = powerflow!, dim = 1, coupling=:antisymmetric)
    swing_vertex = ODEVertex(f! = swing_equation!, dim = 2, sym=[:θ, :ω])
    load_vertex = ODEVertex(f! = dynamic_load!, dim = 1, sym=[:θ])

    nd = network_dynamics(network;
                          load_vertex=load_vertex,
                          generator_vertex=swing_vertex,
                          slack_vertex=swing_vertex,
                          flow_edge=flow_edge)

    # second we generate the parameter tuple for the simulation
    p = get_parameters(network)

    cb_gen = get_callback_generator(network)

    return (;nd, p, cb_gen)
end

function get_parameters(network::MetaGraph)
    set_admittance!(network)
    edge_p = set_coupling!(network)

    vertex_p = Vector{NTuple{3,Float64}}(undef, nv(network))
    for v in 1:nv(network)
        P = get_prop(network, v, :P)
        inertia = has_prop(network, v, :inertia) ? get_prop(network, v, :inertia) : 0.0
        type = get_prop(network, v, :type)
        τ = if type === :gen || type === :slack
            get_prop(network, v, :damping)
        elseif type === :load
            get_prop(network, v, :timescale)
        end
        vertex_p[v] = (P, inertia, τ)
    end
    return (vertex_p, edge_p)
end

function set_admittance!(network::MetaGraph)
    R::Vector{Float64} = get_prop(network, edges(network), :R)
    X::Vector{Float64} = get_prop(network, edges(network), :X)
    Y = 1 ./ (R .+ im .* X)
    set_prop!(network, edges(network), :_Y, Y)
end

function set_coupling!(network::MetaGraph)
    warn = false
    K = Vector{Float64}(undef, ne(network))
    for (i, e) in enumerate(edges(network))
        v1::Float64 = get_prop(network, e.src, :Vm)
        v2::Float64 = get_prop(network, e.dst, :Vm)
        Y::ComplexF64 = get_prop(network, e, :_Y)
        if !iszero(real(Y))
            warn = true
        end
        K[i] = v1 * v2 * imag(Y)
        set_prop!(network, e, :_K, K[i])
    end
    warn && @warn "There have been nonzero R values, this is undefined!"
    return K
end


function calculate_apparant_power(u, p, t, nd::nd_ODE_Static, network::MetaGraph)
    out = zeros(ne(network))
    calculate_apparant_power!(out, u, p, t, nd, network)
    return out
end

function calculate_apparant_power!(S, u, p, t, nd::nd_ODE_Static, network::MetaGraph; gd=nd(u, p, t, GetGD))
    for (i, e) in enumerate(edges(nd.graph))
        Vs::Float64 = get_prop(network, e.src, :Vm)
        Vd::Float64 = get_prop(network, e.dst, :Vm)
        Y::ComplexF64 = get_prop(network, e, :_Y)

        θs = get_vertex(gd, e.src)[1]
        θd = get_vertex(gd, e.dst)[1]
        K = p[2][i]

        if iszero(K)
            S[i] = 0
        else
            sqr = √(Vs^2 + Vd^2 - 2*abs(Vs)*abs(Vd)*cos(θd-θs))
            S[i] = max(Vs, Vd) * abs(Y) * sqr
        end
    end
    nothing
end

function calculate_active_power(u, p, t, nd::nd_ODE_Static, network::MetaGraph; gd=nd(u, p, t, GetGD))
    P = zeros(ne(network))
    for (i, e) in enumerate(edges(nd.graph))
        P[i] = get_edge(gd, i)[1]
        # TODO: check again with ne code
        # is equivalent to
        # X = - 1.0/rts96.susceptance[i]
        # P[i] = Vs*Vd/X * sin(θs - θd) * sign(K) - get_edge(gd, i)[1]
    end
    return P
end

"""
    get_callback_generator(network::MetaGraph)

This function returns a constructor for the overload callback with two kw arguments
  - `trip_lines=true` : toggle whether the CB should kill lines (set K=0)
  - `load_S=nothing` : provide `SavedValues` type for S values
  - `load_P=nothing` : provide `SavedValues` type for P values
  - `failurs=nothing` : provide `SavedValues` type where to save the failed lines
  - `verbose=true` : toogle verbosity (print on line failures)

This is all a bit hacky. I am creating a SavingCallback for the `load` values. However
the SavingCallback is missing some values right befor the discontinuity. Those values
are injected inside the shutdown affect. In order for this to work the affect needs to
bump the `saveiter` counter of the other callback. Very ugly.
"""
function get_callback_generator(network::MetaGraph)
    function gen_cb(;trip_lines=true, load_S=nothing, load_P=nothing, failures=nothing, verbose=true)
        save_S_fun = let _network=network
            (u, t, integrator) -> calculate_apparant_power(u, integrator.p, t, integrator.f.f, _network)
        end
        save_P_fun = let _network=network
            (u, t, integrator) -> calculate_active_power(u, integrator.p, t, integrator.f.f, _network)
        end

        scb_S = load_S === nothing ? nothing : SavingCallback(save_S_fun, load_S)
        scb_P = load_P === nothing ? nothing : SavingCallback(save_P_fun, load_P)

        condition = let _current_load=zeros(ne(network)), _network=network, _rating=get_prop(network,edges(network),:rating)
            (out, u, t, integrator) -> begin
                # upcrossing through zero triggers condition
                calculate_apparant_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                out .= _current_load .- _rating
                nothing
            end
        end

        affect! = let _failures=failures, _verbose=verbose, _load_S=load_S, _load_P=load_P, _scb_S=scb_S, _scb_P=scb_P
            (integrator, event_idx) -> begin
                _verbose && println("Shutdown line $event_idx at t = $(integrator.t)")

                if _failures !== nothing
                    push!(_failures.t, integrator.t)
                    push!(_failures.saveval, event_idx)
                end
                if _scb_S !== nothing
                    _scb_S.affect!.saveiter += 1
                    push!(_load_S.t, integrator.t)
                    push!(_load_S.saveval, save_S_fun(integrator.u, integrator.t, integrator))
                end
                if _scb_P !== nothing
                    _scb_P.affect!.saveiter += 1
                    push!(_load_P.t, integrator.t)
                    push!(_load_P.saveval, save_P_fun(integrator.u, integrator.t, integrator))
                end
                edge_p = integrator.p[2]
                edge_p[event_idx] = 0.0
                nothing
            end
        end

        if trip_lines
            vcb = VectorContinuousCallback(condition, affect!, ne(network);
                                           save_positions=(true,true),
                                           affect_neg! = nothing, #only trigger on upcrossing -> 0
                                           abstol=10eps(),reltol=0)
        else
            vcb = nothing
        end

        return CallbackSet(vcb, scb_P, scb_S)
    end
    return gen_cb
end

function InitialFailCB(idxs, time)
    affect! = (integrator) -> integrator.p[2][idxs] .= 0.0
    PresetTimeCallback(time, affect!)
end

function TerminateSteadyStateAfter(tmin)
    test = (integ, abs, rel) ->
        integ.t > tmin && DiffEqCallbacks.allDerivPass(integ, abs, rel)
    TerminateSteadyState(1e-8, 1e-6, test)
end
