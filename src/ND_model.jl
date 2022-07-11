using Serialization
using NetworkDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using LinearAlgebra
using OrdinaryDiffEq.DiffEqBase

export nd_model, get_parameters, calculate_apparent_power
export steadystate, simulate, SolutionContainer, issteadystate
export terminated

####
#### Models
####

function powerflow!(e, v_s, v_d, K, t)
    e[1] = - K * sin(v_s[1] - v_d[1])
    nothing
end

function swing_equation!(dv, v, edges, (power, M, D), t)
    # dθ = ω
    dv[1] = v[2]
    # dω = (P - γ ω + flowsum)/H
    # dv[2] = power - 0.1 * v[2]
    dv[2] = power - D * v[2]
    for e in edges
        dv[2] += e[1]
    end
    # ω0 = 2π * 50
    # dv[2] = dv[2] * 2ω0/H
    dv[2] = dv[2]/M
    nothing
end

#= Following Bergen & Hill (1981) we model loads dynamically:
D⋅dϕᵢ/dt = Pᵢ + ∑ⱼ Kᵢⱼ⋅sin(ϕⱼ-ϕᵢ)
The parameter D defines a time scale for the dynamics. For D → 0 this model approaches
the behavior of the algebraic load model (⩯ instantanious power balance). =#
function dynamic_load!(dv, v, edges, (power, τ, _), t)
    # dynamic_load_time = 0.1
    dv[1] = power
    for e in edges
        dv[1] += e[1]
    end
    dv[1] = dv[1] / τ
    nothing
end

function algebraic_load!(dv, v, edges, (power, _, _), t)
    dv[1] = power
    for e in edges
        dv[1] += e[1]
    end
    nothing
end

struct SolutionContainer{G,T}
    network::G
    initial_fail::Vector{Int}
    failtime::Float64
    trip_lines::Symbol
    sol::T
    load_S::SavedValues{Float64,Vector{Float64}}
    load_P::SavedValues{Float64,Vector{Float64}}
    failures::SavedValues{Float64,Int}
end

terminated(sc::SolutionContainer) = sc.sol.t[end] < sc.sol.prob.tspan[end]

function project_theta!(nd, x0)
    θidx = idx_containing(nd, "θ")
    for (i, θ) in enumerate(θidx)
        n = (θ + π) ÷ 2π
        x0[i] = θ - n * 2π
    end
end

function steadystate(network; project=false, verbose=false, tol=1e-7, zeroidx=nothing)
    verbose && println("Find steady state...")
    (nd, p) = nd_model(network)
    x0 = zeros(length(nd.syms));
    x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())

    θidx = idx_containing(nd, "θ")
    if zeroidx !== nothing
        offset = x_static[θidx[zeroidx]]
        x_static[θidx] .= x_static[θidx] .- offset
        @assert iszero(x_static[θidx[zeroidx]])
    end

    project && project_theta!(nd, x0)

    ex = extrema(x_static[θidx])
    if ex[1] < -π/2 || ex[2] > π/2
        @warn "Steadystate: θ ∈ $ex, consider using project=true flag!"
    end

    residuum = issteadystate(network, x_static; ndp=(nd, p))
    @assert residuum < tol "No steady state found $residuum"

    return x_static
end

function issteadystate(network, x_static; ndp=nd_model(network))
    (nd, p) = ndp
    dx = similar(x_static)
    nd(dx, x_static, p, 0.0)
    return maximum(abs.(dx))
end

function simulate(network;
                  verbose=true,
                  x_static=steadystate(network; verbose),
                  initial_fail=Int[],
                  failtime=0.1,
                  tspan=(0., 100.),
                  trip_lines=:dynamic,
                  filename=nothing,
                  terminate_steady_state=true,
                  solverargs=(;),
                  warn=true)
    (nd, p, overload_cb) = nd_model(network);
    prob = ODEProblem(nd, copy(x_static), tspan, p);

    failures = SavedValues(Float64, Int);
    load_S = SavedValues(Float64, Vector{Float64});
    load_P = SavedValues(Float64, Vector{Float64});

    cbs = CallbackSet(overload_cb(;trip_lines, load_S, load_P, failures, verbose));
    if !isempty(initial_fail)
        cbs = CallbackSet(InitialFailCB(initial_fail, failtime; failures, verbose), cbs)
    end
    if trip_lines !== :static && terminate_steady_state
        min_t = isempty(initial_fail) ? nothing : failtime+eps(failtime)
        # cbs = CallbackSet(cbs, TerminateSteadyState(1e-8, 1e-6; min_t))
        cbs = CallbackSet(cbs, TerminateSelectiveSteadyState(nd; min_t))
    end

    verbose && println("Run simulation for trips on lines $(initial_fail)")
    sol = solve(prob, AutoTsit5(Rosenbrock23()); callback=cbs, progress=true, solverargs...);
    container = SolutionContainer(network,
                                  initial_fail, failtime, trip_lines,
                                  sol, load_S, load_P, failures)

    if terminate_steady_state
        if sol.t[end] < tspan[2]
            verbose && println("Terminated on steady state at $(sol.t[end])")
        else
            warn && @warn "Did not reach steady state! (lines $trip_lines)"
        end
    end

    if filename !== nothing
        return serialize(filename, container)
    else
        return container
    end
end

function simulate_pdis(network;
                       verbose=true,
                       x_static=steadystate(network; verbose),
                       node=1,
                       disturbance=1u"pu",
                       failtime=0.1,
                       tspan=(0., 100.),
                       terminate_steady_state=true,
                       solverargs=(;),
                       warn=true)
    (nd, p, overload_cb) = nd_model(network);

    # find the fault parameters
    fault = deepcopy(network)
    set_prop!(fault, node, :P, get_prop(fault, node, :P) + disturbance)
    pfault = get_parameters(fault);

    prob = ODEProblem(nd, copy(x_static), tspan, p);

    failures = SavedValues(Float64, Int);
    load_S = SavedValues(Float64, Vector{Float64});
    load_P = SavedValues(Float64, Vector{Float64});

    cbs = CallbackSet(overload_cb(;trip_lines=Int[], load_S, load_P, failures, verbose));
    cbs = CallbackSet(cbs, TerminateSelectiveSteadyState(nd; min_t=failtime+eps(failtime)))
    affect! = (integrator) -> integrator.p = pfault
    cbs = CallbackSet(PresetTimeCallback(failtime, affect!), cbs)

    verbose && println("Run simulation for pertubation $disturbance on lines $(node)")
    sol = solve(prob, AutoTsit5(Rosenbrock23()); callback=cbs, progress=true, solverargs...);
    container = SolutionContainer(network,
                                  [node,], failtime, :none,
                                  sol, load_S, load_P, failures)

    if terminate_steady_state
        if sol.t[end] < tspan[2]
            verbose && println("Terminated on steady state at $(sol.t[end])")
        else
            warn && @warn "Did not reach steady state!"
        end
    end

    return container
end

function nd_model(network::MetaGraph)
    @assert isapprox(sum(describe_nodes(network).P), 0, atol=1e-14) "Power sum not zero!"
    flow_edge = StaticEdge(f=powerflow!, dim=1, coupling=:antisymmetric)

    standardmodels = Dict(
        :gen => :swing,
        :load => :dynload,
        :syncon => :swing,
    )
    for i in 1:nv(network)
        if !has_prop(network, i, :model)
            mod = standardmodels[get_prop(network, i, :type)]
            set_prop!(network, i, :model, mod)
        end
    end

    vertexmodels = Dict(
        :swing => ODEVertex(f=swing_equation!, dim=2, sym=[:θ, :ω]),
        :dynload => ODEVertex(f=dynamic_load!, dim=1, sym=[:θ])
    )

    vertex_array = [vertexmodels[k] for k in get_prop(network, 1:nv(network), :model)]
    vertex_array = [v for v in vertex_array] # try to narrow type

    nd = network_dynamics(vertex_array, flow_edge, SimpleGraph(network))

    # second we generate the parameter tuple for the simulation
    p = get_parameters(network)

    cb_gen = get_callback_generator(network)

    return (;nd, p, cb_gen)
end

function get_parameters(network::MetaGraph)
    set_inertia!(network)
    edge_p = set_coupling!(network)

    vertex_p = Vector{NTuple{3,Float64}}(undef, nv(network))
    for v in 1:nv(network)
        P = get_prop(network, v, :P)
        model = get_prop(network, v, :model)
        if model === :swing
            M = ustrip(u"s^2", get_prop(network, v, :_M))
            D = ustrip(u"s", get_prop(network, v, :damping))
            vertex_p[v] = (P, M, D)
        elseif model === :dynload
            τ = ustrip(u"s", get_prop(network, v, :timeconst))
            vertex_p[v] = (P, τ, 0.0)
        end
    end
    return (vertex_p, edge_p)
end

function set_inertia!(network::MetaGraph)
    for v in vertices(network)
        if has_prop(network, v, :_M)
            continue
        elseif has_prop(network, v, :H)
            ω0 = 2π * 50u"1/s"
            H = get_prop(network, v, :H)
            M = upreferred(H/(2ω0))
            set_prop!(network, v, :_M, M)
        end
    end
end

function set_admittance!(network::MetaGraph)
    R::Vector{Float64} = get_prop(network, edges(network), :R)
    X::Vector{Float64} = get_prop(network, edges(network), :X)
    Y = 1 ./ (R .+ im .* X)
    set_prop!(network, edges(network), :_Y, Y)
end

function set_coupling!(network::MetaGraph)
    if has_prop(network, edges(network), :_K)
        return get_prop(network, edges(network), :_K)
    end

    if !has_prop(network, edges(network), :_Y)
        set_admittance!(network)
    end

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

function set_coupling_sum!(network)
    DynamicCascades.set_admittance!(network)
    DynamicCascades.set_coupling!(network)
    branches = describe_edges(network)
    nodes = describe_nodes(network)

    alledges = edges(network)|>collect

    Ksum = Vector{Float64}(undef, length(nodes.n))
    for n in nodes.n
        eidx = findall(e->n∈(e.src, e.dst), alledges)
        Ksum[n] = sum(branches._K[eidx])
    end
    Ksum
    set_prop!(network, nodes.n, :Ksum, Ksum)
end

function calculate_apparent_power(network::MetaGraph, u)
    nd, p = nd_model(network)
    calculate_apparent_power(u, p, 0.0, nd.f, network)
end

function calculate_apparent_power(u, p, t, nd, network::MetaGraph)
    out = zeros(ne(network))
    calculate_apparent_power!(out, u, p, t, nd, network)
    return out
end

function calculate_apparent_power!(S, u, p, t, nd, network::MetaGraph; gd=nd(u, p, t, GetGD))
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

function calculate_active_power(network::MetaGraph, u)
    nd, p = nd_model(network)
    calculate_active_power(u, p, 0.0, nd, network)
end

function calculate_active_power(u, p, t, nd, network::MetaGraph; gd=nd(u, p, t, GetGD))
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
    function gen_cb(;trip_lines=:dynamic, load_S=nothing, load_P=nothing, failures=nothing, verbose=true)
        save_S_fun = let _network=network
            (u, t, integrator) -> calculate_apparent_power(u, integrator.p, t, integrator.f.f, _network)
        end
        save_P_fun = let _network=network
            (u, t, integrator) -> calculate_active_power(u, integrator.p, t, integrator.f.f, _network)
        end

        scb_S = load_S === nothing ? nothing : SavingCallback(save_S_fun, load_S)
        scb_P = load_P === nothing ? nothing : SavingCallback(save_P_fun, load_P)

        ## static line condition
        if trip_lines == :dynamic
            condition = let _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
                (out, u, t, integrator) -> begin
                    # upcrossing through zero triggers condition
                    calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                    out .= _current_load .- _rating
                    nothing
                end
            end

            affect! = let _failures = failures, _verbose = verbose, _load_S = load_S, _load_P = load_P, _scb_S = scb_S, _scb_P = scb_P
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

            vcb = VectorContinuousCallback(condition, affect!, ne(network);
                save_positions = (true, true),
                affect_neg! = nothing, #only trigger on upcrossing -> 0
                abstol = 10eps(), reltol = 0)

        elseif trip_lines == :static
            ## static line condition
            nd, = nd_model(network)
            condition = getSteadyStateCondition(nd)
            affect! = let _failures = failures, _verbose = verbose, _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
                function (integrator)
                    u = integrator.u
                    t = integrator.t
                    # wait til after the failure which just happend
                    if isempty(_failures.t) || t <= _failures.t[end]
                        return
                    end

                    calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                    _current_load .-= _rating
                    triggered = findall(x->x>0, _current_load)

                    if isempty(triggered)
                        terminate!(integrator)
                    else
                        integrator.p[2][triggered] .= 0.0
                            _verbose && println("Shutdown line ",
                                length(triggered) == 1 ? only(triggered) : triggered,
                                " at t = $(integrator.t)")
                        if !isnothing(failures)
                            for i in triggered
                                push!(failures.t, integrator.t)
                                push!(failures.saveval, i)
                            end
                        end
                    end
                end
            end
            vcb = DiscreteCallback(condition, affect!)

        else
            vcb = nothing
        end

        # static line condition

        return CallbackSet(vcb, scb_P, scb_S)
    end
    return gen_cb
end

function InitialFailCB(idxs, time; failures = nothing, verbose = true)
    affect! = function (integrator)
        # set coupling to zero
        integrator.p[2][idxs] .= 0.0
        verbose && println("Shutdown line ",
            length(idxs) == 1 ? only(idxs) : idxs,
            " at t = $(integrator.t)")
        if !isnothing(failures)
            for i in idxs
                push!(failures.t, integrator.t)
                push!(failures.saveval, i)
            end
        end
    end
    PresetTimeCallback(time, affect!)
end

function ChangePCB(fualtp, time)
end

function getSteadyStateCondition(nd; min_t=nothing)
    idxs = idx_containing(nd, "ω")
    function (u, t, integrator)
        if !isnothing(min_t) && t <= min_t
            return false
        end
        testval = first(get_tmp_cache(integrator))
        DiffEqBase.get_du!(testval, integrator)
        for i in idxs
            if abs(testval[i]) > 1e-6
                return false
            end
        end
        return true
    end
end

function TerminateSelectiveSteadyState(nd; min_t = nothing)
    condition = getSteadyStateCondition(nd; min_t)
    affect! = (integrator) -> terminate!(integrator)
    DiscreteCallback(condition, affect!; save_positions = (true, false))
end
