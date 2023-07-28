using Serialization
using NetworkDynamics
using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using SciMLNLSolve
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

function swing_equation!(dv, v, edges, (_, power, M, _, D), t)
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
function dynamic_load!(dv, v, edges, (_, power, _, τ, _), t)
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

function swing_dynload!(dv, v, edges, (swing_vs_dynload, power, M, τ, D), t)
    if swing_vs_dynload == 1.0 # swing equation node #bool
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

    elseif swing_vs_dynload == 2.0 # dynamical load node
        dv[2] = 0.0 # state not used => shall not be dynamic.
        # v[2] = 0
        dv[1] = power
        for e in edges
            dv[1] += e[1]
        end
        dv[1] = dv[1] / τ
    else
        error()
    end
    nothing
end

# # Implementation of dead band
# function H(ω)
#     if ω < 0.0
#         return 0.0
#     else
#         return 1.0
#     end
# end
#
# function interval(ω, ω_db)
#     abs(H(ω+ω_db) - H(ω_db-ω))
# end
#
# function swing_dynload!(dv, v, edges, (swing_vs_dynload, power, M, τ, D), t)
#     if swing_vs_dynload == 1.0 # swing equation node #bool
#         # dθ = ω
#         dv[1] = v[2]
#         # dω = (P - γ ω + flowsum)/H
#         # dv[2] = power - 0.1 * v[2]
#         ω_db = 0.01
#         dv[2] = power - 0.1 * D * v[2] - 0.9 * D * v[2] * interval(v[2], ω_db)
#         for e in edges
#             dv[2] += e[1]
#         end
#         # ω0 = 2π * 50
#         # dv[2] = dv[2] * 2ω0/H
#         dv[2] = dv[2]/M
#
#     elseif swing_vs_dynload == 2.0 # dynamical load node
#         dv[2] = 0.0 # state not used => shall not be dynamic.
#         # v[2] = 0
#         dv[1] = power
#         for e in edges
#             dv[1] += e[1]
#         end
#         dv[1] = dv[1] / τ
#     else
#         error()
#     end
#     nothing
# end

struct SolutionContainer{G,T}
    network::G
    initial_fail::Vector{Int}
    failtime::Float64
    trip_lines::Symbol
    trip_nodes::Symbol
    trip_load_nodes::Symbol
    sol::T
    frequencies_load_nodes::SavedValues{Float64,Vector{Float64}}
    load_S::SavedValues{Float64,Vector{Float64}}
    load_P::SavedValues{Float64,Vector{Float64}}
    failures::SavedValues{Float64,Int}
    failures_nodes::SavedValues{Float64,Int}
    failures_load_nodes::SavedValues{Float64,Int}
end

terminated(sc::SolutionContainer) = sc.sol.t[end] < sc.sol.prob.tspan[end]

function project_theta!(nd, x0)
    θidx = idx_containing(nd, "θ")
    for (i, θ) in enumerate(θidx)
        n = (θ + π) ÷ 2π
        x0[i] = θ - n * 2π
    end
end

"""
# Arguments
- `zeroidx::Integer=nothing`: If this flag is set, this shifts the phase angle
at all nodes by the phase angle of the node with index `zeroidx`.
"""
function steadystate(network; project=false, verbose=false, tol=1e-7, zeroidx=nothing)
    verbose && println("Find steady state...")
    (nd, p) = nd_model(network)
    x0 = zeros(length(nd.syms));
    # x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())
    x_static = solve(NonlinearProblem(nd, x0, p), NLSolveJL())
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
                  init_pert=:line,
                  P_perturb=0.5,
                  tspan=(0., 100.),
                  trip_lines=:dynamic,
                  trip_nodes=:dynamic,
                  trip_load_nodes=:dynamic,
                  f_min = -2.5,
                  f_max = 1.5,
                  filename=nothing,
                  terminate_steady_state=true,
                  solverargs=(;),
                  warn=true)
    (nd, p, overload_cb) = nd_model(network);
    prob = ODEProblem(nd, copy(x_static), tspan, p);

    # saving for each integration time step the frequencie of each load node
    frequencies_load_nodes = SavedValues(Float64, Vector{Float64});

    #= saving for each integration time step the value of the load (apparent and
    active) power on each line =#
    load_S = SavedValues(Float64, Vector{Float64});
    load_P = SavedValues(Float64, Vector{Float64});

    # saving time and edge that is failing
    failures = SavedValues(Float64, Int);

    # saving time and generator/load node that is failing
    failures_nodes = SavedValues(Float64, Int);
    failures_load_nodes = SavedValues(Float64, Int);

    #= TODO Next line, overload_cb() returns CallbackSet(vccb_lines, vccb_nodes_max, vccb_nodes_min, scb_P, scb_S).
    So it is CallbackSet(CallbackSet(vccb_lines, vccb_nodes_max, vccb_nodes_min, scb_P, scb_S)). Try out if this is
    redundant =#
    cbs = CallbackSet(overload_cb(;trip_lines, trip_nodes, trip_load_nodes, f_min, f_max, frequencies_load_nodes, load_S, load_P, failures, failures_nodes, failures_load_nodes, verbose));
    if !isempty(initial_fail)
        # NOTE add node failures in case of implementing initial node failures
        cbs = CallbackSet(InitialFailCB(network, initial_fail, failtime; init_pert, P_perturb, failures, failures_nodes, verbose), cbs)
    end
    if (trip_lines !== :static || trip_nodes !== :static || trip_load_nodes !== :static) && terminate_steady_state
        min_t = isempty(initial_fail) ? nothing : failtime+eps(failtime)
        # cbs = CallbackSet(cbs, TerminateSteadyState(1e-8, 1e-6; min_t))
        cbs = CallbackSet(cbs, TerminateSelectiveSteadyState(nd; min_t))
    end

    verbose && println("Run simulation for trips on lines/nodes $(initial_fail)")
    # sol = solve(prob, AutoTsit5(Rosenbrock23()); callback=cbs, progress=true, solverargs...); # larger artefact
    # sol = solve(prob, AutoVern9(Rodas5()); callback=cbs, progress=true, solverargs...); # small arefact
    # sol = solve(prob, Rodas5P(); callback=cbs, progress=true, solverargs...); # no artefact
    # sol = solve(prob, Rodas4(); callback=cbs, progress=true, solverargs...);
    sol = solve(prob, Rodas4P(); callback=cbs, progress=true, solverargs...);
    # sol = solve(prob, KenCarp4(); callback=cbs, progress=true, solverargs...);

    container = SolutionContainer(network,
                                  initial_fail, failtime, trip_lines, trip_nodes, trip_load_nodes,
                                  sol, frequencies_load_nodes, load_S, load_P, failures, failures_nodes, failures_load_nodes)

    if terminate_steady_state
        if sol.t[end] < tspan[2]
            verbose && println("Terminated on steady state at $(sol.t[end])")
        else
            warn && @warn "Did not reach steady state! (lines $trip_lines, nodes $trip_nodes, load nodes $trip_load_nodes)"
        end
    end

    if filename !== nothing
        return serialize(filename, container)
    else
        return container
    end
end

# NOTE [2023-04-21 Fr] has to be adapted probably
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

# function nd_model(network::MetaGraph, generator_model::Symbol, load_model::Symbol)
function nd_model(network::MetaGraph)
    @assert isapprox(sum(describe_nodes(network).P), 0, atol=1e-14) "Power sum not zero!"
    flow_edge = StaticEdge(f=powerflow!, dim=1, coupling=:antisymmetric)

    #= model choice hardcoded NOTE: implement as argument in simulate() in case of
    implementing new model =#

    # generator models
    # generator_model = :swing
    generator_model = :swing_dynload

    # load models
    # load_model = :slack
    load_model = :dynload

    vertexmodels = Dict{Symbol, VertexFunction}() # mapping :type to vertex model

    # generator models
    if generator_model === :swing_dynload
        vertexmodels[:gen] = ODEVertex(f=swing_dynload!, dim=2, sym=[:θ, :ω])
        vertexmodels[:syncon] = ODEVertex(f=swing_dynload!, dim=2, sym=[:θ, :ω])
    elseif generator_model === :swing # NOTE implemented only for testing
        vertexmodels[:gen] = ODEVertex(f=swing_equation!, dim=2, sym=[:θ, :ω])
        vertexmodels[:syncon] = ODEVertex(f=swing_equation!, dim=2, sym=[:θ, :ω])
    else
        error("No generator model defined!")
    end

    # load models /slack
    if load_model === :dynload
        vertexmodels[:load] = ODEVertex(f=dynamic_load!, dim=1, sym=[:θ])
    elseif load_model === :slack
        θ_slack = 0.0 # phase angle of slack
        vertexmodels[:load] = StaticVertex(; f=(θ, edges, _, t) -> θ .= θ_slack, dim=1, sym=[:θ])
    else
        error("No load model defined!")
    end

    vertex_array = [vertexmodels[i] for i in get_prop(network, 1:nv(network), :type)]
    #= try to narrow type: if vertex_array used to of e.g Type{Any} as it contained
    two different types but then vertex_array was changed and now only contain elements
    of the same type, then this line results in vertex_array of only a single type.=#
    vertex_array = [v for v in vertex_array]

    nd = network_dynamics(vertex_array, flow_edge, SimpleGraph(network))

    # generate the parameter tuple for the simulation
    p = get_parameters(network)

    cb_gen = get_callback_generator(network,nd)

    return (;nd, p, cb_gen)
end

#= NOTE Vector{NTuple{5,Float64}} generally less performant compared to
Array{Float64}(undef, nv(network), 5) =#
function get_parameters(network::MetaGraph)
    set_inertia!(network)
    edge_p = set_coupling!(network)
    vertex_p = Vector{NTuple{5,Float64}}(undef, nv(network))
    for v in 1:nv(network)
        P = get_prop(network, v, :P)
        τ = ustrip(u"s", get_prop(network, v, :timeconst))
        type = get_prop(network, v, :type)
        if type in (:gen, :syncon)
            M = ustrip(u"s^2", get_prop(network, v, :_M))
            D = ustrip(u"s", get_prop(network, v, :damping))
            swing_vs_dynload = 1.0
        elseif type === :load
            M = NaN # Shall never be used for load. If used => warning (NaN is "infectious")
            D = NaN
            swing_vs_dynload = 2.0
        end
        vertex_p[v] = (swing_vs_dynload, P, M, τ, D)
    end
    return (vertex_p, edge_p)
end

function set_inertia!(network::MetaGraph)
    for v in vertices(network)
        if has_prop(network, v, :_M)
            continue
        elseif has_prop(network, v, :H)
            ω0 = 2π * 50u"1/s" # angular frequency
            H = get_prop(network, v, :H)
            M = upreferred((2H)/ω0)
            # M = upreferred(H/(2ω0)) # old code
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
            sqr = √(Vs^2 + Vd^2 - 2*Vs*Vd*cos(θd-θs))
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
  - `trip_nodes=true` : toggle whether the CB should switch generator to (dynamic) load nodes (and set load to zero)
  - `trip__load_nodes=true` : toggle whether the CB should set the load of (dynamic) load nodes to zero
  - `frequencies_load_nodes=nothing` : provide `SavedValues` type for frequency values of load nodes
  - `load_S=nothing` : provide `SavedValues` type for S values
  - `load_P=nothing` : provide `SavedValues` type for P values
  - `failures=nothing` : provide `SavedValues` type where to save the failed lines
  - `failures_nodes=nothing` : provide `SavedValues` type where to save the failed nodes
  - `failures_load_nodes=nothing` : provide `SavedValues` type where to save the failed laod nodes
  - `verbose=true` : toogle verbosity (print on line and node failures)

This is all a bit hacky. We create a SavingCallback for the `load` and the `frequencies_load_nodes` values. However,
the SavingCallback is missing one value each right before (NB: at?!) the discontinuity. Those values
are injected inside the shutdown affect. In order for this to work the affect needs to
bump the `saveiter` counter of the other (saving) callback. Very ugly.
"""
function get_callback_generator(network::MetaGraph, nd::ODEFunction)
    function gen_cb(;trip_lines=:dynamic, trip_nodes=:dynamic, trip_load_nodes=:dynamic,
                    f_min = -1.0, f_max = 1.0, frequencies_load_nodes=nothing, load_S=nothing, load_P=nothing,
                    failures=nothing, failures_nodes=nothing, failures_load_nodes=nothing, verbose=true)
        save_S_fun = let _network=network
            (u, t, integrator) -> calculate_apparent_power(u, integrator.p, t, integrator.f.f, _network)
        end
        save_P_fun = let _network=network
            (u, t, integrator) -> calculate_active_power(u, integrator.p, t, integrator.f.f, _network)
        end

        scb_S = load_S === nothing ? nothing : SavingCallback(save_S_fun, load_S)
        scb_P = load_P === nothing ? nothing : SavingCallback(save_P_fun, load_P)

        θ_state_idxs= idx_containing(nd, "θ")
        load_node_idxs = findall(x -> x==:load, get_prop(network, 1:nv(network), :type))
        θ_load_node_state_idxs = θ_state_idxs[load_node_idxs]

        # array indicates whether load node is failed: 1: not failed, 0: failed.
        failed_load_nodes = ones(Float64, length(θ_load_node_state_idxs))

        function load_frequencies(u, t, integrator)
            nd = SciMLBase.unwrapped_f(integrator.f.f)
            dx = similar(u)
            nd(dx, u, integrator.p, integrator.t)
            return @views dx[θ_load_node_state_idxs] .* failed_load_nodes
        end

        scb_load_node_frequencies = frequencies_load_nodes === nothing ? nothing : SavingCallback(load_frequencies, frequencies_load_nodes; save_start=true)

        function all_failed_warnings(integrator)
            if length(failures_nodes.t) == count(x -> x==:gen || x==:syncon, get_prop(network, 1:nv(network), :type))
                @warn "All generator nodes have failed, integrator is terminated!"
                terminate!(integrator)
            elseif length(failures_load_nodes.t) == count(x -> x==:load, get_prop(network, 1:nv(network), :type))
                @warn "All load nodes have failed, integrator is terminated!"
                terminate!(integrator)
            elseif length(failures.t) == ne(network)
                @warn "All lines have failed, integrator is terminated!"
                terminate!(integrator)
            end
        end

        if trip_nodes !== :none || trip_load_nodes !== :none
            # frequency bounds
            # Source: https://eepublicdownloads.entsoe.eu/clean-documents/Network%20codes%20documents/NC%20RfG/210412_IGD_Frequency_ranges.pdf
            # f_min = -2.5 # -0.01 # -0.2/(2π) # -1/2π # -0.4
            # f_max = 1.5 # 0.01 # +0.2/(2π) # 1/2π # 0.24
            ω_min = 2π * f_min
            ω_max = 2π * f_max
            # ω_min = -2.5
            # ω_max = 1.5

            ω_state_idxs = idx_containing(nd, "ω") # array: indices of ω-states
            gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs]) # array: indices of vertices that are generators
        end

        ## dynamic line condition
        if trip_lines == :dynamic
            condition = let _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
                (out, u, t, integrator) -> begin
                    # upcrossing through zero triggers condition
                    calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                    out .= _current_load .- _rating
                    nothing
                end
            end

            # function condition(out, u, t, integrator)
            #     current_load = zeros(ne(network))
            #     rating = get_prop(network, edges(network), :rating)
            #     calculate_apparent_power!(current_load, u, integrator.p, t, integrator.f.f, network)
            #     out .= current_load .- rating
            # end

            affect! = let _failures = failures, _verbose = verbose, _load_S = load_S, _load_P = load_P, _scb_S = scb_S, _scb_P = scb_P
                (integrator, event_idx) -> begin
                    _verbose && println("Shutdown line $event_idx at t = $(integrator.t)")
                    #
                    # if _failures !== nothing
                    #     push!(_failures.t, integrator.t)
                    #     push!(_failures.saveval, event_idx)
                    # end
                    # if _scb_S !== nothing
                    #     _scb_S.affect!.saveiter += 1
                    #     push!(_load_S.t, integrator.t)
                    #     push!(_load_S.saveval, save_S_fun(integrator.u, integrator.t, integrator))
                    # end
                    # if _scb_P !== nothing
                    #     _scb_P.affect!.saveiter += 1
                    #     push!(_load_P.t, integrator.t)
                    #     push!(_load_P.saveval, save_P_fun(integrator.u, integrator.t, integrator))
                    # end
                    # edge_p = integrator.p[2]
                    # edge_p[event_idx] = 0.0
                    # # all_failed_warnings(integrator)
                    # auto_dt_reset!(integrator)
                    nothing
                end
            end

            vccb_lines = VectorContinuousCallback(condition, affect!, ne(network); # VectorContinuousCallback (vccb)
                save_positions = (true, true),
                affect_neg! = nothing, #only trigger on upcrossing -> 0
                abstol = 10eps(), reltol = 0)

        else
            vccb_lines = nothing
        end

        ## dynamic node condition
        if trip_nodes == :dynamic
            function node_condition_max(out, u, t, integrator)
                out .= u[ω_state_idxs] .- ω_max
            end

            function node_condition_min(out, u, t, integrator)
                out .= ω_min .- u[ω_state_idxs]
            end

            affect! = let _failures_nodes = failures_nodes, _verbose = verbose
                (integrator, event_idx) -> begin
                    fgen_idx = gen_node_idxs[event_idx] # index of failed generator
                    _verbose && println("Shutdown node $fgen_idx at t = $(integrator.t)")

                    if _failures_nodes !== nothing
                        push!(_failures_nodes.t, integrator.t)
                        push!(_failures_nodes.saveval, fgen_idx)
                    end
                    vertex_p = integrator.p[1]

                    P_adapted = 0.0
                    if has_prop(network, fgen_idx, :P_load) && (trip_load_nodes ==:static || trip_load_nodes ==:none)
                        P_adapted = - ustrip(get_prop(network, fgen_idx, :P_load))
                    end

                    # mutated_tuple = (2.0, 0.0, vertex_p[fgen_idx][3], vertex_p[fgen_idx][4], vertex_p[fgen_idx][5])
                    mutated_tuple = (2.0, P_adapted, vertex_p[fgen_idx][3], vertex_p[fgen_idx][4], vertex_p[fgen_idx][5])
                    vertex_p[fgen_idx] = mutated_tuple
                    integrator.u[ω_state_idxs[event_idx]] = 0.0 # set state to zero to not trigger condition anymore.
                    # all_failed_warnings(integrator)
                    auto_dt_reset!(integrator)
                    nothing
                end
            end

            vccb_nodes_max = VectorContinuousCallback(node_condition_max, affect!, length(ω_state_idxs);
                save_positions = (true, true),
                affect_neg! = nothing, #only trigger on upcrossing -> 0
                abstol = 10eps(), reltol = 0)

            vccb_nodes_min = VectorContinuousCallback(node_condition_min, affect!, length(ω_state_idxs);
                save_positions = (true, true),
                affect_neg! = nothing, #only trigger on upcrossing -> 0
                abstol = 10eps(), reltol = 0)

        else
            vccb_nodes_max = nothing
            vccb_nodes_min = nothing
        end

        ## dynamic load node condition
        if trip_load_nodes == :dynamic
            function load_node_condition_max(out, u, t, integrator)
                out .= load_frequencies(u, t, integrator) .- ω_max
            end

            function load_node_condition_min(out, u, t, integrator)
                out .= ω_min .- load_frequencies(u, t, integrator)
            end

            # affect! = let _failures_load_nodes = failures_load_nodes, _verbose = verbose
            affect! = let _failures_load_nodes = failures_load_nodes, _verbose = verbose, _frequencies_load_nodes = frequencies_load_nodes
                (integrator, event_idx) -> begin
                    fload_node_idx = load_node_idxs[event_idx] # index of failed load node
                    _verbose && println("Shutdown load node $fload_node_idx at t = $(integrator.t)")

                    if _failures_load_nodes !== nothing
                        push!(_failures_load_nodes.t, integrator.t)
                        push!(_failures_load_nodes.saveval, fload_node_idx)
                    end
                    if scb_load_node_frequencies !== nothing
                        scb_load_node_frequencies.affect!.saveiter += 1
                        push!(_frequencies_load_nodes.t, integrator.t)
                        push!(_frequencies_load_nodes.saveval, load_frequencies(integrator.u, integrator.t, integrator))
                    end

                    failed_load_nodes[event_idx] = 0.0
                    vertex_p = integrator.p[1]
                    mutated_tuple = (vertex_p[fload_node_idx][1], 0.0, vertex_p[fload_node_idx][3], vertex_p[fload_node_idx][4], vertex_p[fload_node_idx][5])
                    vertex_p[fload_node_idx] = mutated_tuple
                    # all_failed_warnings(integrator)
                    auto_dt_reset!(integrator)
                    nothing
                end
            end

            vccb_load_nodes_max = VectorContinuousCallback(load_node_condition_max, affect!, length(θ_load_node_state_idxs);
                save_positions = (true, true),
                affect_neg! = nothing,
                abstol = 10eps(), reltol = 0)

            vccb_load_nodes_min = VectorContinuousCallback(load_node_condition_min, affect!, length(θ_load_node_state_idxs);
                save_positions = (true, true),
                # affect_neg! = nothing,
                affect_neg! = affect!,
                abstol = 10eps(), reltol = 0)

        else
            vccb_load_nodes_max = nothing
            vccb_load_nodes_min = nothing
        end

        ## condition for static failures
        if trip_lines == :static || trip_nodes == :static
            ## static line condition
            nd, = nd_model(network)
            condition = getSteadyStateCondition(nd)
            affect! = let _failures = failures, _failures_nodes = failures_nodes, _verbose = verbose, _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
                function (integrator)
                    u = integrator.u
                    t = integrator.t
                    # NOTE add `failures_nodes` when implementing initial node failures
                    # wait til after the failure which just happend
                    if isempty(_failures.t) || t <= _failures.t[end]
                        return
                    end

                    if trip_lines == :static
                        calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                        _current_load .-= _rating
                        triggered_lines = findall(x->x>0, _current_load)
                    else
                        triggered_lines = []
                    end

                    if trip_nodes == :static
                        triggered_nodes = gen_node_idxs[findall(x -> (x > ω_max) || (x < ω_min) , u[ω_state_idxs])]
                    else
                        triggered_nodes = []
                    end

                    if isempty(triggered_lines) && isempty(triggered_nodes)
                        terminate!(integrator)
                    else
                        if trip_nodes == :static && !isempty(triggered_nodes)
                            _verbose && println("Shutdown node ",
                                length(triggered_nodes) == 1 ? only(triggered_nodes) : triggered_nodes,
                                " at t = $(integrator.t)")

                            if !isnothing(failures_nodes)
                                vertex_p = integrator.p[1]
                                for i in triggered_nodes
                                    P_adapted = 0.0
                                    if has_prop(network, i, :P_load) && trip_load_nodes ==:none
                                        P_adapted = - ustrip(get_prop(network, i, :P_load))
                                    end
                                    mutated_tuple = (2.0, P_adapted, vertex_p[i][3], vertex_p[i][4], vertex_p[i][5])
                                    vertex_p[i] = mutated_tuple

                                    integrator.u[ω_state_idxs[findfirst(x -> x == i, gen_node_idxs)]] = 0.0
                                    push!(failures_nodes.t, integrator.t)
                                    push!(failures_nodes.saveval, i)
                                end
                                all_failed_warnings(integrator)
                            end
                        end

                        if trip_lines == :static && !isempty(triggered_lines)
                            _verbose && println("Shutdown line ",
                                length(triggered_lines) == 1 ? only(triggered_lines) : triggered_lines,
                                " at t = $(integrator.t)")

                            integrator.p[2][triggered_lines] .= 0.0

                            if !isnothing(failures)
                                for i in triggered_lines
                                    push!(failures.t, integrator.t)
                                    push!(failures.saveval, i)
                                end
                                all_failed_warnings(integrator)
                            end
                        end
                        auto_dt_reset!(integrator)
                    end
                end
            end
            dcb = DiscreteCallback(condition, affect!)
        else
            dcb = nothing
        end

        # NOTE order of DiscreteCallbacks matters!
        return CallbackSet(vccb_lines, vccb_nodes_max, vccb_nodes_min, vccb_load_nodes_max, vccb_load_nodes_min, dcb, scb_load_node_frequencies, scb_P, scb_S)
    end
    return gen_cb
end

function InitialFailCB(network, idxs, time; init_pert = :line, P_perturb = nothing, failures = nothing, failures_nodes = nothing, verbose = true)
    affect! = function (integrator)

        # init_pert = :node # NOTE TODO Hardcoded => maybe introduce as argument in simulate()
        # init_pert = :power_perturbation
        if init_pert == :line
            # set coupling to zero (corresponds to line failure)
            integrator.p[2][idxs] .= 0.0

            verbose && println("Shutdown line ",
                length(idxs) == 1 ? only(idxs) : idxs,
                " at t = $(integrator.t)")

            if failures !== nothing
                for i in idxs
                    push!(failures.t, integrator.t)
                    push!(failures.saveval, i)
                end
            end
        end

        if init_pert == :node
            nd, = nd_model(network)
            ω_state_idxs = idx_containing(nd, "ω")
            gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
            for i in idxs
                if findfirst(x -> x == i, gen_node_idxs) == nothing
                    error("Index ", i, " is not a generator")
                end
            end

            verbose && println("Shutdown node ",
                length(idxs) == 1 ? only(idxs) : idxs,
                " at t = $(integrator.t)")
            vertex_p = integrator.p[1]
            for i in idxs
                P_adapted = 0.0
                if has_prop(network, i, :P_load)
                    P_adapted = - ustrip(get_prop(network, i, :P_load))
                end

                mutated_tuple = (2.0, P_adapted, vertex_p[i][3], vertex_p[i][4], vertex_p[i][5])
                vertex_p[i] = mutated_tuple

                integrator.u[ω_state_idxs[findfirst(x -> x == i, gen_node_idxs)]] = 0.0
                push!(failures_nodes.t, integrator.t)
                push!(failures_nodes.saveval, i)
            end
        end

        if init_pert == :power_perturbation
            nd, = nd_model(network)
            #= NOTE This allows to perturb both, generator and load nodes.
            The following code block can be removed once we perturb nodes that are not generators.
            =#
            ω_state_idxs = idx_containing(nd, "ω")
            gen_node_idxs = map(s -> parse(Int, String(s)[4:end]), nd.syms[ω_state_idxs])
            for i in idxs
                if findfirst(x -> x == i, gen_node_idxs) == nothing
                    println("Perturbed node ", i, " is not a generator")
                end
            end

            verbose && println("Perturb node(s) ",
                length(idxs) == 1 ? only(idxs) : idxs,
                " at t = $(integrator.t)")
            vertex_p = integrator.p[1]
            for i in idxs
                P_adapted = vertex_p[i][2] + P_perturb
                mutated_tuple = (vertex_p[i][1], P_adapted, vertex_p[i][3], vertex_p[i][4], vertex_p[i][5])
                vertex_p[i] = mutated_tuple
            end
        end
        auto_dt_reset!(integrator)
    end
    PresetTimeCallback(time, affect!)
end

function ChangePCB(fualtp, time) # TODO This function is used nowhere?
end

function getSteadyStateCondition(nd; min_t=nothing)
    idxs = idx_containing(nd, "ω")
    function (u, t, integrator)
        if !isnothing(min_t) && t <= min_t
            return false
        end

        # getting internal cache vector containing elements of `integrator.u`
        testval = first(get_tmp_cache(integrator))
        # writing the current derivative at t into `testval`
        # see https://docs.sciml.ai/DiffEqDocs/dev/basics/integrator/#SciMLBase.get_du!
        DiffEqBase.get_du!(testval, integrator)
        for i in idxs
            if abs(testval[i]) > 1e-6
                return false
            end
        end
        return true
    end
end

# CB terminates integration when in steady state and all faults are dynamically modelled
function TerminateSelectiveSteadyState(nd; min_t = nothing)
    condition = getSteadyStateCondition(nd; min_t)
    affect! = (integrator) -> terminate!(integrator)
    DiscreteCallback(condition, affect!; save_positions = (true, false))
end
