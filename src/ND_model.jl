using NetworkDynamics
using Unitful: @u_str
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using OrdinaryDiffEq
# using OrdinaryDiffEqTsit5
using SciMLNLSolve
using MetaGraphs


export simulate, nd_model_and_CB!, steadystate
export SwingDynLoadModel, SwingDynLoadModel_change_Pmech_only, SwingDynLoadModel_change_to_BH_only


###
### Models
###

@mtkmodel DynLoad begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
    end
    @parameters begin
        τ = 1.0, [description = "dyn load time constant"]
        Pload = 1.0, [description = "Load Power"]
    end
    @equations begin
        Dt(θ) ~ 1/τ * (-Pload + Pel)
        Pinj ~ -Pel
    end
end

@mtkmodel SwingDynLoad begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        ω(t) = 0.0, [description = "Rotor angular frequency"]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
        f(t), [description = "Rotor frequency"]
        ΔP(t), [description = "Power to be compensated at node: Pmech - Pload + Pel"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 1, [description = "Damping"]
        Pmech = 1.0, [description = "Mechanical Power"]
        τ = 1.0, [description = "dyn load time constant"]
        Pload = 1.0, [description = "Load Power"]
        # parameters for callback
        node_swing_stat = 1, [description = """ indicates node type: 
                                                1: swing node
                                                0: full failure (inertia failure and power failure)
                                                2: inertia failure
                                                3: power failure """]
        ωmax = Inf, [description = "Maximum rotor frequency, used in callback"]
    end
    @equations begin
        #= #NOTE depending on `node_swing_stat`, `Pmech` is set to zero in CB. `Pmech` is
        kept in the follwing equations to allow changing `Pmech` to other values than zero 
        without changing these equations.=#
        Dt(θ) ~ ifelse((node_swing_stat == 1) | (node_swing_stat == 3),
            ω,                          # gen mode
            1/τ * (Pmech - Pload + Pel) # load mode, # NOTE Pmech is set to zero in CB
        )
        Dt(ω) ~ ifelse((node_swing_stat == 1) | (node_swing_stat == 3),
            1/M * (Pmech - Pload - D*ω + Pel), # gen mode
            0.0                                # load mode
        )
        Pinj ~ -Pel
        f ~ ω/(2π)
        ΔP ~ Pmech - Pload + Pel
    end
end

@mtkmodel StaticPowerLine begin
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        P(t), [description = "active power flow towards node at dst end", output=true]
        Δθ(t), [description = "voltage angle difference"]
        S(t), [description = "apparent power flow towards node at dst end"]
    end
    @parameters begin
        srcV = 1.0, [description = "voltage magnitude at src end"]
        dstV = 1.0, [description = "voltage magnitude at dst end"]
        # NOTE Its not possible to pass `Y` directly as ModelingToolkit’s @parameters must be real-valued
        K = 3.0, [description = "coupling `K = - srcV * dstV * imag(Y)`"]
        Yabs = 1.0, [description = "admittance magnitude `Yabs = abs(Y)`"]
        # parameters for callback
        line_stat = 1, [description = "If 1, the line is active, if 0, it is tripped"]
        rating = Inf, [description = "apparent power line rating"]
    end
    @equations begin
        Δθ ~ srcθ - dstθ
        P ~ line_stat * K * sin(Δθ) 
        # S ~ line_stat*K*√(2 - 2*cos(Δθ)) # this works only for srcV=dstV=1 (:wattsstrogatz)
        S ~ line_stat * max(srcV, dstV) * Yabs * √(srcV^2 + dstV^2 - 2*srcV*dstV*cos(Δθ))
    end
end

#= TODO Reduce code duplification in `SwingDynLoadModel`, `SwingDynLoadModel_change_Pmech_only` 
and `SwingDynLoadModel_change_to_BH_only` by defining different `affect`-functions.=#
function SwingDynLoadModel(; vidx=nothing, kwargs...)
    model = SwingDynLoad(name = :swing_dyn_load; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)

    # define callback, but only if :ωmax != Inf
    get_default(vm, :ωmax) == Inf && return vm

    cond = ComponentCondition([:ω], [:ωmax]) do u, p, t
        abs(u[:ω]) - p[:ωmax]
    end
    affect = ComponentAffect([:ω], [:node_swing_stat,:Pmech]) do u, p, ctx
        println("Vertex $(ctx.vidx) tripped at t=$(ctx.integrator.t)")
        u[:ω] = 0.0
        p[:node_swing_stat] = 0
        p[:Pmech] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(vm, cb)
    vm
end

# not setting Pmech=0 in CB 
function SwingDynLoadModel_change_to_BH_only(; vidx=nothing, kwargs...)
    model = SwingDynLoad(name = :swing_dyn_load_change_to_BH_only; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)

    # define callback, but only if :ωmax != Inf
    get_default(vm, :ωmax) == Inf && return vm

    cond = ComponentCondition([:ω], [:ωmax]) do u, p, t
        abs(u[:ω]) - p[:ωmax]
    end
    affect = ComponentAffect([:ω], [:node_swing_stat]) do u, p, ctx
        println("Vertex $(ctx.vidx): Node changed to BH (Pmech kept constant) at t=$(ctx.integrator.t)")
        u[:ω] = 0.0
        p[:node_swing_stat] = 2
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(vm, cb)
    vm
end

function SwingDynLoadModel_change_Pmech_only(; vidx=nothing, kwargs...)
    model = SwingDynLoad(name = :swing_dyn_load_change_Pmech_only; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)

    # define callback, but only if :ωmax != Inf
    get_default(vm, :ωmax) == Inf && return vm

    cond = ComponentCondition([:ω], [:ωmax]) do u, p, t
        abs(u[:ω]) - p[:ωmax]
    end
    affect = ComponentAffect([:ω], [:node_swing_stat,:Pmech]) do u, p, ctx
        println("Vertex $(ctx.vidx): Pmech=$(p[:Pmech]) set to Pmech=0 at t=$(ctx.integrator.t)")
        p[:node_swing_stat] = 3
        p[:Pmech] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(vm, cb)
    vm
end

#= #NOTE Maybe redundant to have separate load model.
Instead use `SwingDynLoadModel` with `p[:node_swing_stat] = 0` 
On the other hand it is conceptually more straightforward
to use a different model for a node that can not fail.
=#
function DynLoadModel(; vidx=nothing, kwargs...)
    model = DynLoad(name = :load; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)
    vm
end

function Line(; src=nothing, dst=nothing, kwargs...)
    model = StaticPowerLine(name = :line; kwargs...)
    em = EdgeModel(model, [:srcθ], [:dstθ], AntiSymmetric([:P]))
    !isnothing(src) && set_graphelement!(em, src => dst)

    # define callback, but only if :rating != Inf
    get_default(em, :rating) == Inf && return em

    cond = ComponentCondition([:S], [:rating]) do u, p, t
        u[:S] - p[:rating]
    end
    affect = ComponentAffect([], [:line_stat]) do u, p, ctx
        println("Line $(ctx.eidx) tripped at t=$(ctx.integrator.t)")
        p[:line_stat] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(em, cb)
    em
end

###
### main entry point
###

"""
Choice of ODE solver 
  - `solver = Rodas4P()` see test_new_ND/tests_new_ND.jl for solver choice
  - `solverargs` see test_new_ND/tests_new_ND.jl
  - NOTE With AutoTsit5(Rodas4P()) the CB `TerminateSelectiveSteadyState` does not
    fire. This is probably as `dt` gets too large with higher order stiff solvers (see
    https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems)
"""
function simulate(network;
    verbose=true,
    graph=network.graph,
    gen_model=SwingDynLoadModel,
    x_static=steadystate(network; graph=graph, verbose, zeroidx=1),
    initial_fail=Int[], # multiple failures at once implented, e.g. initial_fail=[1,27]
    failtime=0.1,
    tspan=(0., 100000.),
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    freq_bound = 1.0,
    terminate_steady_state=true,
    solver = Rodas4P(),
    solverargs = (;reltol=1e-8, abstol=1e-6),
    warn=true
    )

    nw = nd_model_and_CB!(network; 
            graph=graph,
            gen_model=gen_model,
            trip_lines=trip_lines,
            trip_nodes=trip_nodes,
            freq_bound=freq_bound)

    ###
    ### checks
    ### 

    #= Check if network is power balanced (this is also checked before creating `nw` when
    importing the MetaGraph network in import_rtsgmlc.jl), s. Schmierzettel 7b.=#
    p = NWParameter(nw)
    @assert isapprox(
        sum(p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pload")), :Pload]),
        sum(p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pmech")), :Pmech]), atol=1e-8)
        
    # check if initial loads exceed rating (this could also be done directly in the CB (s. MM HW [2025-03-24 Mo]))  
    s0 = NWState(nw, x_static, pflat(NWParameter(nw)))
    for i in 1:ne(nw)
        if s0[EIndex(i, :S)] > s0[EIndex(i, :rating)]
            error("At least line $i exceeds the apparent power rating.")
        end
    end

    ###
    ### set initial perturbation CB
    ###
    if !isempty(initial_fail)
        # set initial perturbation CB
        init_perturb = PresetTimeComponentCallback(failtime,
            ComponentAffect([], [:line_stat]) do u, p, ctx
                println("Shutdown line $(ctx.eidx) at t=$(ctx.integrator.t)")
                p[:line_stat] = 0
            end
        )
        # set_callback!(nw[EIndex(initial_fail)], init_perturb)
        for i in initial_fail
            set_callback!(nw[EIndex(i)], init_perturb)
        end
    end

    ###
    ### solve ODE problem
    ###
    min_t = isempty(initial_fail) ? nothing : failtime+eps(failtime)
    prob = ODEProblem(nw, uflat(s0), tspan, pflat(s0), callback=CallbackSet(get_callbacks(nw), TerminateSelectiveSteadyState(nw;min_t)));
    sol = solve(prob, solver; solverargs...);

    if terminate_steady_state
        if sol.t[end] < tspan[2]
            verbose && println("Terminated on steady state at $(sol.t[end])")
        else
            warn && @warn "Did not reach steady state!"
        end
    end

    return sol
end

function simulate(exp_name_date, task_id, initial_fail; kwargs...)

    ###
    ### read in parameters, graphs and steady states
    ###
    df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))

    # load parameters, network (=MetaGraph object with associated parametes) and  graph (=mathematical graph)
    if exp_name_date[1:2] == "WS" 
        _,_,_,_,_,_,_,_,_,_,_,_,freq_bound,trip_lines,trip_nodes,_,_ = get_network_args_stripped(df_config, task_id)

        # load network
        network = import_system_wrapper(df_config, task_id)

        # adjust filepaths 
        df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
        # df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

        #= #HACK `import_system_wrapper` generates `networks`, which is a MetaGraph using `watts_strogatz`.
        `watts_strogatz` generates different graphs for different versions of its package (for the same seeds).
        We need `network` as the power injections are saved in the MetaGraph object.
        Thus, in order to use the same mathematical graph across different package versions for the `watts_strogatz`
        function, we do not use `network.graph` for the graph but load the graph from the file .=#
        # load graph
        graph = loadgraph(df_config[task_id,:filepath_graph])

        # # load steady state # NOTE Be cautious: States are ordered differently in old/new ND.
        # steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
        # x_static = steady_state_dict[:SteadyState]

    elseif exp_name_date[1:3] == "RTS"
        _,_,_,freq_bound,trip_lines,trip_nodes,_,_ = RTS_get_network_args_stripped(df_config, task_id)
        network = RTS_import_system_wrapper(df_config, task_id)
        graph = network.graph
    end

    sol = simulate(network;
        graph=graph,
        initial_fail=initial_fail,
        trip_lines=trip_lines,
        trip_nodes=trip_nodes,
        freq_bound=freq_bound,
        kwargs...
        )

    return sol
end

###
### functions called by `simulate`
###

function nd_model_and_CB!(network::MetaGraph;
        graph=network.graph,
        gen_model=SwingDynLoadModel,
        trip_lines = :none,
        trip_nodes = :none,
        freq_bound = 1.0
        )

    ###
    ### get parameters
    ### 
    set_inertia!(network) # for the WS network it does not add node parameters
    set_admittance!(network) # NOTE this is needed for WS network as well!

    ###
    ### build NetworkDynamics.jl Network
    ###

    # set node failure mode
    trip_nodes == :dynamic ? ωmax = freq_bound*2π : ωmax = Inf # This is a global node property

    # loop over vertices and assign vertex models & parameter
    vm_array = VertexModel[]
    for i in 1:nv(network)
        Pload = ustrip(u"pu", get_prop(network, i, :Pload))
        τ = ustrip(u"s", get_prop(network, i, :timeconst))
        type = get_prop(network, i, :type)
        if type == :gen || type == :syncon
            M = ustrip(u"s^2", get_prop(network, i, :_M))
            γ = ustrip(u"s", get_prop(network, i, :damping))
            Pmech = ustrip(u"pu", get_prop(network, i, :Pmech))
            vm = gen_model(M=M, D=γ,τ=τ,ωmax=ωmax,Pmech=Pmech,Pload=Pload)
        elseif type == :load
            vm = DynLoadModel(τ=τ,Pload=Pload) 
        else
            error("The node type :$type is not defined!")
        end
        # set positions for plotting
        set_position!(vm, get_prop(network, i, :pos))
        push!(vm_array, vm)
    end

    # loop over edges and assign edge parameters
    em_array = EdgeModel[]
    for e in edges(network)
        # set coupling: K = srcV * dstV / X = - srcV * dstV * imag(Y) , with 1/X = - Im(Y), X: reactance
        srcV = ustrip(u"pu", get_prop(network, e.src, :Vm))
        dstV = ustrip(u"pu", get_prop(network, e.dst, :Vm))
        Y = get_prop(network, e, :_Y)
        K = - srcV * dstV * imag(Y)
        set_prop!(network, e, :_K, K)
        Yabs = abs(Y)
        # set line failure mode and rating
        trip_lines == :dynamic ? rating = ustrip(u"pu", get_prop(network, e, :rating)) : rating = Inf 
        em = Line(srcV=srcV, dstV=dstV, K=K, Yabs=Yabs, rating=rating)
        push!(em_array, em)
    end

    # generate `Network` object
    return nw = Network(graph, vm_array, em_array)
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


"""
# Arguments
- `zeroidx::Integer=nothing`: If this flag is set, this shifts the phase angle
at all nodes by the phase angle of the node with index `zeroidx`.
- `relax_init_guess` If set to `true` initial guess (array of zeros) is replaced
by relaxing the network into a steady state, which is then taken as initial guess for 
fixpoint solver.
"""
function steadystate(network;
    verbose=false,
    graph=network.graph,
    zeroidx=nothing,
    res_tol=1e-7,
    relax_init_guess = false,
    problem=SteadyStateProblem,
    solver=NLSolveJL(),
    solverargs=(;)
    )     
    verbose && println("Find steady state...")

    #= The CBs do not affect the steady state of the system.
    This is why the default kwargs of `nd_model_and_CB` are set to `trip_lines = :none` 
    and `trip_nodes = :none` which implies `ωmax = Inf` and `rating = Inf`.=#
    nw = nd_model_and_CB!(network; graph=graph)
    s = NWState(nw)
    x0_init_guess = uflat(s)
    if relax_init_guess == true
        prob = ODEProblem(nw, uflat(s), (0., 10000.), pflat(s), callback=TerminateSelectiveSteadyState(nw; tol=1e-12));
        sol = solve(prob, Rodas4P(); reltol=1e-12, abstol=1e-12);
        x0_init_guess = sol[end]
    end
    p = pflat(s)
    x0 = solve(problem(nw, x0_init_guess, p), solver; solverargs...) # `uflat(s)` creates vector of zeros
    #= `s0` is needed in order to access the symbolic vertex indices that have a θ-state.
    In this model all vertices have a θ-state, however `uflat(s0)` returns the states ordered by 
    the different vertex models/components. So here one cannot modify `x0.u` directly. =#
    s0 = NWState(s, x0.u, p) # `NWState` object with steady state
    x_static = uflat(s0)

    θidx = map(idx -> idx.compidx, vidxs(s, :, "θ"))
    if zeroidx !== nothing
        offset = s0.v[θidx[zeroidx],:θ]
        s0.v[θidx,:θ] .-= offset
        x_static = uflat(s0)
        @assert iszero(s0.v[θidx[zeroidx],:θ])
    end

    #= In order to avoid a n*2π phase shift that leads to the same physical system but allows
    different steady states, the steady state is restricted into [-π,π].=#
    ex = extrema(s0.v[θidx,:θ])
    if ex[1] < -π || ex[2] > π
        error("Steadystate: θ ∈ $ex, consider projecting into [-π,π]!")
    end

    residuum = issteadystate(nw, x_static)
    @assert residuum < res_tol "No steady state found: maximum(abs.(dx))=$residuum"

    return x_static 
end

# checks if LHS of ODE is close to zero.
function issteadystate(nw, x_static)
    dx = similar(x_static)
    nw(dx, x_static, pflat(NWState(nw)), 0.0) # Rate of change (derivative) of each state variable storing it in `dx`.
    return maximum(abs.(dx))
end


using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
function TerminateSelectiveSteadyState(nw; min_t = nothing, tol=1e-6)
    #= `vidxs(nw, :, "ω")` returns symbolic index of vertex IF it has a state ω. This is needed for
    networks with node models that do not have a state \omega.
    For WS (with a single node model that has a ω stat) one can alternatively use 
    `ω_idxs = SII.variable_index.(nw, VIndex(1:nv(nw), :ω))`=#
    ω_idxs = SII.variable_index(nw, vidxs(nw, :, "ω"))

    condition = (u, t, integrator) -> begin
        if !isnothing(min_t) && t <= min_t
            return false
        end

        # getting internal cache vector containing elements of `integrator.u`
        du = first(get_tmp_cache(integrator))
        # writing the current derivative at t into `du`
        # see https://docs.sciml.ai/DiffEqDocs/dev/basics/integrator/#SciMLBase.get_du!
        DiffEqBase.get_du!(du, integrator)
        for i in eachindex(ω_idxs)
            if abs(du[ω_idxs[i]]) > tol
                return false
            end
        end
        return true
    end
    affect! = (integrator) -> terminate!(integrator)
    DiscreteCallback(condition, affect!; save_positions = (true, false))
end
# less performant:
# function TerminateSelectiveSteadyState(; min_t = nothing)
#     condition = (u, t, integrator) -> begin
#         if !isnothing(min_t) && t <= min_t
#             return false
#         end
# 
#         # getting internal cache vector containing elements of `integrator.u`
#         du = first(get_tmp_cache(integrator))
#         # writing the current derivative at t into `du`
#         # see https://docs.sciml.ai/DiffEqDocs/dev/basics/integrator/#SciMLBase.get_du!
#         DiffEqBase.get_du!(du, integrator)
#         nw = NetworkDynamics.extract_nw(integrator)
#         du_state = NWState(nw, du)
#         for i in 1:nv(nw)
#             if abs(du_state.v[i, :ω]) > 1e-6
#                 return false
#             end
#         end
#         return true
#     end
#     affect! = (integrator) -> terminate!(integrator)
#     DiscreteCallback(condition, affect!; save_positions = (true, false))
# end