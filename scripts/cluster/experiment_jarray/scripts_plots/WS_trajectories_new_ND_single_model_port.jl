using NetworkDynamics
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using OrdinaryDiffEq
# using OrdinaryDiffEqTsit5
using MetaGraphs


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
        functional = 1, [description = "If 1, the node is functional (gen + load), if 0, it is only a load"]
        ωmax = Inf, [description = "Maximum rotor frequency, used in callback"]
    end
    @equations begin
        Dt(θ) ~ ifelse(functional > 0 ,
            ω,                  # gen mode
            1/τ * (-Pload + Pel) # load mode
        )
        Dt(ω) ~ ifelse(functional > 0,
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
        # srcV(t) = 1.0, [description = "voltage magnitude at src end"]
        # dstV(t) = 1.0, [description = "voltage magnitude at dst end"]
        S(t), [description = "apparent power flow towards node at dst end"]
    end
    @parameters begin
        K = 3, [description = "coupling constant"]
        # parameters for callback
        active = 1, [description = "If 1, the line is active, if 0, it is tripped"]
        rating = Inf, [description = "active power line rating"]
    end
    @equations begin
        Δθ ~ srcθ - dstθ
        P ~ active*K*sin(Δθ)
        # S ~ active * max(srcV, dstV) * K * √(srcV^2 + dstV^2 - 2*srcV*dstV*cos(Δθ))
        S ~ active*K*√(2 - 2*cos(Δθ))
    end
end

function SwingDynLoadModel(; vidx=nothing, kwargs...)
    model = SwingDynLoad(name = :swing_dyn_load; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)

    # define callback, but only if :ωmax != Inf
    get_default(vm, :ωmax) == Inf && return vm

    cond = ComponentCondition([:ω], [:ωmax]) do u, p, t
        abs(u[:ω]) - p[:ωmax]
    end
    affect = ComponentAffect([:ω], [:functional,:Pmech]) do u, p, ctx
        println("Vertex $(ctx.vidx) tripped at t=$(ctx.integrator.t)")
        u[:ω] = 0.0
        # HACK TODO This is doppeltgemoppelt. Maybe only do `p[:Pmech] = 0`
        p[:functional] = 0
        p[:Pmech] = 0 # NOTE this is necessary for plotting ΔP
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(vm, cb)
    vm
end

# NOTE Maybe redundant. Instead use `SwingDynLoadModel` with `p[:functional] = 0`
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
    affect = ComponentAffect([], [:active]) do u, p, ctx
        println("Line $(ctx.eidx) tripped at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(em, cb)
    em
end

using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
function TerminateSelectiveSteadyState(nw; min_t = nothing)
    # ω_idxs = SII.variable_index(nw, vidxs(1:nv(nw), :ω))
    ω_idxs = SII.variable_index.(nw, VIndex(1:nv(nw), :ω))

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
            if abs(du[ω_idxs[i]]) > 1e-6
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

function simulate_new_ND_single_model_port(exp_data_dir, task_id, initial_fail;
    verbose=true,
    failtime=0.1,
    tspan=(0., 100000.),
    trip_lines=:dynamic,
    trip_nodes=:dynamic,
    terminate_steady_state=true,
    solverargs=(;),
    warn=true)

    ###
    ### read in parameters, graphs and steady states
    ###
    df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

    # load parameters
    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)


    # load network
    network = import_system_wrapper(df_config, task_id)

    # adjust filepaths 
    df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
    df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")

    # HACK 
    #=`import_system_wrapper` generates new MetaGraph using `watts_strogatz` function, which
    generates different graphs for different versions of the package (for the same seeds).
    We need `network` as the power injections are save in the MetaGraph object.
    Thus, we do not use `network.graph` for the graph but load the graph from the file.
    =#
    # load graph
    graph = loadgraph(df_config[task_id,:filepath_graph])

    # # load steady state
    # steady_state_dict  = CSV.File(df_config[task_id,:filepath_steady_state])
    # x_static = steady_state_dict[:SteadyState]

    ###
    ### build NetworkDynamics.jl Network
    ###

    # set failure modes
    trip_nodes == :dynamic ? ωmax = freq_bound*2π : ωmax = Inf
    trip_lines == :dynamic ? rating = α*K : rating = Inf 

    # loop over vertices and assign vertex models & parameter
    vm_array = VertexModel[]
    for i in 1:nv(network)
        # P = P_inj - P_load see `balance_power!`; P_inj = Pmech
        # P_inj = P + P_load
        P = get_prop(network, i, :P)
        Pload = get_prop(network, i, :P_load)
        Pmech = P + Pload

        type = get_prop(network, i, :type)
        if type == :gen
            vm = SwingDynLoadModel(M=M,D=γ,τ=τ,ωmax=ωmax,Pmech=Pmech,Pload=Pload)
        elseif type == :load
            vm = DynLoadModel(τ=τ,Pload=Pload)  
        end
        # set positions for plotting
        set_position!(vm, get_prop(network, i, :pos))
        push!(vm_array, vm)
    end

    # generate `Network` object
    nw = Network(graph, vm_array, Line(K=K,rating=rating); dealias=true)

    # Check if network is power balanced
    nw_state = NWState(nw)
    p = nw_state.p
    @assert isapprox(sum(p.v[1:nv(network), :Pload]), sum(p.v[1:nv(network), :Pmech]))

    # set initial perturbation CB
    init_perturb = PresetTimeComponentCallback(failtime,
        ComponentAffect([], [:active]) do u, p, ctx
            println("Shutdown line $(ctx.eidx) at t=$(ctx.integrator.t)")
            p[:active] = 0
        end
    )
    set_callback!(nw[EIndex(initial_fail)], init_perturb)
    # NOTE
    #= A loop would be sytactically possible but is not the way to go as `TerminateSelectiveSteadyState`
    should not be applied componentwise. =#
    # for i in 1:nv(network)
    #     add_callback!(nw[VIndex(i)], TerminateSelectiveSteadyState(nw))
    # end

    ###
    ### solve ODE problem
    ###
    s0=find_fixpoint(nw)

    # check if initial loads exceed rating (it is also possible to do this check directly in the CB (s. MM HW [2025-03-24 Mo]))
    for i in 1:ne(nw)
        if s0[EIndex(i, :S)] > s0[EIndex(i, :rating)]
            error("At least line $i exceeds the apparent power rating.")
        end
    end

    # prob = ODEProblem(nw, x_static, tspan, pflat(s0), callback=get_callbacks(nw)); # for using saved steady state
    min_t = isempty(initial_fail) ? nothing : failtime+eps(failtime)
    prob = ODEProblem(nw, uflat(s0), tspan, pflat(s0), callback=CallbackSet(get_callbacks(nw), TerminateSelectiveSteadyState(nw;min_t)));
    # prob = ODEProblem(nw, uflat(s0), tspan, pflat(s0), callback=CallbackSet(get_callbacks(nw), TerminateSelectiveSteadyState(nw)));

    sol = solve(prob, Rodas4P(); solverargs...); 

    if terminate_steady_state
        if sol.t[end] < tspan[2]
            verbose && println("Terminated on steady state at $(sol.t[end])")
        else
            warn && @warn "Did not reach steady state!"
        end
    end

    return sol, nw
end


function plot_simulation(simulation_dir, plot_dir, df_config, task_id, initial_fail)
    sol = Serialization.deserialize(joinpath(simulation_dir, "task_id=$task_id,initial_fail=$initial_fail.sol"));
    nw = Serialization.deserialize(joinpath(simulation_dir, "task_id=$task_id,initial_fail=$initial_fail.nw"));

    # calculate indices of failing lines and nodes
    all_failing_nodes_idxs = [i for i in 1:nv(nw) if sol(sol.t[end], idxs=vidxs(i, :functional))[1] == 0]
    all_failing_lines_idxs = [i for i in 1:ne(nw) if sol(sol.t[end], idxs=eidxs(i, :active))[1] == 0]
    node_colors_failing = distinguishable_colors(length(all_failing_nodes_idxs))
    line_colors_failing = distinguishable_colors(length(all_failing_lines_idxs))

    all_nodes_idxs = [i for i in 1:nv(nw)]
    all_lines_idxs = [i for i in 1:ne(nw)]
    node_colors = distinguishable_colors(length(all_nodes_idxs))
    line_colors = distinguishable_colors(length(all_lines_idxs))

    ################################################################################
    ############################ Line and nodes ####################################
    ################################################################################
    fontsize = 35
    titlesize = (fontsize+5)
    linewidth = 3.5
    fig = Figure(size=(3100,1500), fontsize= fontsize)
    # Add a global title in the first row spanning all columns.
    xlim = sol.t[end]/3



    N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)

    # FREQUENCIES ########################################################################
    fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title="Lines & nodes that fail (left column), all lines & nodes (right column):  I=$M,D=$γ,τ=$τ,f_b=$freq_bound,α=$α,K=$K,N=$N,k=$k,β=$β,ensemble element=$ensemble_element", titlealign = :left, titlesize = titlesize)
    for i in all_failing_nodes_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
        color = node_colors_failing[findfirst(x -> x == i, all_failing_nodes_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end

    fig[1,2] = ax = Axis(fig; titlealign = :left, titlesize = titlesize)
    for i in all_nodes_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(i, :f)))))
        color = node_colors[findfirst(x -> x == i, all_nodes_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end
    xlims!(ax, 0, xlim)

    # FLOWS ########################################################################
    fig[2,1] = ax = Axis(fig; xlabel="Time [s]", ylabel="App. power flow [p.u.]")
    for i in all_failing_lines_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
        color = line_colors_failing[findfirst(x -> x == i, all_failing_lines_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end

    fig[2,2] = ax = Axis(fig; xlabel="Time [s]")
    for i in all_lines_idxs
        x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
        color = line_colors[findfirst(x -> x == i, all_lines_idxs)]
        lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
        scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
    end
    xlims!(ax, 0, xlim)

    CairoMakie.save(joinpath(plot_dir, "ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.pdf"),fig)
    CairoMakie.save(joinpath(plot_dir, "ensemble_element=$ensemble_element,I=$M,D=$γ,f_b=$freq_bound,task_id=$task_id,initial_fail=$initial_fail.png"),fig)
    # fig
end