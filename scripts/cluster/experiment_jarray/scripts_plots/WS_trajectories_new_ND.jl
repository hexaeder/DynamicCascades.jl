"""
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))
include(abspath(@__DIR__, "WS_trajectories_new_ND_single_model_port.jl"))

using DynamicCascades
# using Graphs
using MetaGraphs
# using Unitful
# using Statistics
# using GraphMakie
# using Colors
# using DynamicCascades: PLOT_DIR

# using CairoMakie

###
### read in parameters from .csv
###
initial_fail = 78
task_id = 1720
# task_id_array = [1720, 1723, 1728]

# exp_name_date = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,trip_lines,trip_nodes,init_pert,ensemble_element = get_network_args_stripped(df_config, task_id)
network = import_system_wrapper(df_config, task_id)

###
### build NetworkDynamics.jl Network
###

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
        vm = SwingDynLoadModel(M=M,D=γ,τ=τ,ωmax=freq_bound*2π,Pmech=Pmech,Pload=Pload)
    elseif type == :load
        vm = DynLoadModel(τ=τ,Pload=Pload)  
    end
    push!(vm_array, vm)
end

# BUG
###
### TODO Check initial loads!
###

# BUG
###
### TODO Check consistency of RNG across the different Julia versions.
###

# generate `Network` object
nw = Network(network.graph, vm_array, Line(K=K,rating=α*K); dealias=true)

# Check if network is power balanced
nw_state = NWState(nw)
p = nw_state.p
@assert isapprox(sum(p.v[1:nv(network), :Pload]), sum(p.v[1:nv(network), :Pmech]))

# set initial perturbation CB
init_perturb = PresetTimeComponentCallback(0.1,
    ComponentAffect([], [:active]) do u, p, ctx
        println("Shutdown line $(ctx.eidx) at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
)
set_callback!(nw[EIndex(initial_fail)], init_perturb)

###
### solve ODE problem
###
# steady_state_dict  = CSV.File("/home/brandner/nb_data/HU_Master/2122WS/MA/MA_data/results_NB/WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344/k=4,β=0.5/steady_states_graphs/graph_seed=11,distr_seed=11,k=4,β=0.5,ensemble_element=7.csv")
# steady_state = steady_state_dict[:SteadyState]
s0=find_fixpoint(nw)
# prob = ODEProblem(nw, steady_state, (0, 1), pflat(s0), callback=get_callbacks(nw));
prob = ODEProblem(nw, uflat(s0), (0, 2), pflat(s0), callback=get_callbacks(nw)); #get callbacks checkt components und baut vector callback.
sol = solve(prob, Rodas4P());


####################################################################################################

###
### Fixpoint calculation for old ND and new ND
###

## old ND
zeroidx=1
(nd, p) = nd_model(network)
x0 = zeros(length(nd.syms));
x_static = solve(NonlinearProblem(nd, x0, p), NLSolveJL())

θidx = idx_containing(nd, "θ")
offset = x_static[θidx[zeroidx]]
x_static[θidx] .= x_static[θidx] .- offset
@assert iszero(x_static[θidx[zeroidx]])

## new ND
zeroidx=1
nw_state = NWState(nw)

x0 = zeros(dim(nw))
#= This does not work anymore in new DE.jl. Possible reason: The system defined by `nw` is an
ODE-System and not an nonlinear problem (which would be nonlinear algebraic equations). In the 
more recent version of DE.jl this throws an error to avoid an inconsisten use of `NonlinearProblem`,
i.e. trying to solve a nonlinear problem that in fact is not a nonlinear problem. Instead, the
natural way is to define a steady state problem.=#
x_static = solve(NonlinearProblem(nw, x0, pflat(nw_state)), NLSolveJL())

# This works
solve(NonlinearProblem((du,u,p)->nw(du,u,p,0.0), x0, pflat(nw_state)), NLSolveJL())

#= This throws an error, because `_solve_fixpoint` which is invoced by `find_fixpoint` expects
the type `AbstractNonlinearSolveAlgorithm`, which doesn't fit `typeof(NLSolveJL())`. =#
find_fixpoint(nw, alg=NLSolveJL()) 

# This works
solve(SteadyStateProblem(nw, x0, pflat(nw_state)), NLSolveJL())

## Getting rid off phase-shift
# (first try)
nw[EIndex(initial_fail)]
uflat_s0 = uflat(s0)
offset = uflat_s0[1]
uflat_s0[1:2:end] .-= offset
uflat_s0

# (Notes HW)
s0.v[1:100,:ω]  
s0.v[:,:θ]
s0.v[1:100,:θ]
s0.v[1:100,:θ] .- s0.v[1,:θ]


# compare steady state for old and new ND.jl
isapprox(steady_state, uflat_s0, atol=1e-3)