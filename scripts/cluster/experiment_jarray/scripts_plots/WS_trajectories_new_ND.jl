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
        vm = SwingDynLoadModel(M=M,D=γ,τ=τ,ωmax=freq_bound,Pmech=Pmech,Pload=Pload)
    elseif type == :load
        vm = DynLoadModel(τ=τ,Pload=Pload)  
    end
    push!(vm_array, vm)
end

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
s0=find_fixpoint(nw)
prob = ODEProblem(nw, uflat(s0), (0, 1), pflat(s0), callback=get_callbacks(nw)); # TODO warum `pflat(s0)` warum extrahiert man Parameter  nicht aus `nw`?
sol = solve(prob, Rodas4P())

