include("helpers_jarray.jl")

if ON_YOGA
    using Revise
else # if on PIK-HPC or Pool
    Pkg.instantiate()
    # Pkg.precompile()
end

# PARAMETERS ###################################################################
########### Only for adding new simulations to existing experiment #############
complement_to_existing_exp = false
# existing experiment
existing_exp_name = "WS_k=4_exp02_PIK_HPC_K_=3,N_G=32_20240208_000237.814"
# number of jobs that already have been simulated
N_jobs_total_existing = 864
# number of frquency bounds that are added
N_new_freq_bounds = length([2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40])
################################################################################

# Experiment name
name = "WS_k=4_exp04_vary_I_only_lines_and_nodes_"
long_name = "Only variation of inertia I. As in MA but using the same states for all number_of_task_ids_between_graphs. Line and node failure model." # for providing more details
save_graph_and_filepath = true
solver_name = "Rodas4P()"
steadystate_choice = :rootfind # :relaxation

# Graph params #############
N_nodes = 100
k_vals = [4]
β_vals = [0.5]

# MetaGraph params ###############
K_vals = 3 # coupling K

# NOTE see below "inertia_variation - relate inertia and damping"
inertia_values = [0.2, 0.5, 1.0, 3.0, 5.0, 7.5, 10.0, 20.0, 30.0]
relate_inertia_and_damping = false
γ_eq_sq_I = false
γ_eq_I = false
if relate_inertia_and_damping
    γ_vals = NaN # damping swing equation nodes γ
else
    γ_vals = 1 # damping swing equation nodes γ
end
τ_vals = 1 # time constant τ

# Distribution power injections
σ_vals = 1 # standard deviation σ
μ_vals = 0 # mean μ

N_ensemble_size = 32

# Cascading params ##############
init_pert = [:line] # initial perturbation set constant to an initial line failure
α_vals = 0.7 # tuning parameter α, :rating = α*K
monitored_power_flow = :apparent

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
# freq_bounds = [round(i/(2*π), digits=4) for i in [0.01, 0.1, 0.5, 1.0, 5.0]]
# This is frequency not angular frequency
# freq_bounds = [0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 0.200, 0.210, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.800]
freq_bounds = [0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.220, 0.250, 0.300, 0.800]



# failure_modes = [trip_lines, trip_nodes]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
# failure_modes = [[:dynamic, :dynamic], [:dynamic, :none]]
failure_modes = [[:dynamic, :dynamic]]

exp_name_params = "K_=$K_vals,N_G=$N_ensemble_size"
exp_name = string(name, server_string, exp_name_params)
# exp_name = string(name, "PIK_HPC_", exp_name_params)

################################################################################
preprocess(complement_to_existing_exp, existing_exp_name, exp_name, long_name,
    save_graph_and_filepath, solver_name, steadystate_choice, N_ensemble_size, k_vals, β_vals, N_nodes, 
    inertia_values, K_vals, γ_vals, τ_vals, σ_vals, μ_vals,
    failure_modes, init_pert, freq_bounds, α_vals, monitored_power_flow)