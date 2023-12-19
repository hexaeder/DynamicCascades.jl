# PARAMETERS ###################################################################
k = [4, 10]
# beta = [0.1, 0.5, 0.9]
beta = [0.1, 0.5]

#= frequency bounds [narrow bounds, wide bounds] bounds.
The value in numerator of round(0.1/(2*π) is the angular frequency =#
# freq_bounds = [round(0.1/(2*π), digits=2), round(0.5/(2*π), digits=2)]
freq_bounds = [round(0.1/(2*π), digits=2)]

# failure_mode = [trip_lines, trip_nodes]
# failure_mode = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
failure_modes = [[:dynamic, :dynamic]]

# inertia_values = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 3.0, 4.0, 5.1, 6.1, 7.1, 8.0, 9.0, 10.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
# inertia_values = [0.2, 1.0, 5.0, 10.0, 15.0, 20.0]
inertia_values = [0.2, 5.0]

# constant parameters
N_ensemble_size = 2 # 100
N_nodes = 100

K = 2 # coupling K
gamma = 1 # damping swing equation nodes γ
tau = 1 # time constant τ

alpha = 0.7 # tuning parameter α, :rating = α*K
init_pert = [:line] # initial perturbation set constant to an initial line failure

sigma = 1 # standard deviation σ
mu = 0 # mean μ
################################################################################


#= create hyperparameter: the order is chosen such, that with an increasing number
of finished jobs of the job array the size of the ensemble increases equally,
e.g. when half of all jobs are finished, for each ensemble half of the grids
shall be simulated=#
hyperparam = collect(Iterators.product(beta, inertia_values, freq_bounds, failure_modes, k,
    N_nodes, K, gamma, tau, alpha, init_pert, sigma, mu))[:]

# Repeat hyperparam N_ensemble_size times
hyperparam_ensemble = repeat(hyperparam, N_ensemble_size)
# println(hyperparameter[1])

# TODO check, separate the two Symbols: [:dynamic, :dynamic]
df_hyperparam_ensemble = DataFrame(map(idx -> getindex.(hyperparam_ensemble, idx), eachindex(first(hyperparam_ensemble))),
    [:beta, :inertia_values, :freq_bounds, :failure_modes, :k, :N_nodes, :K, :gamma, :tau, :alpha, :init_pert, :sigma, :mu])

# add "ArrayTaskID"
df_hyperparam_ensemble = hcat(DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble)), df_hyperparam_ensemble)


# generate grids
# loop over all combinations of `k`, `beta` and `inertia_values` in `hyperparameter`
# if steady state for inertia value (mybe check by relaxation)
    # save grids to respective directories
    # save directory file path
    # add colums `graph_seed` and `distr_seed` to DataFrame
    # save graph_seed
    # save distr_seed
    # don't save different inertia_values in grids
# else generate new grid by changing seed


# Save to CSV
CSV.write(joinpath(@__DIR__, "config.csv"), df_hyperparam_ensemble)

##############
network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(1 * 1u"s^2"), graph_seed=124, distr_seed=1230, K=K, α=alpha, γ=(gamma * 1u"s"), τ=(tau * 1u"s"), μ=mu, σ=sigma)
x_static = steadystate(network)


for i in [0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    println("i $i \n")
    network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(i * 1u"s^2"), graph_seed=124, distr_seed=1230, K=K, γ=(gamma * 1u"s"), τ=(tau * 1u"s"), μ=mu, σ=sigma)
    savegraph("/home/vollmich/Desktop/delete_soon/test_graphs/test_graph$i.lg", network)
    # x_static = steadystate(network)
end
################
