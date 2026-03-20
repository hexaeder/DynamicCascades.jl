"""
Creating dataset for ML clustering of cascade dynamics (AI school).
"""

include(abspath(@__DIR__, "..", "..", "scripts/helpers_jarray.jl"))

# ENHANCEMENTS
# TODO save graph
# TODO save power injections (see job array framework)

tspan=(0., 35.)
number_of_time_steps = 71     # remove duplicates, keep first occurrence
times = round.(collect(range(tspan[1], tspan[2]; length=number_of_time_steps)), digits=4)

# TODO loop over parameters `I`, `initial_fail`, `freq_bound`, `distr_seed` (probably use job array)
inertia_values = [0.2, 5.0, 10.0, 30.0]
freq_bounds = [0.01, 0.03, 0.15]
distr_seeds = [4,5,8]

# NOTE checking if steady state within res_tol exists
# for I in inertia_values, distr_seed in distr_seeds
#     network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, graph_seed=4, μ=0, σ=1, distr_seed=distr_seed,
#                     K=3, α=0.7, M=I*u"s^2", γ=1u"s", τ=1u"s")
#     x_static=steadystate(network; graph=network.graph, zeroidx=1, res_tol=1e-6) 
#     println("I=$I, distr_seed=$distr_seed, steady state found")
# end

df = DataFrame()
for I in inertia_values, f_b in freq_bounds, distr_seed in distr_seeds #
    network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, graph_seed=4, μ=0, σ=1, distr_seed=distr_seed,
                                    K=3, α=0.7, M=I*u"s^2", γ=1u"s", τ=1u"s")
    x_static=steadystate(network; graph=network.graph, zeroidx=1, res_tol=1e-6)
    
    for initial_fail in 1:2:ne(network) # NOTE  
        sol = simulate(network;
            x_static=x_static,
            initial_fail=Int[initial_fail],
            tspan=tspan,
            terminate_steady_state=false,
            freq_bound = f_b);

        features = Symbol[] # array for all features
        values = [] # array for all corresponding values

        # save scalar parameters 
        append!(features, [:I, :f_b, :initial_fail, :distr_seed])
        append!(values, [I, f_b, initial_fail, distr_seed])

        #= # NOTE this is redundant as we can equivalently save `distr_seed` which defines the power injections
        nw = NetworkDynamics.extract_nw(sol)
        p = NWParameter(nw)
        Pload_array = p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pload")), :Pload]
        Pmech_array = p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pmech")), :Pmech]

        # alternative way
        s = NWState(nw)
        Pload_array == s.p.v[1:100, :Pload] # get parameter :sym of vertex idx
        Pmech_array == s.p.v[1:100, :Pmech]
        =#

        # save trajectories
        for t in times
            # frequencies
            append!(features, [Symbol("ω$(i)_t=$t") for i in 1:nv(network)]) # features: frequencies of all nodes at time t
            append!(values, sol(t; idxs = VIndex(1:nv(network), :ω))) # values: frequencies of all nodes at time t
            # apparent power flows
            append!(features, [Symbol("S$(i)_t=$t") for i in 1:ne(network)]) # features: apparent power flow of all nodes at time t
            append!(values, sol(t; idxs = EIndex(1:ne(network), :S))) # values: apparent power flow of all nodes at time t
        end
        append!(df, DataFrame(reshape(values, 1, :), features)) 
    end
end

# Write out parameters and trajectories to .CSV 
savepath = "data/run_251009_01"
cd(abspath(@__DIR__))
isdir(savepath) || mkpath(savepath)
CSV.write(joinpath(savepath, "data_trajectories_fix.csv"), df)

