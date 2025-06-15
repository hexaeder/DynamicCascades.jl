"""
WS test against old ND. This code produces the same results (up to numerical errors), see
/home/brandner/.julia/dev/worktree-mwe_old_ND_maybe_plots/scripts/archive_nb/porting_to_new_ND/tests_old_ND.jl
"""

graph_seed = 11 # df_config.graph_seed[1728] # 11
distr_seed = 11 # df_config.distr_seed[1728] # 11
ensemble_element = 7 #df_config.ensemble_element[1728] #7
filename = "graph_seed=$graph_seed,distr_seed=$distr_seed,k=4,β=0.5,ensemble_element=$ensemble_element"
network = loadgraph(joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_testnetwork.mg"))), MGFormat())

## Translate power injection definition from mwe_old_ND_maybe_plots
# P = P_inj - P_load => P_inj = Pmech = P + P_load
Pmech = get_prop(network, 1:nv(network), :P) .+ get_prop(network, 1:nv(network), :P_load)
Pload = get_prop(network, 1:nv(network), :P_load)
set_prop!(network, 1:nv(network), :Pmech, Pmech)
set_prop!(network, 1:nv(network), :Pload, Pload)

filepath = joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_steady_state_new_ND.csv")))
# x_static = steadystate_new_ND(network; verbose=true, zeroidx=1)  
# # write steady state to .csv
# CSV.write(filepath, Dict(:SteadyState => x_static))

# load stead state 
x_static = DataFrame(CSV.File(filepath)).SteadyState


simulate_new_ND(network;
    gen_model=SwingDynLoadModel,
    x_static=x_static,
    initial_fail=78,
    failtime=0.1,
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    freq_bound = 0.03,
    terminate_steady_state=true,
    solver = Rodas4P(),
    solverargs = solverargs,
    warn=true
    );
#= 
Shutdown line 78 at t=0.1
Line 199 tripped at t=5.953557751899777
Line 194 tripped at t=7.796003330921655
Vertex 98 tripped at t=8.850676880255643
Line 131 tripped at t=16.77656732677714
Line 147 tripped at t=22.48368039743785
Line 146 tripped at t=26.674165748846708
Line 14 tripped at t=28.7945722732567
Line 106 tripped at t=29.79947508886516
Line 135 tripped at t=30.011951768607673
Vertex 59 tripped at t=30.754460769379662
Line 149 tripped at t=32.69434101378367
Vertex 62 tripped at t=33.20882229554317
Line 148 tripped at t=33.29451619720856
Vertex 61 tripped at t=33.469430774560266
Vertex 60 tripped at t=33.760593498538924
Terminated on steady state at 669.4255648511929
=#

simulate_new_ND(network;
    gen_model=SwingDynLoadModel,
    x_static=x_static,
    initial_fail=78,
    failtime=0.1,
    trip_lines = :dynamic,
    trip_nodes = :none,
    freq_bound = 0.03,
    terminate_steady_state=true,
    solver = Rodas4P(),
    solverargs = solverargs,
    warn=true
    );
#= 
Shutdown line 78 at t=0.1
Line 199 tripped at t=5.953557751899777
Line 194 tripped at t=7.796003330921655
Line 140 tripped at t=9.43993900246847
Line 131 tripped at t=17.06531884887773
Line 147 tripped at t=23.112060924038254
Line 106 tripped at t=31.060847765456277
Line 146 tripped at t=34.38172603143674
Line 14 tripped at t=38.32915069634669
Line 135 tripped at t=39.13784790420183
Line 149 tripped at t=41.02532791381812
Line 148 tripped at t=42.28700839938836
Terminated on steady state at 713.0886641873707
=#

simulate_new_ND(network;
    gen_model=SwingDynLoadModel,
    x_static=x_static,
    initial_fail=78,
    failtime=0.1,
    trip_lines = :none,
    trip_nodes = :dynamic,
    freq_bound = 0.01,
    terminate_steady_state=true,
    solver = Rodas4P(),
    solverargs = solverargs,
    warn=true
    );
#= 
Shutdown line 78 at t=0.1
Vertex 29 tripped at t=2.2241521710185306
Vertex 98 tripped at t=2.529786580940513
Terminated on steady state at 551.3118608877102
=#