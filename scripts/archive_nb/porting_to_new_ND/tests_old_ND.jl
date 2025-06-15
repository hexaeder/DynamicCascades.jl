"""
See scripts/archive_nb/porting_to_new_ND/tests_new_ND.jl
"""
using Revise

###################################################################################
###################################### RTS ########################################
###################################################################################

# Solver: Rodas4P()
solverargs = (;reltol=1e-6, abstol=1e-6)
include("/home/brandner/.julia/dev/worktree-mwe_old_ND_maybe_plots/scripts/cluster/experiment_jarray/helpers_jarray.jl")

###
### stiff parameter setting 1
###
damping = 0.1u"s"
scale_inertia = 0.2 
tconst = 0.01u"s"
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

## lines + nodes
sol = simulate(network;
    x_static = steadystate(network; tol=1e-5, zeroidx=1),
    initial_fail = Int[27],
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    f_min = -0.42,
    f_max = 0.42,
    monitored_power_flow = :apparent,
    solverargs = solverargs,
    verbose = true);
#= 
Run simulation for trips on lines/nodes [27]
Shutdown line 27 at t = 0.1
Shutdown node 21 at t = 0.1064657989094697
Shutdown node 15 at t = 0.11191490706886047
Shutdown node 16 at t = 0.11390874793453909
Shutdown node 18 at t = 0.11957229210780022
Shutdown node 14 at t = 0.12318726784023706
Shutdown node 22 at t = 0.13849432307516196
Shutdown node 23 at t = 0.14103642497346353
Shutdown line 11 at t = 0.1427245931039492
Shutdown node 13 at t = 0.1436414225512485
Shutdown node 71 at t = 0.1442255556992654
Shutdown node 61 at t = 0.1476905989223052
Shutdown node 64 at t = 0.14808923228060047
Shutdown node 42 at t = 0.14847133089922496
Shutdown node 40 at t = 0.1494796606156778
Shutdown line 3 at t = 0.1508492738911071
Shutdown line 24 at t = 0.15160638495271256
Shutdown node 62 at t = 0.1520975015582945
Shutdown node 45 at t = 0.1522289862684344
Shutdown line 38 at t = 0.1530805274931234
Shutdown line 4 at t = 0.15324584210471665
Shutdown line 2 at t = 0.15338340130447106
Shutdown line 5 at t = 0.1535463699634967
Shutdown line 101 at t = 0.15441608337221677
Shutdown node 7 at t = 0.15603605218255256
Shutdown node 39 at t = 0.15614743483639693
Shutdown line 96 at t = 0.1586467283983145
Shutdown line 76 at t = 0.16011976915675907
Shutdown line 84 at t = 0.16020022591385316
Shutdown line 79 at t = 0.16055938707094558
Shutdown line 77 at t = 0.16157340543242127
Shutdown node 38 at t = 0.1616400808283702
Shutdown line 78 at t = 0.1617100910555968
Shutdown node 46 at t = 0.1677363846575386
Shutdown line 71 at t = 0.17000940125275615
Shutdown line 56 at t = 0.17524159144508414
Shutdown line 40 at t = 0.17759527296181948
Shutdown line 44 at t = 0.17808608895770822
Shutdown line 57 at t = 0.18205144106398424
Shutdown node 31 at t = 0.18791081925835212
Shutdown node 26 at t = 0.19295539197104394
Shutdown line 39 at t = 0.19316749159055124
Shutdown node 25 at t = 0.2035395592299696
Shutdown line 55 at t = 0.2171580188559864
Shutdown line 53 at t = 0.21741009385387108
Shutdown node 37 at t = 0.22514819631521668
Shutdown node 47 at t = 0.2295605473478557
Shutdown node 1 at t = 0.27440179024327394
Shutdown node 69 at t = 0.28454079267441684
Shutdown node 66 at t = 0.3063368016804435
Shutdown node 63 at t = 0.31182415570087163
Shutdown line 107 at t = 0.31983575243699464
Shutdown line 103 at t = 0.32222463369401916
Shutdown line 75 at t = 0.3263004177618686
Shutdown node 70 at t = 0.3642990007570975
Shutdown node 2 at t = 0.5181331631478133
Terminated on steady state at 5.249213593722638
=#

## lines
sol = simulate(network;
    x_static = steadystate(network; tol=1e-5, zeroidx=1),
    initial_fail = Int[27],
    trip_lines = :dynamic,
    trip_nodes = :none,
    monitored_power_flow = :apparent,
    solverargs = solverargs,
    verbose = true);

#= 
Run simulation for trips on lines/nodes [27]
Shutdown line 27 at t = 0.1
Shutdown line 29 at t = 0.2279766259580031
Shutdown line 11 at t = 0.2896307813144398
Shutdown line 12 at t = 0.3532388999974855
Shutdown line 37 at t = 0.36769308788456123
Terminated on steady state at 5.964961825948961
=#

## nodes
sol = simulate(network;
    x_static = steadystate(network; tol=1e-5, zeroidx=1),
    initial_fail = Int[27],
    trip_lines = :none,
    trip_nodes = :dynamic,
    f_min = -0.42,
    f_max = 0.42,
    monitored_power_flow = :apparent,
    solverargs = solverargs,
    verbose = true);

#= 
Run simulation for trips on lines/nodes [27]
Shutdown line 27 at t = 0.1
Shutdown node 21 at t = 0.1064657989094697
Shutdown node 15 at t = 0.11191490706886047
Shutdown node 16 at t = 0.11390874793453909
Shutdown node 18 at t = 0.11957229210780022
Shutdown node 14 at t = 0.12318726784023706
Shutdown node 22 at t = 0.13849432307516196
Shutdown node 23 at t = 0.14103642497346353
Shutdown node 13 at t = 0.14365432197165914
Shutdown node 71 at t = 0.1442257141044533
Shutdown node 61 at t = 0.14769674882078324
Shutdown node 64 at t = 0.14809765125649346
Shutdown node 42 at t = 0.14864671775659172
Shutdown node 40 at t = 0.14973346422764275
Shutdown node 62 at t = 0.15212301587877383
Shutdown node 45 at t = 0.1525037057862011
Shutdown node 39 at t = 0.15363708952037416
Shutdown node 7 at t = 0.15453665853983173
Shutdown node 1 at t = 0.1567941333192985
Shutdown node 38 at t = 0.1571553924632905
Shutdown node 2 at t = 0.15716716189425223
Shutdown node 66 at t = 0.15888377384518315
Shutdown node 63 at t = 0.16116246615974916
Shutdown node 46 at t = 0.16131709574391795
Shutdown node 69 at t = 0.16220411030160345
Shutdown node 55 at t = 0.16460216593257415
Shutdown node 47 at t = 0.16473810849527848
Shutdown node 37 at t = 0.16620085035790255
Shutdown node 49 at t = 0.1688395533878113
Shutdown node 70 at t = 0.16889511551961364
Shutdown node 50 at t = 0.16921569063051947
Shutdown node 25 at t = 0.1714509337816712
Shutdown node 31 at t = 0.17163986185895497
Shutdown node 26 at t = 0.1726221332103448
Terminated on steady state at 0.17265207254470774
=#


###
### stiff parameter setting 2
###
damping = 0.1u"s"
scale_inertia = 20.0 
tconst = 0.01u"s"
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

sol = simulate(network;
    x_static = steadystate(network; tol=1e-5, zeroidx=1),
    initial_fail = Int[27],
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    f_min = -0.42,
    f_max = 0.42,
    monitored_power_flow = :apparent,
    solverargs = solverargs,
    verbose = true);
#=
Run simulation for trips on lines/nodes [27]
Shutdown line 27 at t = 0.1
Shutdown line 29 at t = 0.47378384913066857
Shutdown line 11 at t = 0.8996347901175625
Shutdown line 12 at t = 1.4642997853856996
Shutdown line 37 at t = 1.669382919848303
Shutdown node 18 at t = 1.9730690151909214
Shutdown node 7 at t = 2.0029817159070396
Shutdown line 35 at t = 2.0991911486582215
Shutdown line 19 at t = 2.514104264954847
Shutdown node 21 at t = 2.644999699411969
Shutdown line 18 at t = 2.800859527098185
Shutdown line 16 at t = 2.8008618754927657
Shutdown line 3 at t = 2.8037322043144552
Shutdown line 5 at t = 2.8064636637445375
Shutdown line 4 at t = 2.806725224456402
Shutdown line 2 at t = 2.8104009231828453
Shutdown line 20 at t = 2.815448764494649
Shutdown line 6 at t = 2.819087197430302
Shutdown node 14 at t = 2.865035493456324
Shutdown node 16 at t = 2.890701549440099
Shutdown node 15 at t = 2.9426960925432537
Shutdown line 65 at t = 4.37681117212307
Shutdown line 62 at t = 4.47678671292867
Shutdown line 63 at t = 4.682808880279084
Shutdown line 24 at t = 5.205136178502981
Shutdown line 41 at t = 5.238976671485123
Shutdown line 43 at t = 5.32635269226053
Shutdown line 49 at t = 5.328410262852333
Shutdown line 42 at t = 5.332399200106098
Shutdown line 44 at t = 5.332716635142535
Shutdown node 45 at t = 5.70671492351767
Shutdown line 39 at t = 5.713448313551811
Shutdown node 23 at t = 5.739609883515779
Shutdown node 40 at t = 5.742715861115942
Shutdown line 40 at t = 5.897387425357938
Shutdown node 42 at t = 6.6415417958591485
Shutdown line 71 at t = 6.941665964642887
Shutdown line 73 at t = 6.980434447563574
Shutdown node 38 at t = 6.998914523810241
Shutdown line 56 at t = 7.003269598729543
Shutdown line 55 at t = 7.012623988330783
Shutdown line 53 at t = 7.012757326597221
Shutdown node 39 at t = 7.768823859298056
Shutdown node 37 at t = 8.505381986012022
Shutdown node 47 at t = 8.851791361810104
Shutdown node 26 at t = 10.908200569490885
Shutdown line 100 at t = 12.151407017225834
Shutdown line 97 at t = 12.269223447099431
Shutdown node 2 at t = 12.36824308846532
Shutdown line 98 at t = 12.51476570005224
Shutdown node 66 at t = 14.300504413815613
Shutdown node 62 at t = 14.425571607458362
Shutdown node 64 at t = 14.551745557279782
Shutdown node 69 at t = 14.6673032899722
Shutdown node 71 at t = 14.698042225516215
Shutdown node 61 at t = 14.715345869303496
Shutdown line 84 at t = 14.716537846302773
Shutdown line 76 at t = 14.719497167945418
Shutdown line 77 at t = 14.720889651758881
Shutdown line 78 at t = 14.72119246680668
Shutdown line 75 at t = 14.722892542148376
Shutdown line 79 at t = 14.72717020957541
Shutdown node 63 at t = 14.93136227780051
Shutdown node 25 at t = 30.50522873483372
Shutdown node 13 at t = 36.25818482228428
Shutdown node 46 at t = 36.44724097756454
Shutdown node 1 at t = 37.03828151533727
Shutdown node 22 at t = 59.16674739805526
Terminated on steady state at 255.6594950858193
=#



###################################################################################
####################################### WS ########################################
###################################################################################
"""
see for comparison to new ND: test_new_ND/WS_test_against_old_ND.jl
"""

# save network 
# exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
# df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
# task_id = 1728
# network = import_system_wrapper(df_config, task_id)

graph_seed = 11 # df_config.graph_seed[1728] # 11
distr_seed = 11 # df_config.distr_seed[1728] # 11
ensemble_element = 7 #df_config.ensemble_element[1728] #7
filename = "graph_seed=$graph_seed,distr_seed=$distr_seed,k=4,β=0.5,ensemble_element=$ensemble_element"

# savegraph(joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_testnetwork.mg"))), network)

# # save graph as RNG of `wattsstrogatz` changed in the past
# savegraph(joinpath("/home/brandner/.julia/dev/DynamicCascades/test_new_ND/WS_test_data", string(filename,"_graph_testnetwork_old")), network.graph)
# savegraph(joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_graph_testnetwork_old"))), network.graph)

# # write power injections to .csv (saved in main branch )
# filepath1 = joinpath("/home/brandner/.julia/dev/DynamicCascades/test_new_ND/WS_test_data", string(filename,"_power_injections_testnetwork_old.csv"))
# filepath2 = joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_power_injections_testnetwork_old.csv")))
# P = get_prop(network, 1:nv(network), :P)
# P_load = get_prop(network, 1:nv(network), :P_load)
# P_inj = get_prop(network, 1:nv(network), :P_inj)
# CSV.write(filepath1, Dict(:P => P, :P_load => P_load, :P_inj => P_inj))
# CSV.write(filepath2, Dict(:P => P, :P_load => P_load, :P_inj => P_inj))

network = loadgraph(joinpath(abspath(@__DIR__,"WS_test_data", string(filename,"_testnetwork.mg"))), MGFormat())
x_static_old_ND = DataFrame(CSV.File(abspath(@__DIR__,"WS_test_data", string(filename,"_steady_state_old_ND.csv")))).SteadyState

simulate(network;
    verbose=true,
    x_static=x_static_old_ND,
    initial_fail=[78],
    failtime=0.1,
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    f_min = -0.03,
    f_max = 0.03,
    terminate_steady_state=true,
    solverargs = solverargs,
    warn=true
    );
#= 
Run simulation for trips on lines/nodes [78]
Shutdown line 78 at t = 0.1
Shutdown line 199 at t = 5.953557774825518
Shutdown line 194 at t = 7.796003380028169
Shutdown node 98 at t = 8.850676925822473
Shutdown line 131 at t = 16.776567341299007
Shutdown line 147 at t = 22.48368041000715
Shutdown line 146 at t = 26.674165803049135
Shutdown line 14 at t = 28.794572334749347
Shutdown line 106 at t = 29.799475130747435
Shutdown line 135 at t = 30.011951819615504
Shutdown node 59 at t = 30.754460827439868
Shutdown line 149 at t = 32.69434107062934
Shutdown node 62 at t = 33.208822361685144
Shutdown line 148 at t = 33.29451626111776
Shutdown node 61 at t = 33.46943082098263
Shutdown node 60 at t = 33.76059356376504
Terminated on steady state at 669.4255749465186
=#    

# lines
simulate(network;
    verbose=true,
    x_static=x_static_old_ND,
    initial_fail=[78],
    failtime=0.1,
    trip_lines = :dynamic,
    trip_nodes = :none,
    f_min = -0.03,
    f_max = 0.03,
    terminate_steady_state=true,
    solverargs = solverargs,
    warn=true
    );
#= 
Run simulation for trips on lines/nodes [78]
Shutdown line 78 at t = 0.1
Shutdown line 199 at t = 5.953557774825518
Shutdown line 194 at t = 7.796003380028169
Shutdown line 140 at t = 9.439939039593673
Shutdown line 131 at t = 17.06531885630484
Shutdown line 147 at t = 23.112060928377222
Shutdown line 106 at t = 31.060847771996002
Shutdown line 146 at t = 34.381726025591774
Shutdown line 14 at t = 38.32915069233711
Shutdown line 135 at t = 39.13784789668635
Shutdown line 149 at t = 41.02532791894525
Shutdown line 148 at t = 42.28700840247944
Terminated on steady state at 713.0886638569665
=#  

# nodes
simulate(network;
    verbose=true,
    x_static=x_static_old_ND,
    initial_fail=[78],
    failtime=0.1,
    trip_lines = :none,
    trip_nodes = :dynamic,
    f_min = -0.01, #NOTE This parameter has changed
    f_max = 0.01,
    terminate_steady_state=true,
    solverargs = solverargs,
    warn=true
    );
#= 
Run simulation for trips on lines/nodes [78]
Shutdown line 78 at t = 0.1
Shutdown node 29 at t = 2.2241521711760575
Shutdown node 98 at t = 2.5297865925232537
Terminated on steady state at 551.3118784026026
=#


