"""
See scripts/archive_nb/porting_to_new_ND/tests_old_ND.jl
"""

# TODO automate tests (with CI?)

include(abspath(@__DIR__, "..", "scripts/helpers_jarray.jl"))

using Serialization
using CairoMakie
CairoMakie.activate!()
using Dates


###################################################################################
###################################### RTS ########################################
###################################################################################

###
### steady state
###

"""
Uncomment for saving RTS steady state
"""
# network = import_system(:rtsgmlc; damping=0.1u"s", scale_inertia=1, tconst=0.01u"s")
# x_static = steadystate(network;
#                                 verbose=true,
#                                 zeroidx=1,
#                                 res_tol=1e-7,
#                                 relax_init_guess = true,
#                                 solverargs=(;reltol=1e-12, abstol=1e-12)
#                                 )  

# # write steady state to .csv
# filepath = abspath(@__DIR__, "..", "data/RTS-GMLC", string("steady_state_RTS-GMLC" , Dates.format(now(), "_yyyymmdd"), ".csv"))
# CSV.write(filepath, Dict(:SteadyState => x_static))


"""
###
### Choice of ODE solver
###
For choosing solver and `reltol` I set `reltol=1e-10`, `create_reference_solution = true`,
`choose_solver = true` and created the files in the directory `sol_objects_test/RTS_solver_choice/`. 
I then set `create_reference_solution = false` and increased `reltol` until there was a single test that 
was not successful. For the parameters below, this happened at reltol=1e-7. So choose reltol=1e-8.

# varying parameters
gen_model_array = [SwingDynLoadModel, SwingDynLoadModel_change_Pmech_only, SwingDynLoadModel_change_to_BH_only]
failure_modes_array = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
inertia_values = [0.2, 1.0, 20.0]
initial_fail_array = [27,29]

# fixed parameters
freq_bound = 0.42
damping=0.1u"s"
tconst=0.01u"s"
"""


"""
###
### Test for RTS simulations
###
Below few short tests for RTS simulations against a reference solution.

NOTE The referece solution from 20250604 is tested in `test_new_ND/WS_test_against_old_ND.jl` 
against the legacy version in branch `mwe_old_ND_maybe_plots` /scripts/archive_nb/porting_to_new_ND/tests_old_ND.jl 
"""

datetime = "20250604" # date of reference solution 
create_reference_solution = false # `false`: the recent code is compared to previously created reference solution
choose_solver = false
quick_test = true

## solver settings
using OrdinaryDiffEq
solver = Rodas4P(); solver_string = "Rodas4P()"
reltol=1e-8 # NOTE choose this tolerance for simulations
abstol=1e-6
solverargs = (;reltol=reltol, abstol=abstol) # default: reltol=1e-3, abstol=1e-6

## parameters to be tested
# varying parameters
if quick_test == true
    failure_modes_array = [[:dynamic, :dynamic]]
    inertia_values = [0.2]
else
    failure_modes_array = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
    inertia_values = [0.2, 20.0]
end

# fixed parameters
gen_model_array = [SwingDynLoadModel]
initial_fail_array = [27]
freq_bound = 0.42
damping=0.1u"s"
tconst=0.01u"s"

# load stead state 
x_static = DataFrame(CSV.File(abspath(@__DIR__, "..", "data/RTS-GMLC/steady_state_RTS-GMLC_20250604.csv"))).SteadyState

for 
    gen_model in gen_model_array,
    failure_mode in failure_modes_array,
    initial_fail in initial_fail_array,
    scale_inertia in inertia_values

    println("gen_model=$gen_model")
    println("failure_mode=$failure_mode")
    println("initial_fail=$initial_fail")
    println("scale_inertia=$scale_inertia")

    network = import_system(:rtsgmlc; damping=damping, scale_inertia=scale_inertia, tconst=tconst)

    sol = simulate(network;
        gen_model=gen_model,
        x_static=x_static,
        initial_fail=initial_fail,
        trip_lines = failure_mode[1],
        trip_nodes = failure_mode[2],
        freq_bound = freq_bound,
        solver = solver,
        solverargs = solverargs
        );

    filename = "RTS_$datetime,damping=$damping,scale_inertia=$scale_inertia," *
           "tconst=$tconst,gen_model=$gen_model,initial_fail=$initial_fail,trip_lines=$(failure_mode[1])," *
           "trip_nodes=$(failure_mode[2]),freq_bound=$freq_bound,solver=$solver_string"
    # for choosing solver and comparing different solver args, `solverargs` is missing in following string 
    if choose_solver == false
        filename = string("solverargs=$solverargs",filename)
    end
    filepath = joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,".sol")))

    if create_reference_solution == true
        # save reference solution
        Serialization.serialize(filepath, sol)
        save(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_ref.png"))), plot_simulation(sol))
    else
        # test against reference solution
        ref_sol = Serialization.deserialize(filepath)
        failure_times_ref_sol = describe_failures(ref_sol).failure_time
        failure_times_sol = describe_failures(sol).failure_time
        @assert length(failure_times_ref_sol) == length(failure_times_sol) "Different number of failures!"
        atol = 1e-4
        @assert isapprox(failure_times_ref_sol, failure_times_sol, atol=atol) "Failure times vary more than atol=$atol"
        println("=========== RTS Test successful :) =========== \n")
        # plot trajectories of recent code allowing a quick check
        save(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename, Dates.format(now(), "_yyyymmdd_HH"),".png"))), plot_simulation(sol))
    end
end


###################################################################################
####################################### WS ########################################
###################################################################################

"""
###
### Choice of ODE solver
###
For choosing solver and `reltol` I set `reltol=1e-10`, `create_reference_solution = true`,
`choose_solver = true` and created the files in the directory `sol_objects_test/WS_solver_choice/`. 
I then set `create_reference_solution = false` and increased `reltol` until there was a single test that 
was not successful. For the parameters below, this happened at reltol=1e-6. So choose to be save reltol=1e-8.

# CHECK 
 - Before doing larger simulations consider to also check different mathematical graphs (use different graph_seed_values
   and distr_seed_values)
 - Also, SwingDynLoadModel_change_Pmech_only, SwingDynLoadModel_change_to_BH_only have not been checked.
 - Maybe check whether other solvers are more performant.

# varying parameters
gen_model_array = [SwingDynLoadModel]
failure_modes_array = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
inertia_values = [0.2, 30.0]
damping_values = [0.2, 30.0]
β_values = [0.5]
graph_seed_values = distr_seed_values = [4]

# fixed parameters
freq_bound = 0.02
tconst = 1u"s"
initial_fail = 1 #NOTE Different `initial_fail` not needed as we have different mathematical graphs.
"""




"""
###
### Test for WS simulations
###
Below few short tests for WS simulations against a reference solution.
"""
datetime = "20250606" # date of reference solution 
create_reference_solution = false # `false`: the recent code is compared to previously created reference solution
choose_solver = false
quick_test == true

## solver settings
solver = Rodas4P(); solver_string = "Rodas4P()"
reltol=1e-8 # NOTE choose this tolerance for simulations
abstol=1e-6
solverargs = (;reltol=reltol, abstol=abstol) # default: reltol=1e-3, abstol=1e-6)

## parameters to be tested
# varying parameters
gen_model_array = [SwingDynLoadModel]
if quick_test == true
    failure_modes_array = [[:dynamic, :dynamic]]
else
    failure_modes_array = [[:dynamic, :dynamic], [:dynamic, :none], [:none, :dynamic]]
end
inertia_values = [30.0]
damping_values = [0.2]
β_values = [0.5]
graph_seed_values = distr_seed_values = [4]

# fixed parameters
freq_bound = 0.02
tconst = 1u"s"
initial_fail = 1 #NOTE Different `initial_fail` not needed as we have different mathematical graphs.

for 
    gen_model in gen_model_array,
    failure_mode in failure_modes_array,
    inertia in inertia_values,
    damping in damping_values,
    β in β_values,
    graph_seed in graph_seed_values,
    distr_seed in distr_seed_values

    println("gen_model=$gen_model")
    println("failure_mode=$failure_mode")
    println("inertia=$inertia")
    println("damping=$damping")
    println("β=$β")
    println("graph_seed=$graph_seed")
    println("distr_seed=$distr_seed")
    println("f_b=$freq_bound")


    filename = "WS_$datetime,damping=$damping,inertia=$inertia," *
           "tconst=$tconst,gen_model=$gen_model,trip_lines=$(failure_mode[1])," *
           "trip_nodes=$(failure_mode[2]),freq_bound=$freq_bound," *
           "β=$β,graph_seed=$graph_seed,distr_seed=$distr_seed,solver=$solver_string"
    # for choosing solver and comparing different solver args, `solverargs` is missing in following string 
    if choose_solver == false
        filename = string("solverargs=$solverargs",filename)
    end
    filepath = joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,".sol")))

    if create_reference_solution == true
        network = import_system(:wattsstrogatz; N=100, k=4, β=β, graph_seed=graph_seed,
                                        μ=0, σ=1, distr_seed=distr_seed,
                                        K=3, α=0.7, M=inertia*1u"s^2", γ=damping*1u"s", τ=tconst)

        # save referece network
        savegraph(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,".mg"))), network)

        # save reference graph
        savegraph(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_graph"))), network.graph)

        # save reference power injections
        Pmech = get_prop(network, 1:nv(network), :Pmech)
        Pload = get_prop(network, 1:nv(network), :Pload)
        CSV.write(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_power_injections.csv"))), Dict(:Pmech => Pmech, :Pload => Pload))

        # calculate reference steady state
        x_static = steadystate(network; verbose=true, res_tol=1e-5, relax_init_guess=false, zeroidx=1)  

        # save reference steady state to .csv
        CSV.write(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_steady_state.csv"))), Dict(:SteadyState => x_static))

        # calc reference solution
        sol = simulate(network;
            gen_model=gen_model,
            x_static=x_static,
            initial_fail=initial_fail,
            trip_lines = failure_mode[1],
            trip_nodes = failure_mode[2],
            freq_bound = freq_bound,
            solver = solver,
            solverargs = solverargs
            );

        # save reference solution
        Serialization.serialize(filepath, sol)

        # save reference plot
        save(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_ref.png"))), plot_simulation(sol))

    else
        # load referece network
        network = loadgraph(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,".mg"))), MGFormat())
        
        # load reference graph
        graph = loadgraph(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_graph"))))
        
        # load reference power injections
        Pmech = DataFrame(CSV.File(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_power_injections.csv"))))).Pmech
        Pload = DataFrame(CSV.File(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_power_injections.csv"))))).Pload
        set_prop!(network, 1:nv(network), :Pmech, Pmech)
        set_prop!(network, 1:nv(network), :Pload, Pload)
        
        # load reference steady state
        x_static = DataFrame(CSV.File(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename,"_steady_state.csv"))))).SteadyState
        
        # calc solution to be tested
        sol = simulate(network;
            graph=graph,
            gen_model=gen_model,
            x_static=x_static,
            initial_fail=1,
            trip_lines = failure_mode[1],
            trip_nodes = failure_mode[2],
            freq_bound = freq_bound,
            solver = solver,
            solverargs = solverargs
            );
        
        # test against reference solution
        ref_sol = Serialization.deserialize(filepath)
        failure_times_ref_sol = describe_failures(ref_sol).failure_time
        failure_times_sol = describe_failures(sol).failure_time
        @assert length(failure_times_ref_sol) == length(failure_times_sol) "Different number of failures!"
        atol = 1e-4
        @assert isapprox(failure_times_ref_sol, failure_times_sol, atol=atol) "Failure times vary more than atol=$atol"
        println("=========== WS Test successful :) =========== \n")
        # plot trajectories of recent code allowing a quick check
        save(joinpath(abspath(@__DIR__,"sol_objects_test", string(filename, Dates.format(now(), "_yyyymmdd_HH"),".png"))), plot_simulation(sol))
    end
end

println("======================================================================================================= \n")
println("=================================== RTS and WS successfully tested :) ================================= \n")
println("======================================================================================================= \n")