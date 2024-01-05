
Generate new network for ArrayTaskID=2 for the parameters failure_mode=[:dynamic, :dynamic],freq_bound=0.02,N=100,k=4,β=0.1,graph_seed=1,μ=0,σ=1,distr_seed=1,K=6,α=0.7,M=5.0 s^2,γ=1 s,τ=1 s,init_pert=line
ArrayTaskID=2 with parameters failure_mode=[:dynamic, :dynamic],freq_bound=0.02,N=100,k=4,β=0.1,graph_seed=201,μ=0,σ=1,distr_seed=201,K=9,α=0.7,M=0.5 s^2,γ=1 s,τ=1 s,init_pert=line.

No static solution found within tolerance for ArrayTaskID=5 with parameters failure_mode=[:dynamic, :dynamic],freq_bound=0.02,N=100,k=4,β=0.5,graph_seed=1,μ=0,σ=1,distr_seed=1,K=6,α=0.7,M=0.2 s^2,γ=1 s,τ=1 s,init_pert=line.
: No static solution found: ArrayTaskID=29 with parameters N=100,k=4,β=0.5,graph_seed=1,μ=0,σ=1,distr_seed=1,K=6,M=0.7,γ=1 s,τ=1 s.

,N=100,k=4,β=0.5,graph_seed=228,μ=0,σ=1,distr_seed=228,K=6,α=0.7,M=0.2 s^2,γ=1 s,τ=1 s,init_pert=line
N=100,k=4,β=0.9,graph_seed=10181,μ=0,σ=1,distr_seed=10181,K=6,M=0.2,γ=1 s,τ=1 s


N=100,k=4,β=0.9,graph_seed=6,μ=0,σ=1,distr_seed=6,K=6,M=0.2 s^2,γ=1 s,τ=1 s





savegraph("graph.lg", network)

network

bla = loadgraph("graph.lg")



cd("/home/vollmich/.julia/dev/DynamicCascades")
pwd()


g = path_graph(5)
mg = MetaGraph(g, 3.0)
set_prop!(mg, :description, "This is a metagraph.")

i=1.0
savegraph("/home/vollmich/Desktop/delete_soon/test_graph.lg", network)


mg1 = loadgraph("/home/vollmich/Desktop/delete_soon/test_graph.lg",MGFormat())

rng = MersenneTwister(1234);

uuid4()
UUID("7a052949-c101-4ca3-9a7e-43a2532b2fa8")

@time for i in [0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    println("i $i \n")
    network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(i * 1u"s^2"), graph_seed=124, distr_seed=1230, K=1, γ=1u"s", τ=1u"s", σ=1.0)
    savegraph("/home/vollmich/Desktop/delete_soon/test_graphs/test_graph$i.lg", network)
    # x_static = steadystate(network)
end


lines = 1:2
gen_τs   = 0.1:0.1:0.2
slack_τs = 0.1:0.1:0.2
load_τs  = 0.1:0.1:0.2
hyperparameter = collect(Iterators.product(lines, gen_τs, slack_τs, load_τs))[:]


################################################################################

#SBATCH --qos=priority
#SBATCH --time=1-00:00:00
#SBATCH --job-name=231214_01_WS_C_test_lines+nodes_jarray
#SBATCH --output=231214_01_WS_C_test_lines+nodes_jarray.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2


#bl SBATCH --nodes=1
#bl SBATCH --mail-type=begin        # send email when job begins
#bl SBATCH --mail-type=end          # send email when job ends
#bl SBATCH --mail-user=brandner@pik-potsdam.de



##################################

savegraph("graph.lgz", g)

g = loadgraph(string(directories[1],"/graph.lgz"))

s. MM Julia coding vom [2022-11-06 So]

#SBATCH --qos=medium
#SBATCH --job-name=N020
#SBATCH --account=icone
#SBATCH --output=%x-%j-%N.out
#SBATCH --error=%x-%j-%N.err
#SBATCH --ntasks=100
#SBATCH --time=4-0


echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia
julia script.jl $SLURM_NTASKS




DataFrame(
    ArrayTaskID = 1:length(hyperparam_ensemble),
    Beta = [x[1] for x in hyperparam_ensemble],
    Inertia_Values = [x[2] for x in hyperparam_ensemble],
    Freq_Bounds = [x[3] for x in hyperparam_ensemble],
    Failure_Modes = [x[4] for x in hyperparam_ensemble],
    K = [x[5] for x in hyperparam_ensemble],
    N_Nodes = [x[6] for x in hyperparam_ensemble],
    Coupling_K = [x[7] for x in hyperparam_ensemble],
    Gamma = [x[8] for x in hyperparam_ensemble],
    Tau = [x[9] for x in hyperparam_ensemble],
    Alpha = [x[10] for x in hyperparam_ensemble],
    Init_Pert = [x[11] for x in hyperparam_ensemble],
    Sigma = [x[12] for x in hyperparam_ensemble],
    Mu = [x[13] for x in hyperparam_ensemble]
)


df_hyperparam_ensemble = DataFrame()


df_hyperparam_ensemble.ArrayTaskID = 1:length(hyperparam_ensemble)
df_hyperparam_ensemble


df_hyperparam_ensemble[!, :ArrayTaskID] = 1:length(hyperparam_ensemble)

df_hyperparam_ensemble

new_column_df = DataFrame(ArrayTaskID = 1:length(hyperparam_ensemble))





using Distributions
μ = 0
σ = 1.0
N = 20
d = Distributions.Normal(μ, σ)
x = rand(d, 2N)

inertia = 2.0u"s^2"
inertia
bla = inertia * 1u"s^2"
network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(inertia * 1u"s^2"), graph_seed=124, distr_seed=1230, K=6, γ=(1 * 1u"s"), τ=(1 * 1u"s"), μ=0, σ=1)
steadystate(network)

network
get_prop(network, 2, :_M)

i = 10
set_prop!(network, 1:nv(network), :_M, i * 1u"s^2")


M=(i * 1u"s^2")

savegraph("/home/vollmich/Desktop/delete_soon/test_graphs/test_graph$i.lg", network)

##############
sdir(dir) || mkdir(dir)

for i in [0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
    println("i $i \n")
    network = import_system(:wattsstrogatz; N=100, k=4, β=0.7, M=(i * 1u"s^2"), graph_seed=124, distr_seed=1230, K=K, γ=(gamma * 1u"s"), τ=(tau * 1u"s"), μ=mu, σ=sigma)
    savegraph("/home/vollmich/Desktop/delete_soon/test_graphs/test_graph$i.lg", network)
    # x_static = steadystate(network)
end
################

string(2.0)

while task_id <= 10
    string(2.0)
    task_id += 1
end


for i = 1:10
   if i % 3 != 0
       continue
   end
   println(i)
end

i = 0
while i <= 10
   i += 1
   if i % 3 != 0
       continue
   end
   println(i)
end


i = 0
while i <= 10
    i += 1
    println(i)

    for j in 1:5
        try
            network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, M=0.7u"s^2", graph_seed=1, distr_seed=1, K=6.0, γ=1u"s", τ=1u"s", α= 0.7, σ=1.0, μ=0)
            x_static = steadystate(network)
        catch
            println("didnt work")
            continue
        end

        if i % 3 != 0
           continue
        end

end

network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, M=0.7u"s^2", graph_seed=4, distr_seed=4, K=6.0, γ=1u"s", τ=1u"s", α= 0.7, σ=1.0, μ=0)
x_static = steadystate(network)


####################
task_id = 5

string_args = string_metagraph_args(df_hpe, task_id)
println("Generate new MetaGraph:ArrayTaskID=$task_id with parameters $string_args")

graph_seed = 1; distr_seed = 1
#
#
df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
network = import_system_wrapper(df_hpe, task_id)

max_trials = 20
trial_counter = 1
steady_state_for_all_inertia_values = false
while steady_state_for_all_inertia_values == false
    try
        # Try all inertia values
        for i in inertia_values
            set_prop!(network, 1:nv(network), :_M, i * 1u"s^2")
            println("Test steady state with inertia I=$i...")
            steadystate(network)
        end

        # If (in case of no error in previous for loop) a steady state exists for all inertia values
        trial_counter = 1 # reset trial_counter
        steady_state_for_all_inertia_values = true

    # If a there is at least one inertia value for which no steady state exists
    catch
        trial_counter += 1
        if trial_counter > max_trials
            error("Tried $max_trials different values for `graph_seed` and `distr_seed`. Exiting...")
        end

        N,k,β,graph_seed,μ,σ,distr_seed,K,_,_,γ,τ,_,_,_ = get_network_args(df_hpe, task_id)
        M = ustrip(u"s^2", get_prop(network, 1, :_M))
        @warn "No static solution found: ArrayTaskID=$task_id with parameters N=$N,k=$k,β=$β,graph_seed=$graph_seed,μ=$μ,σ=$σ,distr_seed=$distr_seed,K=$K,M=$M,γ=$γ,τ=$τ."

        # generate new grid by increasing seeds by one
        graph_seed += 1; distr_seed += 1
        df_hpe[task_id,:graph_seed] = graph_seed; df_hpe[task_id,:distr_seed] = distr_seed
        string_args = string_network_args(df_hpe, task_id)
        println("Generate new network with changed seeds for ArrayTaskID=$task_id with parameters $string_args")
        network = import_system_wrapper(df_hpe, task_id)
    end
end

#################
i = 0
while i <= 10
   i += 1
   println("$i \n")
   if i % 3 == 0
       break
   end
   print("Hallo \n")
end

steady_state_for_inertia_values = false
i = 0
while steady_state_for_inertia_values != true
   i += 1
   println("$i \n")
   if i == 3
       steady_state_for_inertia_values = true
       break
   end
   print("Hallo \n")
end


number_of_task_ids_between_graphs = length(inertia_values) * length(freq_bounds) * length(failure_modes)

# Loop over each ArrayTaskID:
task_id = 1
while task_id <= length(df_hpe.ArrayTaskID)
    #= For every new parameter configuration (excluding the variation of the inertia
    parameter), generate a new network:=#
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        string_args = string_network_args(df_hpe, task_id)
        println("Generate new network for ArrayTaskID=$task_id with parameters $string_args")
        network = import_system_wrapper(df_hpe, task_id)

        #= Check for every inertia value of a single parameter configuration if steady
        state within tolerance exists. If for one inertia value no steady state exists,
        go back and generate new network.=#
        trial_counter = 1
        steady_state_for_all_inertia_values = false
        while steady_state_for_all_inertia_values != true
            if steady_state_for_all_inertia_values == true
                break
            end

            i = 0
            while i < length(inertia_values)
                i += 1
                set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2")

                try
                    steadystate(network) # test if steady state exists
                    trial_counter = 1

                    #     # save grid with parameter configuration if last steady state for last inertia value exists
                    #     if ((task_id-1) % length(inertia_values)) == (length(inertia_values) - 1)
                    #
                        # create folder if not already existing
                        # isdir(dir_path) || mkdir(dir_path)
                    #         #
                    #         # # save network
                    #         # println("Save network for jobs with ArrayTaskIDs $(task_id-length(inertia_values)+1) to $task_id.")
                    #         # savegraph(filepath, network)
                    #
                        # graph_seed += 1; distr_seed += 1 NOTE man braucht keinen neuen seed, weil sich für jede task_id die parameter ändern
                    #     end
                    #     task_id += 1
                    # # in case `steadystate(network)` throws exeption

                    if i == length(inertia_values)
                        steady_state_for_inertia_values = true
                    end

                catch
                    string_args = string_network_args(df_hpe, task_id)
                    println("No static solution found within tolerance for ArrayTaskID=$task_id with parameters $string_args.")

                    # generate new grid by changing seeds
                    df_hpe[task_id,:graph_seed] += 1; df_hpe[task_id,:distr_seed] += 1
                    string_args = string_network_args(df_hpe, task_id)
                    println("Generate new network with changed seeds for ArrayTaskID=$task_id with parameters $string_args")
                    network = import_system_wrapper(df_hpe, task_id)

                    # # start again with first inertia value
                    # task_id -= ((task_id-1) % length(inertia_values)) # TODO dann wieder zurück zu Zeile 65
                    i = 0


                    if trial_counter == 20
                        error("Tried 20 grids => exiting")
                    end
                    trial_counter += 1 # TODO somehow while loop in here until steady state found. For each `catch` maximum tries
                    println("trial_counter=$trial_counter")
                    exit("Bis hier")
                end
            end
        # Save/overwrite filepath, graph_seed, distr_seed in df
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)

        # Create paths
        t=now()
        datetime = Dates.format(t, "yyyymmdd_HHMMSS.s")
        dir_params = "failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β,"
        exp_path = joinpath(@__DIR__, exp_name)
        ispath(exp_path) || mkdir(exp_path)
        dir_path = joinpath(@__DIR__, exp_name, dir_params)
        ispath(dir_path) || mkdir(dir_path)
        graph_params = "graph_seed=$graph_seed,distr_seed=$distr_seed,failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β,"
        filepath = joinpath(dir_path, "graphs", string(graph_params,datetime,".lg"))

        # Assign filepath to df
        df_hpe[task_id,:filepath] = filepath
        end

    # else
    # set_prop!(network, 1:nv(network), :_M, df_hpe[task_id,:inertia_values])
    end

    task_id += 1
end

steady_state_for_inertia_values = false
trial_counter = 1
i = 0
while i < 4
    i += 1
    try
        trial_counter = 1
        sd
        println("try $i")
        if i == 4
            steady_state_for_inertia_values = true
        end
    catch
        trial_counter += 1
        if trial_counter == 20
            error("Tried 20 grids => exiting")
        end
        i = 0
        println("trial_counter $trial_counter")
        error("Tried 20 grids => exiting")
    end
end

# for i in 1:24
#     println(i)
#     network = import_system_wrapper(df_hpe, i)
#     steadystate(network)
# end

# set_prop!(network, 1:nv(network), :_M, 10)

# df_hpe[1,:inertia_values]

# network = import_system_wrapper(df_hpe, 5)
# steadystate(network)


# # For adjusting how many lines/colums are printed
# ENV["COLUMNS"]=200
# ENV["LINES"]=100



################################################################################
################################################################################

# GENERATION OF NETWORKS
number_of_task_ids_between_graphs = length(inertia_values) * length(freq_bounds) * length(failure_modes)

# Loop over each ArrayTaskID:
task_id = 1
while task_id <= length(df_hpe.ArrayTaskID)
    #= For every new parameter configuration (excluding the variation of the inertia
    parameter), generate a new network:=#
    if ((task_id-1) % number_of_task_ids_between_graphs) == 0
        string_args = string_network_args(df_hpe, task_id)
        println("Generate new network for ArrayTaskID=$task_id with parameters $string_args")
        network = import_system_wrapper(df_hpe, task_id)

        #= Check for every inertia value of a single parameter configuration if steady
        state within tolerance exists. If for one inertia value no steady state exists,
        go back and generate new network.=#
        exit_id = 1
        steady_state_for_all_inertia_values = false
        while steady_state_for_all_inertia_values != true
            i = 0
            while i < length(inertia_values)
                i += 1
                set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2")
                try
                    steadystate(network) # test if steady state exists
                    exit_id = 1

                    #     # save grid with parameter configuration if last steady state for last inertia value exists
                    #     if ((task_id-1) % length(inertia_values)) == (length(inertia_values) - 1)
                    #
                        # create folder if not already existing
                        # isdir(dir_path) || mkdir(dir_path)
                    #         #
                    #         # # save network
                    #         # println("Save network for jobs with ArrayTaskIDs $(task_id-length(inertia_values)+1) to $task_id.")
                    #         # savegraph(filepath, network)
                    #
                        # graph_seed += 1; distr_seed += 1 NOTE man braucht keinen neuen seed, weil sich für jede task_id die parameter ändern
                    #     end
                    #     task_id += 1
                    # # in case `steadystate(network)` throws exeption
                    if i == length(inertia_values)
                        steady_state_for_inertia_values = true
                    end
                catch
                    if exit_id == 20
                        error("Tried 20 grids => exiting")
                    end
                    exit_id += 1 # TODO somehow while loop in here until steady state found. For each `catch` maximum tries
                    println("exit_id=$exit_id")
                    string_args = string_network_args(df_hpe, task_id)
                    println("No static solution found within tolerance for ArrayTaskID=$task_id with parameters $string_args.")

                    # generate new grid by changing seeds
                    df_hpe[task_id,:graph_seed] += 1; df_hpe[task_id,:distr_seed] += 1
                    string_args = string_network_args(df_hpe, task_id)
                    println("Generate new network with changed seeds for ArrayTaskID=$task_id with parameters $string_args")
                    network = import_system_wrapper(df_hpe, task_id)
                    # # start again with first inertia value
                    # task_id -= ((task_id-1) % length(inertia_values)) # TODO dann wieder zurück zu Zeile 65
                    i = 0
                end
            end
        # Save/overwrite filepath, graph_seed, distr_seed in df
        N,k,β,graph_seed,μ,σ,distr_seed,K,α,M,γ,τ,freq_bound,failure_mode,init_pert = get_network_args(df_hpe, task_id)

        # Create paths
        t=now()
        datetime = Dates.format(t, "yyyymmdd_HHMMSS.s")
        dir_params = "failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β,"
        exp_path = joinpath(@__DIR__, exp_name)
        ispath(exp_path) || mkdir(exp_path)
        dir_path = joinpath(@__DIR__, exp_name, dir_params)
        ispath(dir_path) || mkdir(dir_path)
        graph_params = "graph_seed=$graph_seed,distr_seed=$distr_seed,failure_mode=$failure_mode,freq_bound=$freq_bound,k=$k,β=$β,"
        filepath = joinpath(dir_path, "graphs", string(graph_params,datetime,".lg"))

        # Assign filepath to df
        df_hpe[task_id,:filepath] = filepath
        end
    # else
    # set_prop!(network, 1:nv(network), :_M, df_hpe[task_id,:inertia_values])
    end
    task_id += 1
end

# TODO unterschiedlicher Seed für alle N_G grids

# Save to CSV
CSV.write(joinpath(@__DIR__, string(exp_name, "/config.csv")), df_hpe)



# TODO change inertia value of graph
# TODO use correct inertia value via set_prop!(network, 1:nv(network), :_M, inertia_values[i] * 1u"s^2")

# scale_inertia_values = [0.2, 1.0, 2.1, 3.0, 4.0, 5.0, 7.2, 11.0, 15.0, 21.0]
scale_inertia_values = [0.2, 1.0] # varying parameter
df_all_failures = DataFrame()
df_all_failures_nodes = DataFrame()
@time for scale_inertia in scale_inertia_values
    network = import_system(:wattsstrogatz; N=N, k=4, β=0.7, M=(inertia * 1u"s^2"), graph_seed=124, distr_seed=1230, K=1, γ=1u"s", τ=1u"s", σ=1.0)
    x_static = steadystate(network)
    println("Scaling of inertia $scale_inertia \n steady state \n $x_static \n ")
    number_failures = Float64[]
    number_failures_nodes = Float64[]
    # for i in 1:ne(network)
    for i in 1:2 # TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODOREMOVE
        sol = simulate(network;
                       initial_fail = Int[i],
                       init_pert = :line,
                       tspan = (0, 10),
                       trip_lines = :dynamic,
                       trip_nodes = :dynamic,
                       trip_load_nodes = :none,
                       f_min = -freq_bound,
                       f_max = freq_bound,
                       solverargs = (;dtmax=0.01),
                       verbose = true);
        push!(number_failures, length(sol.failures.saveval)-1) # `-1` as we don't want to count the initial failure
        push!(number_failures_nodes, length(sol.failures_nodes.saveval))
    end
    df_all_failures[!, string(scale_inertia)] = number_failures
    df_all_failures_nodes[!, string(scale_inertia)] = number_failures_nodes
end


a = 1:5


b = ["e","d","b","c","a"]


for (a,b) in zip(a,b)
    println("a=$a, b=$b \n")
end




network = import_system(:wattsstrogatz; N=20, k=4, β=0.7, M=2.1u"s^2", graph_seed=124, distr_seed=1230, K=1.0, γ=1u"s", τ=1u"s", α= 0.7, σ=1.0, μ=0)
x_static = steadystate(network)
# AssertionError: No steady state found 5.896376410463589e-7

################################################################################
# Different way:




condition = let _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
    (out, u, t, integrator) -> begin
        # upcrossing through zero triggers condition
        calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
        # calculate_active_power(_current_load, u, integrator.p, t, integrator.f.f, _network)
        out .= _current_load .- _rating
        nothing


if trip_lines == :dynamic
    if monitored_power_flow == :apparent
        condition = let _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
            (out, u, t, integrator) -> begin
                # upcrossing through zero triggers condition
                calculate_apparent_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                out .= _current_load .- _rating
                nothing
            end
        end
    elseif monitored_power_flow == :active
        condition = let _current_load = zeros(ne(network)), _network = network, _rating = get_prop(network, edges(network), :rating)
            (out, u, t, integrator) -> begin
                # upcrossing through zero triggers condition
                calculate_active_power!(_current_load, u, integrator.p, t, integrator.f.f, _network)
                out .= _current_load .- _rating
                nothing
            end
        end
    end

monitored_power_flow =:apparent
