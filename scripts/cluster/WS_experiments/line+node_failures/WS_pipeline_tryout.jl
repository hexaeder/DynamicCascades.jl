################################################################################
network = import_system(:wattsstrogatz; N=N, k=4, β=0.7, M=1.0u"s^2", graph_seed=124, distr_seed=1230, K=1, γ=1u"s", τ=1u"s", σ=1.0)
x_static = steadystate(network)

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
