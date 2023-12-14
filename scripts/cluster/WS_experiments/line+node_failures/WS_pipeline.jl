# maybe use hyperparameter

# see scripts/archive/cascades_pool.jl
lines = 1:2
gen_τs   = [1.2,3.0]
slack_τs = 0.1:0.1:0.2
load_τs  = 0.1:0.1:0.2
hyperparameter = collect(Iterators.product(lines, gen_τs, slack_τs, load_τs))[:]

println(hyperparameter[1])



# TODO evtl. scripts/WS_local/line+node_failures/231211_01_WS_L_test_lines+nodes.jl nochal zerlegen

# TODO make structure usable for PC-Pool, PIK-HPC, Local PC
# TODO use HWs structure
