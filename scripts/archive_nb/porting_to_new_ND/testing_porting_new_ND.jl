"""
See scripts/archive_nb/porting_to_new_ND/testing_old_ND.jl
"""

###
### RTS: Test new ND w/o node failures against old ND w/o node failures (there was a bug in node failure CB)
###

## new ND
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))

damping = 0.1u"s"
scale_inertia = 1.2
tconst = 1.0u"s"
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

# get network that was created with old network to make sure to use the same network across the 2 versions
filepath = "/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND.mg"
network = loadgraph(filepath, MGFormat())

# steady_state_dict  = CSV.File("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND_steady_state.csv")
# x_static = steady_state_dict[:SteadyState]

sol = simulate(network;
    gen_model=SwingDynLoadModel,
    initial_fail=27,
    tspan=(0., 3.),
    trip_lines = :dynamic,
    trip_nodes = :none,
    freq_bound = 1.5,
    solverargs = (;dtmax=0.01),
    verbose = true);



###
### RTS: Test parameters
###

# p_new
nw = NetworkDynamics.extract_nw(sol)
p_new = describe_vertices(nw)
p_new = mapcols(col -> 
    # if the column’s element type can hold NaN, replace missings
    eltype(col) <: Union{Missing, Float64} ? coalesce.(col, NaN) : col,
    p_new
)

# p_old
p_old = Serialization.deserialize("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND_parameters_old_ND")

# NOTE Test parameters one by one by commenting out the other parameters. Snippet doesn't work otherwise. 
atol = 1e-16
for i in 1:nv(network)
    # power
    P_new = p_new.Pmech[i] - p_new.Pload[i]
    P_old = p_old[1][i][2]
    @assert isapprox(P_new, P_old, atol=atol) "i=$i,P_new=$P_new,P_old=$P_old"
    # inertia
    I_new = p_new.M[i]
    I_old = p_old[1][i][3]
    # @assert isapprox(I_new, I_old, atol=atol) "i=$i,I_new=$I_new,I_old=$I_old"
    @assert I_new === I_old "i=$i,I_new=$I_new,I_old=$I_old"
    # τ
    τ_new = p_new.τ[i]
    τ_old = p_old[1][i][4]
    @assert isapprox(τ_new, τ_old, atol=atol) "i=$i,τ_new=$τ_new,τ_old=$τ_old"
    # D
    D_new = p_new.D[i]
    D_old = p_old[1][i][5]
    # @assert isapprox(D_new, D_old, atol=atol) "i=$i,D_new=$D_new,D_old=$D_old"
    @assert D_new === D_old "i=$i,D_new=$D_new,D_old=$D_old"
end

e_new = NetworkDynamics.describe_edges(nw)

atol = 1e-8
for i in 1:ne(network)
    # K
    K_new = e_new.K[i]
    K_old = abs(p_old[2][i][1])
    @assert isapprox(K_new, K_old, atol=atol) "i=$i,K_new=$K_new,K_old=$K_old"
    # line rating
    rating_new = e_new.rating[i]
end

# rating
e_new.rating == ustrip(get_prop(network, edges(network), :rating))


###
### RTS: Plot trajectories
###

## inspector
nw = NetworkDynamics.extract_nw(sol)

vindices = [i for i in 1:nv(nw)]
eindices = [i for i in 1:ne(nw)]

inspect(sol; restart=true, reset=true)
set_sol!(sol) # optional if after inspect(sol)
set_graphplot!(; nstate=[:ω], estate=[:S], nstate_rel=false, estate_rel=false)
define_timeseries!([
    (; selcomp=[VIndex(i) for i in vindices], states=[:ω, :ωmax], rel=false),
    (; selcomp=[EIndex(i) for i in [27,29,37,11,12]], states=[:S, :rating], rel=false),
])

## static plot
using CairoMakie
using Colors
nw = NetworkDynamics.extract_nw(sol)

all_nodes_idxs = [i for i in 1:nv(nw)]
all_lines_idxs = [i for i in 1:ne(nw)]
node_colors = distinguishable_colors(length(all_nodes_idxs))
line_colors = distinguishable_colors(length(all_lines_idxs))


fontsize = 35
titlesize = (fontsize+5)
linewidth = 3.5
fig = Figure(size=(3100,1500), fontsize=fontsize)

# FREQUENCIES 
fig[1,1] = ax = Axis(fig; ylabel="Frequency [Hz]", title="Bla", titlealign = :left, titlesize = titlesize)
for i in all_nodes_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=vidxs(nw, i, "ω")))))
    color = node_colors[findfirst(x -> x == i, all_nodes_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end

# FLOWS 
fig[1,2] = ax = Axis(fig; xlabel="Time [s]")
for i in all_lines_idxs
    x, y = remove_zero_tail!(deepcopy(sol.t), deepcopy(map(first, sol(sol.t, idxs=eidxs(i, :S)))))
    color = line_colors[findfirst(x -> x == i, all_lines_idxs)]
    lines!(ax, x, y; label="Node $i", color=color, linewidth=linewidth) # consider plotdensity=Int(1e5)
    scatter!(ax, (x[end], y[end]); color=color, marker=:star5, markersize=25)
end
fig




###
### WS
###

## before adapting new ND to RTS
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))
include("/home/brandner/.julia/dev/for_testing_delete_soon/WS_trajectories_new_ND_single_model_port.jl")
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
sub_dir = "trajectories_MA_after_up3"
dir = joinpath(exp_data_dir, sub_dir)
task_id = 1728
initial_fail = 78
sol = simulate_new_ND_single_model_port(exp_data_dir, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);


## wrapper after adapting new ND to RTS
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
task_id = 1728
initial_fail = 78

sol = simulate(exp_name_date, task_id, initial_fail;
    gen_model=SwingDynLoadModel,
    tspan=(0., 40.),
    solverargs = (;dtmax=0.01),
    verbose = true);


###
### Fixpoint calculation for old ND and new ND
###

#= #NOTE The states in the backend of old and new ND have a different order for RTS:
In the new ND the states are ordered by components.=#

## old ND
zeroidx=1
(nd, p) = nd_model(network)
x0 = zeros(length(nd.syms));
x_static = solve(NonlinearProblem(nd, x0, p), NLSolveJL())

θidx = idx_containing(nd, "θ")
offset = x_static[θidx[zeroidx]]
x_static[θidx] .= x_static[θidx] .- offset
@assert iszero(x_static[θidx[zeroidx]])

## new ND
zeroidx=1
nw_state = NWState(nw)

x0 = zeros(dim(nw))
#= This does not work anymore in new DE.jl. Possible reason: The system defined by `nw` is an
ODE-System and not an nonlinear problem (which would be nonlinear algebraic equations). In the 
more recent version of DE.jl this throws an error to avoid an inconsistent use of `NonlinearProblem`,
i.e. trying to solve a nonlinear problem that in fact is not a nonlinear problem. Instead, the
natural way is to define a steady state problem.=#
x_static = solve(NonlinearProblem(nw, x0, pflat(nw_state)), NLSolveJL())

# This works
solve(NonlinearProblem((du,u,p)->nw(du,u,p,0.0), x0, pflat(nw_state)), NLSolveJL())

#= This throws an error, because `_solve_fixpoint` which is invoced by `find_fixpoint` expects
the type `AbstractNonlinearSolveAlgorithm`, which doesn't fit `typeof(NLSolveJL())`. =#
find_fixpoint(nw, alg=NLSolveJL()) 

# This works
solve(SteadyStateProblem(nw, x0, pflat(nw_state)), NLSolveJL())

## Getting rid off phase-shift
# (first try)
nw[EIndex(initial_fail)]
uflat_s0 = uflat(s0)
offset = uflat_s0[1]
uflat_s0[1:2:end] .-= offset
uflat_s0

# (Notes HW)
s0.v[1:100,:ω]  
s0.v[:,:θ]
s0.v[1:100,:θ]
s0.v[1:100,:θ] .- s0.v[1,:θ]

# compare steady state for old and new ND.jl
isapprox(steady_state, uflat_s0, atol=1e-3)

###
### Checking if new df_hpe is correct
###
# old
exp_name_date = "WS_k=4_exp05_2_I_over_Dsq_nodes_PIK_HPC_K_=3,N_G=32_20250125_125951.601"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
df_hpe_old = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
# new
exp_name_date = "WS_k=4_exp05_2_I_over_Dsq_nodes_PIK_HPC_K_=3,N_G=32_20250326_230941.02"
exp_data_dir = joinpath(RESULTS_DIR, "preprocessing_to_include_graphs", exp_name_date)
df_hpe_new = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))

for task_id in df_hpe_old.ArrayTaskID
    if string_network_args(df_hpe_old, task_id) != string_network_args(df_hpe_new, task_id)
        println("error in task_id=$task_id")
    end
end


###
### RTS get number of GFM and GFL nodes
### 
n_gen = 0
n_load = 0
for i in 1:nv(network)
    type = get_prop(network, i, :type)
    if type == :gen || type == :syncon
        n_gen += 1
    elseif type == :load
        n_load += 1
    else
        error("The node type :$type is not defined!")
    end
end
n_gen
n_load


#= Testing load node frequencies
see worktree-mwe_old_ND_maybe_plots/scripts/archive_nb/porting_to_new_ND/testing_load_node_frequencies.jl =#
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))

using CairoMakie
CairoMakie.activate!()

damping = 0.1u"s"
scale_inertia = 1.0
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = 0.01u"s")

sol = simulate(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    verbose=true,
    initial_fail=23,
    failtime=0.1,
    tspan=(0., 0.4),
    trip_lines = :none,
    trip_nodes = :dynamic,
    freq_bound = 1.0/(2π),
    terminate_steady_state=true,
    ODE_solver = AutoTsit5(Rodas5P()),
    solverargs = (;dtmax=0.01),
    warn=true
    );

# Find steady state...
# Shutdown line 23 at t=0.1
# Vertex 13 tripped at t=0.1469172276792874
# Vertex 23 tripped at t=0.31099962157331357
# Vertex 16 tripped at t=0.31833691057498176
# Vertex 42 tripped at t=0.32940018058770976
# Vertex 40 tripped at t=0.33317837677715806
# Vertex 45 tripped at t=0.33965519086272294
# Vertex 15 tripped at t=0.33991868504932626
# Vertex 14 tripped at t=0.34081334127495905
# Vertex 18 tripped at t=0.3459729641972824
# Vertex 21 tripped at t=0.34823223061975483
# Vertex 39 tripped at t=0.35065112403711446
# Vertex 38 tripped at t=0.3565746652845902
# Vertex 46 tripped at t=0.35980634486184143
# Vertex 7 tripped at t=0.3646663577750275
# Vertex 1 tripped at t=0.3652271316001347
# Vertex 22 tripped at t=0.36562555225807164
# Vertex 2 tripped at t=0.3664463031760093
# Vertex 71 tripped at t=0.368984831098414
# Vertex 47 tripped at t=0.3707490098684115
# Vertex 37 tripped at t=0.37254020096656265
# Vertex 64 tripped at t=0.3739680143292981
# Vertex 61 tripped at t=0.3743719190059114
# Vertex 66 tripped at t=0.37631382796702656
# Vertex 62 tripped at t=0.38003618118166643
# Vertex 69 tripped at t=0.38064777150439044
# Vertex 31 tripped at t=0.3808744685984422
# Vertex 25 tripped at t=0.38120897967242
# Vertex 26 tripped at t=0.3829854130302008
# Vertex 63 tripped at t=0.38603474224362244
# Vertex 70 tripped at t=0.39197160113918633
# Vertex 55 tripped at t=0.3929696292582295
# Vertex 49 tripped at t=0.3966172650202103
# Vertex 50 tripped at t=0.39728176953161987
# ┌ Warning: Did not reach steady state!

plot_simulation(sol)