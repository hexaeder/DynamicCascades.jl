"""
See scripts/archive_nb/porting_to_new_ND/testing_old_ND.jl
"""

###
### WS test
###
## `simulate_new_ND()` without wrapper after adapting new ND to RTS
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
# get graph
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
df_config = DataFrame(CSV.File(joinpath(RESULTS_DIR, exp_name_date, "config.csv")))
df_config[!, :filepath_graph] = replace.(df_config[!, :filepath_graph],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
df_config[!, :filepath_steady_state] = replace.(df_config[!, :filepath_steady_state],"/home/brandner" => "/home/brandner/nb_data/HU_Master/2122WS/MA")
task_id = 1728
graph = loadgraph(df_config[task_id,:filepath_graph])

network = import_system(:wattsstrogatz; N=100, k=4, β=0.5, graph_seed=11,
        μ=0, σ=1, distr_seed=11, K=3, α=0.7, M=30*1u"s^2", γ=1*1u"s", τ=1*1u"s")

simulate_new_ND(network;
    graph=graph, 
    gen_model=SwingDynLoadModel,
    initial_fail=78,
    failtime=0.1,
    tspan=(0., 40.),
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    freq_bound = 0.03,
    terminate_steady_state=true,
    ODE_solver = AutoTsit5(Rodas5P()),
    warn=true
    );
#= 
Shutdown line 78 at t=0.1
Line 199 tripped at t=5.953723915979403
Line 194 tripped at t=7.796300228395314
Vertex 98 tripped at t=8.850977737823337
Line 131 tripped at t=16.776726663702437
Line 147 tripped at t=22.48384108354435
Line 146 tripped at t=26.674319991326886
Line 14 tripped at t=28.794766730477853
Line 106 tripped at t=29.79965093880108
Line 135 tripped at t=30.012141236675394
Vertex 59 tripped at t=30.754639398734433
Line 149 tripped at t=32.6945286496743
Vertex 62 tripped at t=33.20901295100614
Line 148 tripped at t=33.29470516561573
Vertex 61 tripped at t=33.46960161580322
Vertex 60 tripped at t=33.76078662307933
┌ Warning: Did not reach steady state!
└ @ Main ~/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl:374
=#



###
### RTS test
###
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
# Rodas4P()
damping = 0.1u"s"
scale_inertia = 1.2
tconst = 1.0u"s"
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

# # For loading the same graph
# filepath = "/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND.mg"
# network = loadgraph(filepath, MGFormat())

sol = simulate_new_ND(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    verbose=true,
    initial_fail=27,
    failtime=0.1,
    tspan=(0., 10.),
    trip_lines = :dynamic,
    trip_nodes = :none,
    freq_bound = 1.5,
    terminate_steady_state=true,
    solverargs = (;dtmax=0.01),
    warn=true
    );
#= 
Shutdown line 27 at t=0.1
Line 29 tripped at t=0.4124845045417539
Line 37 tripped at t=0.7975984254665589
Line 11 tripped at t=1.3627955214633518
Line 12 tripped at t=1.5007722036394446
┌ Warning: Did not reach steady state!
=#

## smaller τ
# Rodas4P()
damping = 0.1u"s"
scale_inertia = 1.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

sol = simulate_new_ND(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    verbose=true,
    initial_fail=27,
    failtime=0.1,
    tspan=(0., 10.),
    trip_lines = :dynamic,
    trip_nodes = :none,
    freq_bound = 1.5,
    terminate_steady_state=true,
    solverargs = (;dtmax=0.0001),
    warn=true
    );
#= 
Shutdown line 27 at t=0.1
Line 29 tripped at t=0.19341022145383527
Line 11 tripped at t=0.30445308669178334
Line 12 tripped at t=0.44383125872221846
Line 37 tripped at t=0.4880584753358242
Line 71 tripped at t=0.796309218838602
Line 56 tripped at t=0.9961127774919499
Line 59 tripped at t=0.9967503604894374
Line 58 tripped at t=1.0152130243535804
Line 41 tripped at t=1.1016889483165229
Line 49 tripped at t=1.1037063853209952
Line 43 tripped at t=1.1067824074479864
Line 42 tripped at t=1.1085461956188258
Line 44 tripped at t=1.1128117794825694
Line 57 tripped at t=1.1185219507148567
Line 100 tripped at t=1.1941321724462979
Line 97 tripped at t=1.2390820131019546
Line 98 tripped at t=1.2431974978517988
Line 73 tripped at t=1.2441125129491277
Line 104 tripped at t=1.2971048969582084
Line 35 tripped at t=1.463662221800599
Line 19 tripped at t=1.631494783605627
Line 21 tripped at t=1.713719828515491
Line 22 tripped at t=1.714505711798458
Line 3 tripped at t=1.7171793009352656
Line 4 tripped at t=1.7196344838095585
Line 5 tripped at t=1.719847148868116
Line 2 tripped at t=1.7221897352577358
Line 65 tripped at t=2.2485824909972756
Line 62 tripped at t=2.284167845716608
Line 63 tripped at t=2.3133969817263895
Line 23 tripped at t=2.4563202737365555
Line 40 tripped at t=2.7463364057801196
┌ Warning: Did not reach steady state!
=#


## different solver
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
damping = 0.1u"s"
scale_inertia = 1.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

sol = simulate_new_ND(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    # x_static=x_static,
    verbose=true,
    initial_fail=27,
    failtime=0.1,
    tspan=(0., 3.),
    trip_lines = :dynamic,
    trip_nodes = :none,
    freq_bound = 1.5,
    terminate_steady_state=true,
    ODE_solver = AutoTsit5(Rodas5P()),
    solverargs = (;reltol=1e-4),
    warn=true
    );

##### delete soon


nw_state = NWState(nw);
x0 = solve(SteadyStateProblem(nw, uflat(nw_state), pflat(nw_state)), NLSolveJL());
x_static = x0.u;


nw_state = NWState(nw)
s0 = NWState(nw, x_static, pflat(nw_state))


p = NWState(nw; ufill=0).p
p = NWState(nw).p
pflat(NWState(nw; ufill=0).p)
NWParameter(nw, p)


pflat(NWState(nw)) == pflat(NWParameter(nw)) # true
NWState(nw, x_static, pflat(NWState(nw))) == NWState(nw, x_static, pflat(NWParameter(nw))) # false

s1 = NWState(nw, x_static, pflat(NWState(nw)))
s2 = NWState(nw, x_static, pflat(NWParameter(nw)))

uflat(s1) == uflat(s2) && pflat(s1) == pflat(s2)

@which NWState(nw, x_static, pflat(NWState(nw))) == NWState(nw, x_static, pflat(NWParameter(nw)))


##
nw_state = NWState(nw)
p = pflat(nw_state)
x0 = solve(SteadyStateProblem(nw, uflat(nw_state), p), NLSolveJL())
x_static = x0.u

θidx = map(idx -> idx.compidx, vidxs(nw_state, :, "θ"))
s0 = NWState(nw_state, x_static, p)
zeroidx = 1;
offset = s0.v[θidx[zeroidx],:θ];
s0.v[θidx,:θ] .-= offset
s0
x_static = uflat(s0)
@assert iszero(s0.v[θidx[zeroidx],:θ])

s0.v[θidx,:θ] 
  0.0
 -0.000785761136132046
 -0.0005388966093823333
 -0.03880242157325059
 -0.039641903178045834
 -0.07589221567429288

s = NWState(nw)
p = pflat(s)
x0 = solve(SteadyStateProblem(nw, uflat(s), p), NLSolveJL()) # `uflat(s)` creates vector of zeros
#= `s0` is needed in order to access the symbolic vertex indices that have a θ-state.
In this model all vertices have a θ-state, however `uflat(s0)` returns the states ordered by 
the different vertex models/components. So here one cannot modify `x0.u` directly. =#
s0 = NWState(s, x0.u, p) # `NWState` object with steady state
x_static = uflat(s0)

dx = similar(x_static)
nw(dx, x_static, p, 0.0) # Rate of change (derivative) of each state variable storing it in `dx`.
maximum(abs.(dx))


#####

#= 
Shutdown line 27 at t=0.1
Line 29 tripped at t=0.193410177423005
Line 11 tripped at t=0.30445307930191867
Line 12 tripped at t=0.44383125504138254
Line 37 tripped at t=0.48805828010108876
Line 71 tripped at t=0.796308935388333
Line 56 tripped at t=0.9961123058421334
Line 59 tripped at t=0.9967495191082278
Line 58 tripped at t=1.015212693029165
Line 41 tripped at t=1.1016884928582986
Line 49 tripped at t=1.1036957953499746
Line 43 tripped at t=1.1067764938921625
Line 42 tripped at t=1.1085463045429322
Line 44 tripped at t=1.1128129636774686
Line 57 tripped at t=1.1185209337918582
Line 100 tripped at t=1.1941321492806816
Line 97 tripped at t=1.2390811075702792
Line 98 tripped at t=1.243197568029987
Line 73 tripped at t=1.244111946812969
Line 104 tripped at t=1.2971046528823473
Line 35 tripped at t=1.4636521959157789
Line 19 tripped at t=1.6314837991573508
Line 21 tripped at t=1.7137054781429202
Line 22 tripped at t=1.714508954854923
Line 3 tripped at t=1.717176165015607
Line 4 tripped at t=1.719642085100675
Line 5 tripped at t=1.7198389928225513
Line 2 tripped at t=1.722186040906692
Line 65 tripped at t=2.248531849634343
Line 62 tripped at t=2.284125004892754
Line 63 tripped at t=2.313321006329557
Line 23 tripped at t=2.456214087970374
Line 40 tripped at t=2.7462421180650254
┌ Warning: Did not reach steady state!
└ @ Main ~/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl:417
=#

# another stiff parameter setting
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
damping = 0.1u"s"
scale_inertia = 0.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

sol = simulate_new_ND(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    verbose=true,
    initial_fail=27,
    failtime=0.1,
    tspan=(0., 3.),
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    freq_bound = 0.42,
    terminate_steady_state=true,
    # solverargs = (;dtmax=0.01),
    #= # NOTE: with `AutoTsit5(Rodas5P())` the order and number of failures doesnt change with increased precision, 
    i.e. reltol=1e-4 instead of default value reltol=1e-3
    =#
    solverargs = (;reltol=1e-4), 
    warn=true
    );

## include node failures
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
damping = 0.1u"s"
scale_inertia = 1.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)
 
sol = simulate_new_ND(network;
    graph=network.graph, 
    gen_model=SwingDynLoadModel,
    verbose=true,
    initial_fail=27,
    failtime=0.1,
    tspan=(0., 10.),
    trip_lines = :dynamic,
    trip_nodes = :dynamic,
    freq_bound = 0.42,
    terminate_steady_state=true,
    solverargs = (;dtmax=0.01),
    warn=true
    );
#=
Shutdown line 27 at t=0.1
Line 29 tripped at t=0.19341026750998908
Vertex 16 tripped at t=0.20180015464622822
Vertex 18 tripped at t=0.22239078524106
Vertex 15 tripped at t=0.240791597245135
Vertex 14 tripped at t=0.24682570450857158
Line 35 tripped at t=0.247901506834905
Line 11 tripped at t=0.2499682340181457
Line 19 tripped at t=0.2546956849834631
Line 2 tripped at t=0.25700020384263866
Line 6 tripped at t=0.2582471937991574
Line 3 tripped at t=0.259470994911436
Line 20 tripped at t=0.26337072496770986
Line 1 tripped at t=0.3076831557926869
Vertex 7 tripped at t=0.32084364155892076
Vertex 1 tripped at t=0.8900273983023308
=#

###
### RTS: Test new ND w/o node failures against old ND w/o node failures (there was a bug in node failure CB)
###

## new ND
include(abspath(@__DIR__, "..", "..", "cluster/experiment_jarray/helpers_jarray.jl"))
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
damping = 0.1u"s"
scale_inertia = 1.2
tconst = 1.0u"s"
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

# get network that was created with old network to make sure to use the same network across the 2 versions
filepath = "/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND.mg"
network = loadgraph(filepath, MGFormat())

# steady_state_dict  = CSV.File("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_test_graph_old_ND_steady_state.csv")
# x_static = steady_state_dict[:SteadyState]

sol = simulate_new_ND(network;
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
include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
task_id = 1728
initial_fail = 78

sol = simulate_new_ND(exp_name_date, task_id, initial_fail;
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


###
### API collection
### 
using SymbolicIndexingInterface

i = 1 # vertex index
vidxs(nw, i, "ω") # returns symbolic index of vertex IF it has a state ω
vidxs(nw, i, :ω)  # returns symbolic index of vertex REGARDLESS if it has a ω state or not

vidxs(1:10, :ω) 

#= Returns index of :ω internal STATE of vertex i. Here one needs to keep in mind that the states are NOT ordered by
the vertices but by the different vertex models. For the RTS this means first the SwingDynLoad nodes then
the DynLoad nodes. =#
SII.variable_index(nw, VIndex(i, :ω)) 

s0.v[i,:ω] # returns ω-value of vertex i


# Get index instead of symbolic index. Sollte man nicht in hot loops verwenden wg. performance
map(idx -> idx.compidx, vidxs(nw, :, "ω"))

vidxs(nw, :, "θ")
vidxs(nw, :, "ω")
vpidxs(nw, :, "Pmech")

vidxs(nw, 13, "ω")[1].compidx # 13
vidxs(nw, 13, "ω")[1].subidx # ω

p.v[map(idx -> idx.compidx, vidxs(nw, :, "ω")), :Pmech]


using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
ω_idxs = Int[] 
for i in 1:nv(nw)
    try
        idx = SII.variable_index(nw, VIndex(i, :ω))
        push!(ω_idxs, idx)
    catch
        nothing
    end
end
ω_idxs




############################################################################################################
###########################################################################################################
############################################################################################################
# tryout

include("/home/brandner/.julia/dev/DynamicCascades/src/ND_model_new_ND.jl")
damping = 0.1u"s"
scale_inertia = 1.1 # NOTE this parameter has changed
tconst = 0.01u"s" # NOTE this parameter has changed
network = import_system(:rtsgmlc; damping, scale_inertia, tconst = tconst)

graph=network.graph
gen_model=SwingDynLoadModel
verbose=true
initial_fail=1
failtime=0.1
tspan=(0., 10.)
trip_lines = :dynamic
trip_nodes = :dynamic
freq_bound = 1.5
terminate_steady_state=true
solverargs=(;)


#= (#HACK) For the RTS network this adds missing parametes to `network`. 
For the WS network it does not add node parameters and it adds :_K 
and few line parameters that are not needed.
# TODO maybe put this in import_rtsgmlc.jl Update: Probably not possible with current setup=#
set_inertia!(network)
set_admittance!(network)

###
### build NetworkDynamics.jl Network
###

# set node failure mode
trip_nodes == :dynamic ? ωmax = freq_bound*2π : ωmax = Inf # This is a global node property

# loop over vertices and assign vertex models & parameter
vm_array = VertexModel[]
for i in 1:nv(network)
    # P = P_inj - P_load see `balance_power!`; P_inj = Pmech
    # P_inj = P + P_load
    P = ustrip(u"pu", get_prop(network, i, :P))
    Pload = ustrip(u"pu", get_prop(network, i, :P_load))
    Pmech = P + Pload
    τ = ustrip(u"s", get_prop(network, i, :timeconst))
    type = get_prop(network, i, :type)
    if type == :gen || type == :syncon
        M = ustrip(u"s^2", get_prop(network, i, :_M))
        γ = ustrip(u"s", get_prop(network, i, :damping))
        vm = gen_model(M=M,D=γ,τ=τ,ωmax=ωmax,Pmech=Pmech,Pload=Pload)
    elseif type == :load
        vm = DynLoadModel(τ=τ,Pload=Pload) 
    else
        error("The node type :$type is not defined!")
    end
    # set positions for plotting
    set_position!(vm, get_prop(network, i, :pos))
    push!(vm_array, vm)
end

# loop over edges and assign edge parameters
em_array = EdgeModel[]
for e in edges(network)
    srcV = ustrip(u"pu", get_prop(network, e.src, :Vm))
    dstV = ustrip(u"pu", get_prop(network, e.dst, :Vm))
    #= #NOTE `abs()` is used here because originally (old ND-version) `:_K` is negative 
    and the power flow is defined as `e[1] = - K * sin(v_s[1] - v_d[1])` see 
    org/DynamicCascades.jl/Why Coupling K is negative. =#
    Y = get_prop(network, e, :_Y)
    K = srcV * dstV * imag(Y) 
    Yabs = abs(Y)
    # K = abs(get_prop(network, e, :_Y))
    # set line failure mode
    trip_lines == :dynamic ? rating = ustrip(u"pu", get_prop(network, e, :rating)) : rating = Inf 
    em = Line(srcV=srcV, dstV=dstV, K=K, Yabs=Yabs, rating=rating)
    push!(em_array, em)
end

# generate `Network` object
# nw = Network(graph, vm_array, Line(K=K,rating=rating); dealias=true)
nw = Network(graph, vm_array, em_array)

# Check if network is power balanced (this is also checked before creating `nw` when importing the MetaGraph network in import_rtsgmlc.jl)
p = NWParameter(nw)
@assert isapprox(
    sum(p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pload")), :Pload]),
    sum(p.v[map(idx -> idx.compidx, vpidxs(nw, :, "Pmech")), :Pmech]))

# set initial perturbation CB
init_perturb = PresetTimeComponentCallback(failtime,
    ComponentAffect([], [:active]) do u, p, ctx
        println("Shutdown line $(ctx.eidx) at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
)
set_callback!(nw[EIndex(initial_fail)], init_perturb)
# NOTE
#= A loop would be syntactically possible but is not the way to go as `TerminateSelectiveSteadyState_new_ND`
should not be applied componentwise. =#
# for i in 1:nv(network)
#     add_callback!(nw[VIndex(i)], TerminateSelectiveSteadyState_new_ND(nw))
# end

###
### solve ODE problem
###

zeroidx = 1 
s0 = NWState(nw, steadystate_new_ND(nw; verbose=true, zeroidx=zeroidx), pflat(nw_state))
θidx = map(idx -> idx.compidx, vidxs(nw_state, :, "θ"))
x_static_new = s0.v[θidx,:θ]

steady_state_dict  = CSV.File("/home/brandner/.julia/dev/for_testing_delete_soon/RTS_new_stready_state_func_old_ND.csv")
x_static_old = steady_state_dict[:SteadyState]

@assert isapprox(x_static_new, x_static_old, atol=1e-8)


# TODO I-wo in neuen Code
# NOTE Don't use with steady states from old ND



