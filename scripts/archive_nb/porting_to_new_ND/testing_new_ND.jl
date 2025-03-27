###
### Fixpoint calculation for old ND and new ND
###

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
more recent version of DE.jl this throws an error to avoid an inconsisten use of `NonlinearProblem`,
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