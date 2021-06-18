using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using NetworkDynamics
using DataFrames

using DynamicCascades

@info "Find steady state..."
gen_τ   = 10.0
slack_τ = 1.0
load_τ  = 0.1

network = import_system(:kaiser2020)
network2 = import_system(:rts96prepared)
(nd, p) = nd_model(network; gen_τ, slack_τ, load_τ);

x0 = zeros(length(nd.syms));
x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())

# function project_theta(θ)
#     n = (θ + π) ÷ 2π
#     return θ - n * 2π
# end
# θidx = idx_containing(nd, "θ")
# x_static[θidx] = project_theta.(x_static[θidx])
# is_static_state(nd, x_static, p)

initial_fail = [55]
failtime = 1.0

(nd, p, overload_cb) = nd_model(network; gen_τ, slack_τ, load_τ);
tspan = (0., 100.)
prob = ODEProblem(nd, copy(x_static), tspan, p);

line_failures = SavedValues(Float64, Int);
S_values = SavedValues(Float64, Vector{Float64});
P_values = SavedValues(Float64, Vector{Float64});
cbs = CallbackSet(initial_fail_cb(initial_fail, failtime),
                  overload_cb(trip_lines=false, load_S=S_values, load_P=P_values, failures=line_failures));

@time sol = solve(prob, AutoTsit5(Rosenbrock23()), callback=cbs, progress=true);

inspect_solution(network, nd, sol, S_values, P_values)
