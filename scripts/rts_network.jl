#!/usr/bin/env julia

using OrdinaryDiffEq
using DiffEqCallbacks
using SteadyStateDiffEq
using DataFrames

using DynamicCascades

@info "Find steady state..."
network = import_system(:rts96prepared; gen_γ=1.0, slack_γ=1.0, load_τ=0.1, losses=false)
(nd, p) = nd_model(network);

x0 = zeros(length(nd.syms));
x_static = solve(SteadyStateProblem(nd, x0, p), SSRootfind())
is_static_state(nd, x_static, p)
# x_static = solve(SteadyStateProblem(nd, x0, p), DynamicSS(AutoTsit5(Rosenbrock23())))

# tspan = (0., 2000.);
# prob = ODEProblem(nd, x0, tspan, p);
# @time solinit = solve(prob, Rosenbrock23(), callback=TerminateSteadyState(), save_everystep=false);

# if solinit.t[end] == tspan[2]
#     @warn "Simulation time ended without reaching steady state!"
# else
#     @info "Found steadystate at t = $(solinit.t[end])..."
# end
# x_static = copy(solinit[end]);

initial_fail = [27]
failtime = 1.0

(nd, p, overload_cb) = nd_model(network);
tspan = (0., 100.)
prob = ODEProblem(nd, copy(x_static), tspan, p);

line_failures = SavedValues(Float64, Int);
S_values = SavedValues(Float64, Vector{Float64});
P_values = SavedValues(Float64, Vector{Float64});
cbs = CallbackSet(initial_fail_cb(initial_fail, failtime),
                  overload_cb(trip_lines=true, load_S=S_values, load_P=P_values, failures=line_failures));

@time sol = solve(prob, AutoTsit5(Rosenbrock23()), callback=cbs, progress=true);

inspect_solution(network, nd, sol, S_values, P_values)
