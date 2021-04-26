using OrdinaryDiffEq
using DiffEqCallbacks
using BenchmarkTools
using Plots
using NetworkDynamics
using IEEE_RTS_96

rts96 = load_rts96()

(nd, p, load_cb, saving_cb) = nd_model(rts96);

x0 = zeros(rts96.num_of_buses + rts96.num_of_generators);
tspan = (0., 2000.)
prob = ODEProblem(nd, x0, tspan, p)

@time sol = solve(prob, Rosenbrock23());
@time sol = solve(prob, Rosenbrock23(), callback=cb);

plot(sol)

x_static = copy(sol[end])

# for line in 1:rts96.num_of_lines
    line = 27
    pnew = deepcopy(p)
    pnew[2][line] = 0.0
    tspan = (0., 50.)
    prob = ODEProblem(nd, x_static, tspan, pnew)
    saved_values = SavedValues(Float64, Vector{Float64})
    @time sol = solve(prob, Rosenbrock23(), callback=CallbackSet(load_cb, saving_cb(saved_values)));
# end


plot_apparent_power(saved_values, [29,83,37,26,10,84,11,6,2,12,48], limits=rts96, xlims=(0, 5))
plot_apparent_power(saved_values, [26], limits=rts96, xlims=(0, 5))

plot_apparent_power(saved_values, 1:100)

plot(sol, xlims=(0,10), ylims=(-2,2))
