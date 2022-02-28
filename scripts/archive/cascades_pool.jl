#!/usr/bin/env julia

@assert VERSION == v"1.6.0-rc2"
using Distributed
const ON_POOL = occursin("pool", gethostname())

@info "Initialize environment on main process"
using Pkg
PKG_DIR = ON_POOL ? "/data/scratch/wuerfel/ieee-rts-96_workdir" : joinpath(@__DIR__, "..")
Pkg.activate(PKG_DIR)
Pkg.instantiate()
Pkg.precompile()

@info "Say hello to my workers!"
if ON_POOL
    using SlurmClusterManager
    addprocs(SlurmManager())
else
    newprocs = 4-length(procs())
    addprocs(newprocs)
end
@everywhere println("Hello from worker $(myid()) @ $(gethostname())!")

@everywhere using Pkg
@everywhere Pkg.activate($PKG_DIR)
@everywhere begin
    using OrdinaryDiffEq
    using DiffEqCallbacks
    using NetworkDynamics
    using DataFrames
    using CSV
    using Dates: now
    using IEEE_RTS_96
end

@everywhere rts96 = load_rts96()

@info "Find steady state on main worker..."
(nd, p, load_cb, saving_cb) = nd_model(rts96);
x0 = zeros(rts96.num_of_buses + rts96.num_of_generators);
tspan = (0., 2000.);
prob = ODEProblem(nd, x0, tspan, p);
sol = solve(prob, Rosenbrock23(), callback=TerminateSteadyState());
tmax = sol.t[end]
@info "Found steadystate at t = $tmax..."

@everywhere x_static = copy($sol[end])

lines = 1:rts96.num_of_lines
gen_τs   = 0.1:0.1:0.5
slack_τs = 0.1:0.1:0.5
load_τs  = 0.1:0.1:0.5
hyperparameter = collect(Iterators.product(lines, gen_τs, slack_τs, load_τs))[:]

results = @sync @distributed vcat for hyperp in hyperparameter
    println("Worker $(myid()): $(hyperp)")

    (line, gen_τ, slack_τ, load_τ) = hyperp
    (nd, plocal, load_cb, saving_cb) = nd_model(rts96; gen_τ, slack_τ, load_τ);

    plocal[2][line] = 0.0
    tspan = (0., 100.)
    prob = ODEProblem(nd, copy(x_static), tspan, plocal)
    linebreaks = SavedValues(Float64, Int)
    sol = solve(prob, Rosenbrock23(), callback=load_cb(linebreaks));

    df = DataFrame()
    df[!, :initial_fault] = [line for i in 1:length(linebreaks.t)]
    df[!, :t] = linebreaks.t
    df[!, :fault] = linebreaks.saveval
    df[!, :gen_τ]   .= plocal[1][findfirst(isequal(1), rts96.bustype)][3]
    df[!, :slack_τ] .= plocal[1][findfirst(isequal(2), rts96.bustype)][3]
    df[!, :load_τ]  .= plocal[1][findfirst(isequal(3), rts96.bustype)][3]
    df[!, :theta] = sol.(df.t)
    df
end
results

dir = ON_POOL ? "/data/scratch/wuerfel/ieee-rts-96_data" : "/Users/hw/tmp/ieee-rts96_data"
isdir(dir) || mkdir(dir)

CSV.write(joinpath(dir, "$(now())_results.csv"), results)
