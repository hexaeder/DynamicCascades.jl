using IEEE_RTS_96
using CSV
using DataFrames
using Plots
using GraphPlot, Cairo, Compose
using LightGraphs
using Random: MersenneTwister

rts96 = load_rts96();
all = CSV.read("../results/2021-03-15T18:28:42.305_results.csv", DataFrame);

groups = groupby(all, [:gen_τ, :load_τ, :slack_τ]);

df = groups[1];


cascades = groupby(df, :initial_fault)

for c in cascades
    initial = c.initial_fault[1]
    p = cascadeplot(rts96, initial, c.t, c.fault);
    draw(PNG("cascade_fault_$initial.png", 26cm, 26cm), p)
end
