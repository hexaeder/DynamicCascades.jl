"""
Plotting scripts: Braessness of individual lines.
"""

using DynamicCascades
using DynamicCascades: PLOT_DIR

using CairoMakie
using Statistics
CairoMakie.activate!()
using Serialization


################################################################################
############################## Plotting scripts ################################
################################################################################

###
### preprocess and load data
###

######
###### WS
######

# full failure
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
model = "full_failure"

exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

# choose inertia
I = 7.5
freq_bounds = [0.03, 0.15]
# freq_bounds = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.22, 0.25, 0.3, 0.8]

# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date) # NOTE For WS set `res_tol=1e-5,` in `steadystate`

# choose `f_b_narrow` and `f_b_wide`
f_b_narrow = 0.03
f_b_wide = 0.15
exclude_non_triggering = true
N_nodes = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:N_nodes]

###
### Braessness histogram
###

fontsize = 28
titlesize = (fontsize-5)
markersize = 8

braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I; exclude_non_triggering=exclude_non_triggering)


fig = Figure(size = (700,500), fontsize = fontsize)
yscale = Makie.pseudolog10  # maps 0→0, then log10(count+1)
fig[1, 1] = ax = Axis(fig; yscale=yscale, 
    xlabel = L"Braessness $B_i$",
    ylabel = L"Count")
# fig[1, 1] = ax = Axis(fig; yscale=yscale, xlabel = L"\mathsf{Frequency\ bound} Braessness $B_i$", ylabel = "Count")


# bins
all_breaessnesses = vcat(braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes)
max_braessness = maximum(all_breaessnesses)
min_braessness = minimum(all_breaessnesses)
bins= (min_braessness-0.5):1:(max_braessness+0.5)

μ = round(mean(braessness_lines .+ braessness_nodes), digits=3)
# hist!(ax, braessness_lines .+ braessness_nodes; label=L"$\left< B \right> = 0.524$", color="#8C8C8CFF", bins=bins, strokewidth=0.2)
# hist!(ax, braessness_lines .+ braessness_nodes; label=L"$\left< B \right> = 0.524$", color="#8C8C8CFF", bins=bins, strokewidth=0.2)
# hist!(ax, braessness_lines .+ braessness_nodes; label="⟨B⟩ = 0.524", color="#8C8C8CFF", bins=bins, strokewidth=0.2)
hist!(ax, braessness_lines .+ braessness_nodes; label=L"$\left< B \right> = 0.524$", color="#8C8C8CFF", bins=bins, strokewidth=0.2)

axislegend(ax; position = :rt, labelsize=fontsize)
ax.yticks = [0, 1, 10, 100, 1000, 5000]

CairoMakie.save(joinpath(MA_DIR, "braessness", "histograms_braessness_lines_and_nodes_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "braessness", "histograms_braessness_lines_and_nodes_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig

