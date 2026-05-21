"""
Plotting scripts: Braessness of individual lines.
"""

using DynamicCascades
 

using CairoMakie
using Statistics
CairoMakie.activate!()
using Serialization


################################################################################
############################## Plotting scripts ################################
################################################################################
include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"))

###
### preprocess and load data
###

######
###### RTS
######
exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)

# choose inertia
I = 3 # scaling factor
freq_bounds = [0.01, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30,
    0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54,
    0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74,
    0.8, 0.85, 0.90, 1.00, 1.2, 1.4, 1.6, 1.8, 2.0]

# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)

# choose `f_b_narrow` and `f_b_wide`
f_b_narrow = 0.3
f_b_wide = 1.0
exclude_non_triggering = true
N_nodes = 73 # nv(import_system(:rtsgmlc))

###
### Braessness histogram
###

fontsize = 24
titlesize = (fontsize-5)
markersize = 8

braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I; exclude_non_triggering=exclude_non_triggering)


fig = Figure(size = (700,1300), fontsize = fontsize)
yscale = Makie.pseudolog10  # maps 0→0, then log10(count+1)


# bins
all_breaessnesses = vcat(braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes)
max_braessness = maximum(all_breaessnesses)
min_braessness = minimum(all_breaessnesses)
bins= (min_braessness-0.5):1:(max_braessness+0.5)

## node failures
fig[1, 1] = ax11 = Axis(fig; yscale=yscale,
    title = "Power grid",
    # titlefont = :regular,
    # xlabel = L"Braessness node failures $B^n_i$",
    # ylabel = L"\textrm{Count}"
    )

μ = round(mean(braessness_nodes), digits=3) # 2.133
# hist!(ax11, braessness_nodes; label=L"$\left< B^n \right> = 0.34$", color="#8C8C8CFF", bins=bins, strokewidth=0.2)
hist!(ax11, braessness_nodes; label=L"Node failures, $\left< B \right>=2.13$", color=:red, bins=bins, strokewidth=0.2)
axislegend(ax11; position = :rt)

## line failures
fig[2, 1] = ax21 = Axis(fig; yscale=yscale,
    # title = "Watts-Strogatz networks",
    # titlefont = :regular,
    # xlabel = L"Braessness line failures $B^l_i$",
    ylabel = L"\textrm{Count}")

μ = round(mean(braessness_lines), digits=3) # 4.40
hist!(ax21, braessness_lines; label=L"Line failures, $\left< B \right>=4.40$", color=:orange, bins=bins, strokewidth=0.2)
axislegend(ax21; position = :rt)

## node + line failures
fig[3, 1] = ax31 = Axis(fig; yscale=yscale,
    # title = "Watts-Strogatz networks",
    # titlefont = :regular,
    xlabel = L"Braessness $B_i$",
    # ylabel = L"\textrm{Count}"
    )

μ = round(mean(braessness_lines .+ braessness_nodes), digits=3) # 6.53
hist!(ax31, braessness_lines .+ braessness_nodes; label=L"Node and line failures, $\left< B \right> = 6.53$", color="#8C8C8CFF", bins=bins, strokewidth=0.2)
axislegend(ax31; position = :rt)
axislegend(ax31; position = :rt)
ax11.yticks = ax21.yticks = ax31.yticks = [0, 1, 3, 5, 8, 1000, 5000]
ax11.xticks = ax21.xticks = ax31.xticks = [-10, 0, 10, 20, 30, 40, 50]
ax11.xticklabelsvisible = ax21.xticklabelsvisible = false

# link the marginal histograms to the scatter
linkxaxes!(ax11, ax21, ax31); linkyaxes!(ax11, ax21, ax31)
Label(fig[1, 1, TopLeft()], fontsize=fontsize+labellettersize, "b", font = :bold, padding = (-55, 0, 5, 0))
Label(fig[2, 1, TopLeft()], fontsize=fontsize+labellettersize, "d", font = :bold, padding = (-55, 0, 5, 0))
Label(fig[3, 1, TopLeft()], fontsize=fontsize+labellettersize, "f", font = :bold, padding = (-55, 0, 5, 0))

model = "full_failure"
CairoMakie.save(joinpath(MA_DIR, "histograms_braessness_RTS_lines,nodes_separately_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "histograms_braessness_RTS_lines,nodes_separately_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig

