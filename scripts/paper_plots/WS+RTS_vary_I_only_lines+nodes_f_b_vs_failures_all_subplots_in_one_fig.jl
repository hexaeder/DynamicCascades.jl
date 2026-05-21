"""
Watts-Strogatz-Network-Ensemble: Using job array framework. Transition that appears
when varying the frequency bounds. Line and node failures summed.
"""

include(abspath(@__DIR__, "..", "helpers_jarray.jl"))

using GraphMakie
using Colors, ColorSchemes
using CairoMakie
CairoMakie.activate!()

include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"))

# general plotting parameters
create_posprocessing_data = false # set to `false` for fast plotting
opacity = 0.30
# markers
marker = (:circle, ":circle")
markersize = 12
linewidth = 4
patchsize = (60, 30)

fontsize = 50
labellettersize = 4

fig = Figure(size=(1600,1200),fontsize = fontsize)

################################################################################
############################### Watts-Strogatz #################################
################################################################################

# exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
left_out_β_values = []
i = 30

###################### Calculate mean and standard error #######################
if create_posprocessing_data == true
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

################################## Plotting ####################################

freq_bounds = exp_params_dict[:freq_bounds]
filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))


# fig_WS_summed
sum_lines_nodes = true
normalize = false

ax_WS_summed = Axis(fig[1, 1], xticklabelrotation=π/2,
title = "Watts-Strogatz networks",
titlefont = :regular,
ylabel = L"# Failures $\left< F \right>$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_WS_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)



# fig_WS_separate
sum_lines_nodes = false
normalize = true

ax_WS_separate = Axis(fig[2, 1], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
ylabel = normalize ? "Fractions failing elements" : L"Averaged failures $N_{fail}$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_WS_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)


# ax_WS_summed.xticks = filtered_freq_bounds
# ax_WS_summed.xlabelpadding = 15
axislegend(ax_WS_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)

# ax_WS_separate.xticks = filtered_freq_bounds[1:2:end]
ax_WS_separate.xticks =  ax_WS_summed.xticks = [0.005, 0.025, 0.045, 0.065, 0.085, 0.110, 0.130, 0.150]

ax_WS_summed.xlabelpadding = 15
# lines!(ax_WS_summed, [NaN], [NaN]; label=L"Damping $D=1$ [$s$]", color=:white)
# axislegend(ax_WS_summed, position = :rt, labelsize=labelsize)


################################################################################
##################################### RTS ######################################
################################################################################

exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
# left_out_frequencies = [0.01, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.26, 0.28, 0.30,
#     0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.51, 0.52, 0.53, 0.54, 0.55,
#     0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.72, 0.74,
#     0.8, 0.85, 0.90, 1.00, 1.2, 1.4, 1.6, 1.8, 2.0] # combined runs
left_out_frequencies = [0.01, 0.08, 1.8, 2.0]
i = 8

###################### Calculate mean and standard error #######################
if create_posprocessing_data == true
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

################################## Plotting ####################################

freq_bounds = exp_params_dict[:freq_bounds]
filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))


# fig_RTS_summed
sum_lines_nodes = true
normalize = false

ax_RTS_summed = Axis(fig[1, 2], xticklabelrotation=π/2,
# title = sum_lines_nodes ? "Summed line and node failures" : "Line and node failures",
title = "Power grid",
titlefont = :regular,

# ylabel = normalize ? "Normalized average of failures" : L"Averaged failures $N_{fail}$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_RTS_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

# fig_RTS_separate
sum_lines_nodes = false
normalize = true

ax_RTS_separate = Axis(fig[2, 2], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
# ylabel = normalize ? "Normalized average of failures" : L"Averaged failures $N_{fail}$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_RTS_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)


ax_RTS_separate.xlabelpadding = 15
# lines!(ax_RTS_separate, [NaN], [NaN]; label=L"Damping $D=0.1$ [$s$]", color=:white)
axislegend(ax_RTS_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)


ax_RTS_separate.xticks = ax_RTS_summed.xticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]

ax_RTS_separate.xlabelpadding = 15
# lines!(ax_RTS_summed, [NaN], [NaN]; label=L"Damping $D=0.1$ [$s$]", color=:white)
# axislegend(ax_RTS_summed, position = :rt, labelsize=labelsize)


# link x-axes
linkxaxes!(ax_WS_separate,ax_WS_summed)
linkxaxes!(ax_RTS_separate,ax_RTS_summed)
ax_WS_summed.xticklabelsvisible = false; ax_RTS_summed.xticklabelsvisible = false; 

# Label(fig[1, 1, TopLeft()], "a", fontsize = 44, font = :bold, padding = (5, 5, 5, 5))
# Label(fig[1, 2, TopLeft()], "b", fontsize = 44, font = :bold, padding = (5, 5, 5, 5))
# Label(fig[2, 1, TopLeft()], "c", fontsize = 44, font = :bold, padding = (5, 5, 5, 5))
# Label(fig[2, 2, TopLeft()], "d", fontsize = 44, font = :bold, padding = (5, 5, 5, 5))

offset = (80, -20)
text!(ax_WS_summed, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_RTS_summed, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_WS_separate, 0, 1; text = "c", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_RTS_separate, 0, 1; text = "d", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)


CairoMakie.save(joinpath(MA_DIR, "WS+RTS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "WS+RTS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig.png"),fig)
fig
