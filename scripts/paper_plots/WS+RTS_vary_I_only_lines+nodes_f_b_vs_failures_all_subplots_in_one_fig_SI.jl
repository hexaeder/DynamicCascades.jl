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
markersize = 15
linewidth = 5
patchsize = (60, 30)
offset = (80, -20)

fontsize = 50
labellettersize = 4



################################################################################
############################### Watts-Strogatz #################################
################################################################################

fig = Figure(size=(1600*1.2,1200*1.2),fontsize = fontsize)

# exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250126_012357.344"
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
left_out_frequencies = [
    0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.8]
left_out_β_values = []


###################### Calculate mean and standard error #######################
if create_posprocessing_data == false
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

################################## Plotting ####################################

freq_bounds = exp_params_dict[:freq_bounds]
filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))

#################################### i = 5 #####################################
i = 5

# fig_WS_summed
sum_lines_nodes = true
normalize = false

ax_summed = Axis(fig[1, 1], xticklabelrotation=π/2,
title= L"Inertia $I=5.0$ $s^2$ (WS)",
titlefont = :regular,
ylabel = L"# Failures $\left< F \right>$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

# fig_WS_separate
sum_lines_nodes = false
normalize = true

ax_separate = Axis(fig[2, 1], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
ylabel = normalize ? "Fractions failing elements" : L"Averaged failures $N_{fail}$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

axislegend(ax_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)
ax_summed.xticks = ax_separate.xticks = [0.005, 0.025, 0.045, 0.065, 0.085, 0.110, 0.130, 0.150]
ax_summed.xlabelpadding = 15

# link x-axes
linkxaxes!(ax_separate,ax_summed)
ax_summed.xticklabelsvisible = false

text!(ax_summed, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_separate, 0, 1; text = "c", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

#################################### i = 20 #####################################
i = 20

# fig_WS_summed
sum_lines_nodes = true
normalize = false

ax_summed = Axis(fig[1, 2], xticklabelrotation=π/2,
title= L"Inertia $I=20.0$ $s^2$ (WS)",
titlefont = :regular,
# ylabel = L"# Failures $\left< F \right>$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

# fig_WS_separate
sum_lines_nodes = false
normalize = true

ax_separate = Axis(fig[2, 2], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
# ylabel = normalize ? "Fractions failing elements" : L"Averaged failures $N_{fail}$",
)

WS_f_b_vs_failures!(exp_name_date, left_out_frequencies, left_out_β_values,  ax_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)


axislegend(ax_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)
ax_summed.xticks = ax_separate.xticks = [0.005, 0.025, 0.045, 0.065, 0.085, 0.110, 0.130, 0.150]
ax_summed.xlabelpadding = 15

# link x-axes
linkxaxes!(ax_separate,ax_summed)
ax_summed.xticklabelsvisible = false

text!(ax_summed, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_separate, 0, 1; text = "d", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)


CairoMakie.save(joinpath(MA_DIR, "WS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig,I=5,20.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "WS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig,I=5,20.png"),fig)
fig



################################################################################
##################################### RTS ######################################
################################################################################
fig = Figure(size=(1600*1.2,1200*1.2),fontsize = fontsize)

exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
exp_data_dir = joinpath(RESULTS_DIR, exp_name_date)
# left_out_frequencies = [0.01, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.26, 0.28, 0.30,
#     0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.49, 0.51, 0.52, 0.53, 0.54, 0.55,
#     0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.66, 0.68, 0.72, 0.74,
#     0.8, 0.85, 0.90, 1.00, 1.2, 1.4, 1.6, 1.8, 2.0] # combined runs
left_out_frequencies = [0.01, 0.08, 1.8, 2.0]

###################### Calculate mean and standard error #######################
if create_posprocessing_data == false
    postprocess_jarray_data(exp_name_date)
end
df_config = DataFrame(CSV.File(joinpath(exp_data_dir, "config.csv")))
exp_params_dict = Serialization.deserialize(joinpath(exp_data_dir, "exp.params"))

################################## Plotting ####################################

freq_bounds = exp_params_dict[:freq_bounds]
filtered_freq_bounds = filter!(x->x ∉ left_out_frequencies, deepcopy(freq_bounds))


#################################### i = 5 #####################################
i = 1

# fig_WS_summed
sum_lines_nodes = true
normalize = false

ax_summed = Axis(fig[1, 1], xticklabelrotation=π/2,
title= L"Inertia $I=1.0$ $s^2$ (PG)",
titlefont = :regular,
ylabel = L"# Failures $\left< F \right>$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

# fig_WS_separate
sum_lines_nodes = false
normalize = true

ax_separate = Axis(fig[2, 1], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
ylabel = normalize ? "Fractions failing elements" : L"Averaged failures $N_{fail}$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

axislegend(ax_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)
ax_summed.xticks = ax_separate.xticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
ax_summed.xlabelpadding = 15

# link x-axes
linkxaxes!(ax_separate,ax_summed)
ax_summed.xticklabelsvisible = false

text!(ax_summed, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_separate, 0, 1; text = "c", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

#################################### i = 20 #####################################
i = 5

# fig_WS_summed
sum_lines_nodes = true
normalize = false

ax_summed = Axis(fig[1, 2], xticklabelrotation=π/2,
title= L"Inertia $I=5.0$ $s^2$ (PG)",
titlefont = :regular,
# ylabel = L"# Failures $\left< F \right>$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_summed, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

# fig_WS_separate
sum_lines_nodes = false
normalize = true

ax_separate = Axis(fig[2, 2], xticklabelrotation=π/2,
xlabel = L"Frequency bound $f_b$ [Hz]",
# ylabel = normalize ? "Fractions failing elements" : L"Averaged failures $N_{fail}$",
)

RTS_f_b_vs_failures!(exp_name_date, left_out_frequencies,  ax_separate, i, sum_lines_nodes, normalize; opacity=opacity,linewidth=linewidth,markersize=markersize)

axislegend(ax_separate, position = :rt, labelsize=labelsize, patchsize=patchsize)
ax_summed.xticks = ax_separate.xticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
ax_summed.xlabelpadding = 15

# link x-axes
linkxaxes!(ax_separate,ax_summed)
ax_summed.xticklabelsvisible = false

text!(ax_summed, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax_separate, 0, 1; text = "d", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

CairoMakie.save(joinpath(MA_DIR, "RTS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig,I=1,5.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "RTS_vary_I_only_lines+nodes_f_b_vs_failures_all_subplots_in_one_fig,I=1,5.png"),fig)
fig










