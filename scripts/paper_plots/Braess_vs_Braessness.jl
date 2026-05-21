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
include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"));

###
### preprocess and load data
###

######
###### WS
######

exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
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
model = SwingDynLoadModel_change_Pmech_only


### 
### Braessness vs Braessness
###

## full failure
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I; exclude_non_triggering=exclude_non_triggering);
exp_nr_full = exp_name_date[11:12]

## inertia failure
exp_name_date = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I; exclude_non_triggering=exclude_non_triggering);
exp_nr_BH = exp_name_date[11:12]

## power failure
exp_name_date = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_b_PIK_HPC_K_=3,N_G=32_20250718_160303.054"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I; exclude_non_triggering=exclude_non_triggering);
exp_nr_Pmech = exp_name_date[11:12]

prefix = exp_name_date[1:10]
exp_nrs = "$prefix,$exp_nr_full,$exp_nr_BH,$exp_nr_Pmech"

###
### 2D-Plots with colorscaling for z-Axis leaving out scatterpoints
###

# data
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes
failures = "line+node"

### 
### Braessness vs Braessness: Scatter using colorscaling for counts
###


fontsize = 34
markersize = 12
axislabelsize = fontsize

### 
### lines and nodes
###
failures = "line & node"
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes

################################################################################
########################## B^Power_i vs. B^{Inertia}_i #########################
################################################################################
fig = Figure(size=(1000,500),fontsize=fontsize)
fig[1,1] = ax = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B^{\text{ Inertia}}_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, ys)
sc = scatter!(ax, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc; label = "log₁₀(counts)", width = 30, spinewidth = 0)

ax.xticks = [-50, -25, 0, 25, 50]
offset = (20, -10)
text!(ax, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_vs_Bpower,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_vs_Bpower,exclude_non_triggering=$exclude_non_triggering.png"),fig)
# fig


################################################################################
### B_i vs. B^{Inertia}_i vs. B^{Power}_i in 2D ####
################################################################################


fig = Figure(size=(1000,500),fontsize=fontsize)
# No scatterpoints left out
fig[1,1] = ax = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B^{\text{ Inertia}}_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
color, colormap, colorrange, strokewidth = colorswitcher(zs; fancy_colors=true)
sc11 = scatter!(ax, xs, ys; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokecolor=:black, strokewidth=strokewidth)
Colorbar(fig[1,2], sc11; label = L"B_i", width = 30)

ax.xticks = [-60, -40, -20, 0, 20, 40, 60, 80]
ax.yticks = [-60, -40, -20, 0, 20, 40, 60, 80]
text!(ax, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_power_vs_inertia_vs_full_2D,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_vs_Bpower,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_power_vs_inertia_vs_full_2D,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_vs_Bpower,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig


################################################################################
####################### "B^{Inertia},B^{Power}_i vs B_i, #######################
################################################################################
fig = Figure(size=(1000,1000),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel=L"$B^{\text{ Inertia}}_i$",ylabel=L"Braessness $B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize,
    title = "Watts-Strogatz networks",
    titlefont = :regular,
    )
logc = get_log_count(ys, zs)
sc11 = scatter!(ax11, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
# Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)

fig[2,1] = ax21 = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"Braessness $B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, zs)
sc21 = scatter!(ax21, xs, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1:2,2], sc21; label = L"\log_{10}(\textrm{counts})", width = 30, spinewidth = 0)

linkxaxes!(ax11, ax21)
linkyaxes!(ax11, ax21)

ax21.xticks = [-50, -25, 0, 25, 50]
# Force upper axes to show no numbers, but keep tick marks
ticks = ax21.xticks[]
ax11.xticks = (ticks, fill("", length(ticks)))

offset = (20, -10)
text!(ax11, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)
text!(ax21, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+labellettersize, font = :bold, color = :black)

CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_Bpower_vs_Bi,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,Binertia_Bpower_vs_Bi,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig


### 
### lines
###
failures = "line"
xs = braessness_lines_Pmech
ys = braessness_lines_BH
zs = braessness_lines

fig = Figure(size=(1000,1500),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; title="Line failures", xlabel=L"$B^{\text{ Inertia}}_i$",ylabel=L"$B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize,
    )
logc = get_log_count(ys, zs)
sc11 = scatter!(ax11, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))

fig[2,1] = ax21 = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, zs)
sc21 = scatter!(ax21, xs, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1:3,2], sc21; label = "log₁₀(counts)", width = 30, height = Relative(2/3), spinewidth = 0)

fig[3,1] = ax31 = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B^{\text{ Inertia}}_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, ys)
sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))

linkxaxes!(ax11, ax21, ax31)
linkyaxes!(ax11, ax21, ax31)

ax31.xticks = [-10, 0, 10, 20, 30, 40]
# Force upper axes to show no numbers, but keep tick marks
ticks = ax31.xticks[]
ax11.xticks = (ticks, fill("", length(ticks)))
ax21.xticks = (ticks, fill("", length(ticks)))

offset = (20, -10)
text!(ax11, 0, 1; text = "b", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)
text!(ax21, 0, 1; text = "d", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)
text!(ax31, 0, 1; text = "f", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)
CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig

### 
### nodes
###
failures = "node"
xs = braessness_nodes_Pmech
ys = braessness_nodes_BH
zs = braessness_nodes

fig = Figure(size=(1000,1500),fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; title="Node failures", xlabel=L"$B^{\text{ Inertia}}_i$",ylabel=L"$B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize,
    )
logc = get_log_count(ys, zs)
sc11 = scatter!(ax11, ys, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))

fig[2,1] = ax21 = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, zs)
sc21 = scatter!(ax21, xs, zs; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1:3,2], sc21; label = "log₁₀(counts)", width = 30, height = Relative(2/3), spinewidth = 0)

fig[3,1] = ax31 = Axis(fig; xlabel=L"$B^{\text{ Power}}_i$",ylabel=L"$B^{\text{ Inertia}}_i$", xlabelsize=axislabelsize, ylabelsize=axislabelsize)
logc = get_log_count(xs, ys)
sc31 = scatter!(ax31, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))

linkxaxes!(ax11, ax21, ax31)
linkyaxes!(ax11, ax21, ax31)

ax31.xticks = [-50, -40, -30, -20, -10, 0, 10, 20]
# Force upper axes to show no numbers, but keep tick marks
ticks = ax31.xticks[]
ax11.xticks = (ticks, fill("", length(ticks)))
ax21.xticks = (ticks, fill("", length(ticks)))

offset = (20, -10)
text!(ax11, 0, 1; text = "a", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)
text!(ax21, 0, 1; text = "c", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)
text!(ax31, 0, 1; text = "e", space = :relative, align = (:left, :top), offset = offset, fontsize = fontsize+20, font = :bold, color = :black)

CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,exclude_non_triggering=$exclude_non_triggering.png"),fig)
fig

