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
# exp_name_date = "WS_k=4_exp10_vary_I_only_lines_and_nodes_N=20_PIK_HPC_K_=3,N_G=16_20250410_185732.694"
exp_name_date = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_142151.671"
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

# f_b_narrow = 0.03
# f_b_wide = 0.035

N_nodes = Serialization.deserialize(joinpath(RESULTS_DIR, exp_name_date, "exp.params"))[:N_nodes]
model = SwingDynLoadModel_change_Pmech_only

###
### Braessness histograms
###
## Braessness vs rho
fig = plot_braessness_vs_rho_scatter_and_histograms(exp_name_date, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

## Braessness vs distance
fig = plot_braessness_vs_dist_histograms(exp_name_date, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


### 
### Braessness vs Braessness
###

## full failure
exp_name_date = "WS_k=4_exp04_vary_I_only_lines_and_nodes_PIK_HPC_K_=3,N_G=32_20250321_171511.976"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);
exp_nr_full = exp_name_date[11:12]

## inertia failure
exp_name_date = "WS_k=4_exp11_vary_I_only_lines_and_nodes_change_to_BH_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_145304.04"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_BH, x_dist_BH = get_network_measures(exp_name_date);
exp_nr_BH = exp_name_date[11:12]

## power failure
exp_name_date = "WS_k=4_exp12_vary_I_only_lines_and_nodes_change_Pmech_complement_f_bPIK_HPC_K_=3,N_G=32_20250420_142151.671"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_Pmech, x_dist_Pmech = get_network_measures(exp_name_date);
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

# xs = braessness_lines_Pmech
# ys = braessness_lines_BH
# zs = braessness_lines
# failures = "lines"

# xs = braessness_nodes_Pmech
# ys = braessness_nodes_BH
# zs = braessness_nodes
# failures = "nodes"

fig = plot_braessness_power_vs_inertia_vs_full_2D(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.png"),fig)


### 
### Braessness vs Braessness: Scatter using colorscaling for counts
###

# lines and nodes
failures = "line & node"
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

# lines 
failures = "line"
xs = braessness_lines_Pmech
ys = braessness_lines_BH
zs = braessness_lines
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_lines_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

# nodes
failures = "node"
xs = braessness_nodes_Pmech
ys = braessness_nodes_BH
zs = braessness_nodes
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


# plot power f. + inertia f. (sum braessnesses) vs full f.
failures = "line & node"
xs = braessness_lines_BH .+ braessness_nodes_BH .+ braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines .+ braessness_nodes

fontsize = 25
titlesize = (fontsize-5)
markersize = 10

fig = Figure(size=(1200,900), fontsize=fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure + power failure",ylabel="full failure", title="Braessness power f. + inertia f. (sum braessnesses) vs full f $failures failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)
fig
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_inertia_f_+_power_f_vs_full_f_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_inertia_f_+_power_f_vs_full_f_lines_and_nodes_logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


###
### 3D scatter plot counts
###
using GLMakie # for interactivity
GLMakie.activate!()

# lines and nodes
xs = braessness_lines_BH .+ braessness_nodes_BH
ys = braessness_lines_Pmech .+ braessness_nodes_Pmech
zs = braessness_lines .+ braessness_nodes

fig = plot_braessness_power_vs_inertia_vs_full_3D(xs, ys, zs, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

###
### 3D scatter plot only scatterpoints with count=1
###
freqs3D = Dict{NTuple{3,Int}, Int}()
for (x,y,z) in zip(xs, ys, zs)
    freqs3D[(x,y,z)] = get(freqs3D, (x,y,z), 0) + 1
end

fig = Figure(size=(1470,1050), fontsize=fontsize)
fig[1,1] = ax11 = Axis3(fig; xlabel="inertia failure", ylabel="power failure", zlabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)

xs_outlier = Int64[]; ys_outlier = Int64[]; zs_outlier = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if freqs3D[(x,y,z)] < 2
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
end
color, colormap, colorrange, strokewidth = colorswitcher(zs_outlier; fancy_colors)
sc11 = scatter!(ax11, xs_outlier, ys_outlier, zs_outlier; markersize=markersize, color=color, colormap=colormap, colorrange=colorrange, strokewidth=strokewidth)
Colorbar(fig[1,2], sc11; label = "full failure", width = 30)
fig
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_vs_braessness_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


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
N_nodes = 73 # nv(import_system(:rtsgmlc))


###
### Braessness histograms
###
## Braessness vs rho
fig = plot_braessness_vs_rho_scatter_and_histograms(exp_name_date, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "scatter_histograms_and_rho_braessness_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

## Braessness vs distance
fig = plot_braessness_vs_dist_histograms(exp_name_date, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "dist_braessness_histogram_$model,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)



###
### Braessness of single lines that interconnect large subgrids
###
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
braessness_lines_plus_nodes = braessness_lines .+ braessness_nodes
x_ρ, x_dist = get_network_measures(exp_name_date);

# TODO maybe check rating
# TODO check trajectories
## line 37
line = 37
edge = collect(edges(import_system(:rtsgmlc).graph))[line] # Edge 21(gen) => 73(load)
braessness_lines_plus_nodes[line] # 0
x_ρ[line] # 3.98
x_dist[line] # 17

## line 73
line = 73
edge = collect(edges(import_system(:rtsgmlc).graph))[line] # Edge 47(gen) => 66(gen)
braessness_lines_plus_nodes[line] # 0
x_ρ[line] # 1.00
x_dist[line] # 17 

## line 108
line = 108
edge = collect(edges(import_system(:rtsgmlc).graph))[line] # Edge 71(gen) => 73(load)
braessness_lines_plus_nodes[line] # 0
x_ρ[line] # -7.06
x_dist[line] # 17


"""
 - TODO Check rating and initial load of these lines. (Lines may be weekly loaded in initial state)
"""


### 
### Braessness vs Braessness
###

## full failure
exp_name_date = "RTS_exp04_variation_frequency+inertia_PIK_HPC_20250616_213442.721"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines, braessness_nodes = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ, x_dist = get_network_measures(exp_name_date);
exp_nr_full = exp_name_date[8:9]

## inertia failure
exp_name_date = "RTS_exp05_variation_frequency+inertia_inertia_failure_PIK_HPC_20250701_174359.744"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_BH, braessness_nodes_BH = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_BH, x_dist_BH = get_network_measures(exp_name_date);
exp_nr_BH = exp_name_date[8:9]

## power failure
exp_name_date = "RTS_exp06_variation_frequency+inertia_power_failure_PIK_HPC_20250701_175708.899"
# write_failures_f_b_to_df(exp_name_date, freq_bounds, I)
# write_network_measures_to_df(exp_name_date)
braessness_lines_Pmech, braessness_nodes_Pmech = get_braessness(exp_name_date, f_b_narrow, f_b_wide, I);
x_ρ_Pmech, x_dist_Pmech = get_network_measures(exp_name_date);
exp_nr_Pmech = exp_name_date[8:9]

prefix = exp_name_date[1:7]
exp_nrs = "$prefix,$exp_nr_full,$exp_nr_BH,$exp_nr_Pmech"

###
### 2D-Plots with colorscaling for z-Axis leaving out scatterpoints
###

# data
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes
failures = "line+node"

# xs = braessness_lines_Pmech
# ys = braessness_lines_BH
# zs = braessness_lines
# failures = "lines"

# xs = braessness_nodes_Pmech
# ys = braessness_nodes_BH
# zs = braessness_nodes
# failures = "nodes"

fig = plot_braessness_power_vs_inertia_vs_full_2D(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_power_vs_inertia_vs_full_2D_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_power_vs_inertia_vs_full_2D_failures=$failures,logcount_colorscale,,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,dim=$dimensions,fancy=$fancy_colors,center=$plot_center_only.png"),fig)


### 
### Braessness vs Braessness: Scatter using colorscaling for counts
###

# lines and nodes
failures = "line & node"
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines .+ braessness_nodes
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcountcolorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

# lines 
failures = "line"
xs = braessness_lines_Pmech
ys = braessness_lines_BH
zs = braessness_lines
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)

# nodes
failures = "node"
xs = braessness_nodes_Pmech
ys = braessness_nodes_BH
zs = braessness_nodes
fig = plot_braessness_vs_braessness_2D_count(xs, ys, zs, failures, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "exp=$exp_nrs,braessness_vs_braessness_2D_count_failures=$failures,logcount_colorscale,N=$N_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


###
### 3D scatter plot counts
###
using GLMakie # for interactivity
GLMakie.activate!()

# lines and nodes
xs = braessness_lines_BH .+ braessness_nodes_BH
ys = braessness_lines_Pmech .+ braessness_nodes_Pmech
zs = braessness_lines .+ braessness_nodes

fig = plot_braessness_power_vs_inertia_vs_full_3D(xs, ys, zs, f_b_narrow, f_b_wide, I)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_power_vs_inertia_vs_full_3D_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(MA_DIR, "braessness", "braessness_power_vs_inertia_vs_full_3D_lines_and_nodes_3D_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)


################################################################################
########################## Old plotting scripts ################################
################################################################################

###
### first scatter plot
###
fontsize = 25
titlesize = (fontsize-5)
markersize = 8
model = "full failure"
# model = exp_params_dict[:node_failure_model]

fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. ρ_Pmech_Pflow: $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,1] = ax31 = Axis(fig; xlabel="ρ_Pmech_Pflow", ylabel="N_wide-N_narrow: node failures")
fig[1,2] = ax12 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. dist(src,dst): $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,2] = ax22 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,2] = ax32 = Axis(fig; xlabel="dist(src,dst)", ylabel="N_wide-N_narrow: node failures")

scatter!(ax11, x_ρ, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
scatter!(ax21, x_ρ, braessness_lines; color=:orange, markersize=markersize)
scatter!(ax31, x_ρ, braessness_nodes; color=:red, markersize=markersize)
scatter!(ax12, x_dist, braessness_lines .+ braessness_nodes; color= (:blue, 0.1), markersize=markersize+10)
scatter!(ax22, x_dist, braessness_lines; color= (:orange, 0.1), markersize=markersize+10)
scatter!(ax32, x_dist, braessness_nodes; color = (:red, 0.1), markersize=markersize+10)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "test_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "test_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

###
### density plots Braessness vs ρ_Pmech_Pflow
###
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="Density", title="Braessness line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; xlabel="ρ_Pmech_Pflow", ylabel="Braessness")
fig[2,2] = ax22 = Axis(fig; xlabel="Density")

# link their axes
linkxaxes!(ax21, ax11)
linkyaxes!(ax21, ax22)

labels = ["line & node failures", "line failures", "node failures"]
ys     = [braessness_lines .+ braessness_nodes, braessness_lines, braessness_nodes]
cols   = (:blue, :orange, :red)

for (lbl, y, col) in zip(labels, ys, cols)
    density!(ax11, x_ρ; color=col)
    scatter!(ax21, x_ρ, y; label=lbl, color=col, markersize=markersize)
    density!(ax22, y; direction=:y, color=col)
end

hidexdecorations!(ax11)
hideydecorations!(ax22)
axislegend(ax21)
fig

###
### violin plot
###
scale = :count
gap = -0.2
ymin = -3
ymax = 3
# model = exp_params_dict[:node_failure_model]

fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; ylabel="N_wide-N_narrow: node & line failures", title="Braessness vs. dist(src,dst): $model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I,scale=$scale", titlealign = :left, titlesize = titlesize)
fig[2,1] = ax21 = Axis(fig; ylabel="N_wide-N_narrow: line failures")
fig[3,1] = ax31 = Axis(fig; xlabel="dist(src,dst)", ylabel="N_wide-N_narrow: node failures")

scale = :count
violin!(ax11, x_dist, braessness_lines .+ braessness_nodes; scale=scale, gap=gap, color=:blue)
violin!(ax21, x_dist, braessness_lines; scale=scale, gap=gap, color=:orange)
violin!(ax31, x_dist, braessness_nodes; scale=scale, gap=gap, color=:red)
# ylims!(ax11, ymin, ymax); ylims!(ax21, ymin, ymax); ylims!(ax31, ymin, ymax)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_dist_braessness_$model,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig


### 
### Braessness vs Braessness: Scatter without colorscaling (old)
###
fig = Figure(size=(2100,1500), fontsize= fontsize)

# lines and nodes
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness node & line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_lines_BH .+ braessness_nodes_BH, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_lines_Pmech .+ braessness_nodes_Pmech, braessness_lines .+ braessness_nodes; color=:blue, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_lines_Pmech .+ braessness_nodes_Pmech, braessness_lines_BH .+ braessness_nodes_BH; color=:blue, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# lines 
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness line failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_lines_BH, braessness_lines; color=:orange, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_lines_Pmech, braessness_lines; color=:orange, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_lines_Pmech, braessness_lines_BH; color=:orange, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

# nodes
fig = Figure(size=(2100,1500), fontsize= fontsize)
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="full failure", title="Braessness node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left, titlesize = titlesize)
scatter!(ax11, braessness_nodes_BH, braessness_nodes; color=:red, markersize=markersize)
fig[2,1] = ax21 = Axis(fig; xlabel="power failure",ylabel="full failure",)
scatter!(ax21, braessness_nodes_Pmech, braessness_nodes; color=:red, markersize=markersize)
fig[1,2] = ax12 = Axis(fig; xlabel="power failure",ylabel="inertia failure",)
scatter!(ax12, braessness_nodes_Pmech, braessness_nodes_BH; color=:red, markersize=markersize)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_nodes,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig

###
### (Pmech=-1, BH=-1) has high count
### 
xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines.+ braessness_nodes

freqs = Dict{Tuple{Int,Int}, Int}()
for (x,y) in zip(xs, ys)
    freqs[(x,y)] = get(freqs, (x,y), 0) + 1
end

freqs[(0,0)] # 5747
freqs[(1,0)] # 21
freqs[(0,1)] # 8
freqs[(1,1)] # 21
freqs[(1,-1)] # 22
freqs[(-1,-1)] # 156

freqs = Dict{NTuple{3,Int}, Int}()
for (x,y,z) in zip(xs, ys, zs)
    freqs[(x,y,z)] = get(freqs, (x,y,z), 0) + 1
end

freqs[(0,0,0)] # 5731
freqs[(1,0,0)] # 0
freqs[(0,1,0)] # 1
freqs[(1,1,0)] # 0
freqs[(1,-1,0)] # 0
freqs[(-1,-1,0)] # 5

xs_outlier = Int64[]; ys_outlier = Int64[]; zs_outlier = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if freqs[(x,y,z)] < count
        push!(xs_outlier, x); push!(ys_outlier, y); push!(zs_outlier, z)
    end
end

###
### Looking at the center (small Braessness) 
###
markersize = 20

xs = braessness_lines_Pmech .+ braessness_nodes_Pmech
ys = braessness_lines_BH .+ braessness_nodes_BH
zs = braessness_lines.+ braessness_nodes


xs_wo_center = Int64[]; ys_wo_center = Int64[]; zs_wo_center = Int64[]
for (x,y,z) in zip(xs, ys, zs)
    if x != 0 || y != 0
        push!(xs_wo_center, x); push!(ys_wo_center, y); push!(zs_wo_center, z)
    end
end

xs = xs_wo_center
ys = ys_wo_center
zs = zs_wo_center

## count
fig = Figure()
fig[1,1] = ax11 = Axis(fig; xlabel="inertia failure",ylabel="power failure", title="Braessness line & node failures: f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I", titlealign = :left)
logc = get_log_count(xs, ys)
sc11 = scatter!(ax11, xs, ys; markersize=markersize, color=logc, colormap=:viridis, colorrange=(minimum(logc), maximum(logc)))
Colorbar(fig[1,2], sc11; label = "log₁₀(counts)", width = 30)
lim = 10.5
xlims!(ax11, -lim, lim); ylims!(ax11, -lim, lim)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes_center_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.pdf"),fig)
# CairoMakie.save(joinpath(exp_data_dir, "braessness_lines", "braessness_vs_braessness_lines_and_nodes_center_logcount_colorscale,f_b_n=$f_b_narrow,f_b_w=$f_b_wide,I=$I.png"),fig)
fig