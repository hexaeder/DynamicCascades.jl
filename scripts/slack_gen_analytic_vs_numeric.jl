using Revise
using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()
# display(inspect_solution(sol))


################################################################################
################################# parameters ###################################
################################################################################
M = 1.0
D = 0.5
K = 1.0
P0 = 0.0
P1 = P_perturb = 0.93
failtime = 1.0
Z = 4*M*K - D^2

################################################################################
############################ numeric vs. analytic model ########################
################################################################################
# NOTE NOTE NOTE NOTE in ND_model.jl uncomment `load_model = :slack`
network = import_system(:nadir_sim; M=float(M)u"s^2", γ=float(D)u"s", K=K, tconst=0.1u"s")
sol = simulate(network;
    initial_fail = [1],
    failtime = failtime,
    init_pert = :power_perturbation,
    P_perturb = P_perturb,
    trip_lines = :none,
    trip_nodes = :none,
    trip_load_nodes = :none,
    tspan = (0.0, 50.0),
    solverargs = (; dtmax = 0.01));

# plot solution
fig = Figure(resolution=(1800,500), fontsize=45)
# fig = Figure(fontsize=35)

# phase | flow #################################################################
# fig[1,1] = ax = Axis(fig; xlabel="time t in s", ylabel="phase angle θ | apparent power flow in p.u.", title="Phase | Power flow")
fig[1,1] = ax = Axis(fig; xlabel= L"Time $t$ [s]", ylabel=L"$\theta$ [rad]")
# phase of gen
t = sol.sol.t
y_numeric = [sol.sol.u[i][1] for i in 1:length(sol.sol.t)]
# lines!(ax, t, y_numeric; label="numeric phase angle", color=:purple, linewidth=3)
lines!(ax, t, y_numeric; label="numeric phase", linewidth=3)
# # phase of slack
# y = [sol.sol.u[i][3] for i in 1:length(sol.sol.t)]
# lines!(ax, t, y; label="slack", linewidth=3)

# phase analytic solution
t = sol.sol.t
y_analytic = [phase(t[i], failtime, M, D, K, P0, P1) for i in 1:length(sol.sol.t)]
# lines!(ax, t, y_analytic; label="analytic phase angle", color=:orange, linewidth=3)
lines!(ax, t, y_analytic; label="analytic phase", linewidth=3)

# power flow numeric
t = sol.load_S.t
y_powerflow_numeric = [sol.load_S.saveval[i][1] for i in 1:length(sol.load_S.t)]
# lines!(ax, t, y_powerflow_numeric; label="numeric power flow", color=:purple, linewidth=3, linestyle=:dash)
lines!(ax, t, y_powerflow_numeric; label="numeric flow", color=:red, linewidth=3, linestyle=:dash)

# power flow analytic
t = sol.sol.t
y_powerflow_analytic = [K * phase(t[i], failtime, M, D, K, P0, P1) for i in 1:length(sol.sol.t)]
# lines!(ax, t, y_powerflow_analytic; label="analytic power flow", color=:pink, linewidth=3, linestyle=:dash)
lines!(ax, t, y_powerflow_analytic; label="analytic flow", color=:black, linewidth=3, linestyle=:dash)
# scatter!(ax, t, y_powerflow_analytic; label="power flow analytic", marker=:cross)

# add parameters/coefficients to legend
R_2 = round(R_squared(y_analytic, y_numeric), digits = 3)
lines!(ax, [NaN], [NaN]; label="R²=$R_2, P'=$P1", color=:white, linewidth=3)
# ylims!(-1, 15)
# xlims!(0, 32)
# axislegend(ax, position = :rb, nbanks=2, labelsize=32)
axislegend(ax, position = :rt, nbanks=1, labelsize=32)
fig
# frequency ####################################################################
fig[1,2] = ax = Axis(fig; xlabel=L"Time $t$ [s]", ylabel=L"$\omega$ [rad/s]")
# frequency of gen
t = sol.sol.t
y_numeric = [sol.sol.u[i][2] for i in 1:length(sol.sol.t)]
lines!(ax, t, y_numeric; label="numeric", linewidth=3)

# frequency analytic solution
y_analytic = [frequency(t[i], failtime, M, D, K, P0, P1) for i in 1:length(sol.sol.t)]
lines!(ax, t, y_analytic; label="analytic", linewidth=3)

# # frequency of slack
# t = sol.frequencies_load_nodes.t
# y = [sol.frequencies_load_nodes.saveval[i][1] for i in 1:length(sol.frequencies_load_nodes.t)]
# lines!(ax, t, y; label="slack", linewidth=3)

# add parameters/coefficients to legend
R_2 = round(R_squared(y_analytic, y_numeric), digits = 3)
lines!(ax, [NaN], [NaN]; label="R²=$R_2, P'=$P1", color=:white, linewidth=3)
axislegend(ax, labelsize=32)
# xlims!(0, 32)
# axislegend(ax, position=(1.0,.35), labelsize=32)
fig

# # RoCoF ########################################################################
# fig[1,3] = ax = Axis(fig; xlabel="time t [s]", ylabel="RoCoF [rad/s^2]", title="Rate of change of angular frequency")
# t = sol.sol.t
# y_numeric = [sol.sol(t[i], Val{1})[2] for i in 1:length(sol.sol.t)]
# lines!(ax, t, y_numeric; label="numeric", linewidth=3)
#
# # RocoF analytic solution
# y_analytic = [RocoF(t[i], failtime, M, D, K, P0, P1) for i in 1:length(sol.sol.t)]
# lines!(ax, t, y_analytic; label="analytic", linewidth=3)
#
# # add parameters/coefficients to legend
# R_2 = round(R_squared(y_analytic, y_numeric), digits = 3)
# lines!(ax, [NaN], [NaN]; label="R^2=$R_2", color=:white, linewidth=3)
# axislegend()
#
# # parameters = "Parameters: inertia M=$M [s^2], damping D=$D [s], coupling K=$K, power perturbation P_perturb=$P_perturb [p.u.]"
# # supertitle = Label(fig[0, :], parameters)
# # sideinfo = Label(fig[3, 1:2], parameters)
# fig

save(joinpath(MA_DIR, "nadir_sim_M=$M,D=$D,K=$K,P_perturb=$P_perturb.pdf"), fig)
save(joinpath(MA_DIR, "nadir_sim_M=$M,D=$D,K=$K,P_perturb=$P_perturb.png"), fig)

################################################################################
################### Testing if nadirs are computed correctly ###################
################################################################################

############################## phase nadir plot ################################
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="t | M", ylabel="phase | phase nadir", xticks=1:1:12)
x = range(0.0, 12, length=1200)
y = phase_nadir.(x, 0.5, K, P0, P1)
lines!(ax, x, y; label="phase nadir")
y = phase.(x, 0, M, D, K, P0, P1)
lines!(ax, x, y; label="phase with M=$M")
vlines!(ax, M; color=:black, linewidth=1, label="M=$M")
hlines!(ax, phase_nadir.(M, 0.5, K, P0, P1); color=:pink, linewidth=1)
axislegend()
fig
save(joinpath(PLOT_DIR, "test_phase_nadir_M=$M,D=$D,K=$K,P_perturb=$P_perturb.png"), fig)

############################## frequency nadir #################################
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="t | M", ylabel="frequency | frequency nadir", xticks=1:1:12)
x = range(0.0, 12, length=1200)
y = frequency_nadir.(x, 0.5, K, P0, P1)
lines!(ax, x, y; label="frequency nadir")
y = frequency.(x, 0, M, D, K, P0, P1)
lines!(ax, x, y; label="frequency with M=$M")
vlines!(ax, M; color=:black, linewidth=1, label="M=$M")
hlines!(ax, frequency_nadir.(M, 0.5, K, P0, P1); color=:pink, linewidth=1)
axislegend()
fig
save(joinpath(PLOT_DIR, "test_frequency_nadir_M=$M,D=$D,K=$K,P_perturb=$P_perturb.png"), fig)

################################################################################
############################### plot nadirs ####################################
################################################################################
# phase nadir
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel= L"Inertia I [$s^2$]", ylabel= L"$\Delta \theta$ [rad]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
for D in D_values
    y = phase_nadir.(x, D, K, P0, P1)
    lines!(ax, x, y; label="D=$D", linewidth=4)
end
axislegend(ax, position=:rb, nbanks = 2)
fig
save(joinpath(MA_DIR, "phase_nadir_different_dampings.pdf"), fig)
save(joinpath(MA_DIR, "phase_nadir_different_dampings.png"), fig)

fig = Figure(fontsize = 30)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel="max. absolute phase deviation ΔΘ", xticks=1:1:12)
x = range(0.0, 12, length=1200)
M_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for M in M_values
    y = phase_nadir.(M, x, K, P0, P1)
    lines!(ax, x, y; label="M=$M", linewidth=4)
end
axislegend()
fig
save(joinpath(MA_DIR, "phase_nadir_different_inertia.pdf"), fig)

# frequency nadir
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel= L"Inertia I [$s^2$]", ylabel= L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
for D in D_values
    y = frequency_nadir.(x, D, K, P0, P1)
    lines!(ax, x, y; label="D=$D", linewidth=4)
end
axislegend()
fig

############# analytic Δω
#############

# Define the analytical expression for the maximum frequency deviation
function frequency_nadir_analytic(I, D, K, P0, P1)
    # Guard against nonpositive inertia
    if I <= 0
        return NaN
    end
    # Compute the power step
    deltaP = P1 - P0
    # Natural frequency
    ωₙ = sqrt(K / I)
    # Damping ratio
    ζ = D / (2 * sqrt(K * I))
    # If the system is overdamped (ζ >= 1) return NaN (or handle appropriately)
    if ζ >= 1
        return NaN
    end
    # Analytical expression:
    # Δω_max = (ΔP/K)*ωₙ * exp[ - (ζ*arccos(ζ)) / sqrt(1-ζ^2) ]
    return (deltaP / K) * ωₙ * exp( - (ζ * acos(ζ)) / sqrt(1 - ζ^2) )
end

# Create a figure and axis
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel = L"Inertia I [$s^2$]", ylabel = L"$\Delta \omega$ [rad/s]", xticks = 1:1:12)

# Define the range of inertia values (avoid I = 0 to prevent division by zero)
x = range(0, 12, length = 1200)

# Define a set of damping values to test
D_values = [0.1, 0.5, 1.0, 2.0]

# Compute and plot the frequency nadir for each damping value
for D in D_values
    y = frequency_nadir_analytic.(x, D, K, P0, P1)
    lines!(ax, x, y; label = "D = $D", linewidth = 4)
end

axislegend(ax)
fig

#############
#############

save(joinpath(MA_DIR, "frequency_nadir_different_dampings.pdf"), fig)
save(joinpath(MA_DIR, "frequency_nadir_different_dampings.png"), fig)

fig = Figure(fontsize = 30)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel="max. absolute frequency deviation Δf", xticks=1:1:12)
x = range(0.0, 12, length=1200)
M_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for M in M_values
    y = frequency_nadir.(M, x, K, P0, P1)
    lines!(ax, x, y; label="M=$M", linewidth=4)
end
axislegend()
fig
save(joinpath(MA_DIR, "frequency_nadir_different_inertia.pdf"), fig)

################################################################################
################################# M-D heatmaps #################################
################################################################################
# https://docs.makie.org/stable/examples/plotting_functions/heatmap/

################################# phase angle ##################################
fig = Figure(fontsize=30)
fig[1,1] = ax = Axis(fig; xlabel="inertia M [s^2]", ylabel="damping D [s]",
                    title="Phase angle deviation, K=$K, P_perturb=$P_perturb [p.u.]",
                    titlesize=25)
xs = range(start=0.0, step=0.1, stop = 6.0)
ys = range(start=0.0, step=0.1, stop = 5.0)
zs = [phase_nadir(x, y, K, P0, P1) for x in xs, y in ys]
hm = Makie.heatmap!(ax, xs, ys, zs, colormap = Reverse(:deep))
#= Alternative way using requiring kwargs: phase_nadir(M, D; K=K, P0=P0, P1=P1)
hm = Makie.heatmap!(ax, 0:0.1:6, 0:0.1:6, phase_nadir, colormap = Reverse(:deep))=#
Colorbar(fig[:, end+1], hm, label="maximal absolute phase angle deviation")
fig

save(joinpath(MA_DIR, "heatmap_phase_K=$K,P_perturb=$P_perturb.pdf"), fig)
save(joinpath(MA_DIR, "heatmap_phase_K=$K,P_perturb=$P_perturb.png"), fig)

################################# frequency ####################################
fig = Figure(fontsize=30)
fig[1,1] = ax = Axis(fig; xlabel="inertia M [s^2]", ylabel="damping D [s]",
                    title="Ang. frequency deviation, K=$K, P_perturb=$P_perturb [p.u.]",
                    titlesize=25)
xs = range(start=0.0, step=0.05, stop = 1.0)
ys = range(start=0.0, step=0.05, stop = 2.0)
zs = [frequency_nadir(x, y, K, P0, P1) for x in xs, y in ys]
hm = Makie.heatmap!(ax, xs, ys, zs, colormap = Reverse(:deep))
Colorbar(fig[:, end+1], hm, label="maximal absolute ang. freq. dev.")
fig

save(joinpath(MA_DIR, "zoomed_heatmap_frequency_K=$K,P_perturb=$P_perturb.pdf"), fig)
save(joinpath(MA_DIR, "zoomed_heatmap_frequency_K=$K,P_perturb=$P_perturb.png"), fig)

##################### border complex and real eigenvalues ######################
fig = Figure(fontsize=30)
ax = Axis(fig[1, 1]; xlabel=L"Inertia I [$s^2$]", ylabel=L"Damping D [$s$]",
                    # title="Complex (oscillating) vs. real eigenvalues (overdamped), K=$K",
                    # titlesize=20
                    xgridvisible = false,
                    ygridvisible = false,
                    )
x = range(0.0, 6, length=600)
y = border_complex_real.(x, K)
lines!(ax, x, y, linewidth=4, label = L"$D=2\sqrt{IK}$")
text!(ax, 1.7, 4., text = "Real-valued (overdamped)", align = (:center, :center), textsize=30)
text!(ax, 3.5, 2., text = "Complex-valued (oscillating)", align = (:center, :center), textsize=30)
axislegend(ax, position=(1.,0.75))
fig
save(joinpath(MA_DIR, "complex_vs_real_eigenvalues_K=$K.pdf"), fig)
save(joinpath(MA_DIR, "complex_vs_real_eigenvalues_K=$K.png"), fig)

function border_complex_real(D, K)
    return 2*sqrt(D*K)
end

############################### nadir functions ################################
function phase_nadir(M, D, K, P0, P1)
    if (4*M*K - D^2) < 0
        return NaN
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    t = t_θ_nadir(M, D, K, P0, P1)
    phase(t, 0, M, D, K, P0, P1)
end

function t_θ_nadir(M, D, K, P0, P1)
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    t = π / w
    second_deriv = RocoF(t, 0, M, D, K, P0, P1)
    if second_deriv > 0.0
        @warn("Second derivative is $second_deriv > 0.0 so this is a minimum.")
    end
    return t
end

function phase_nadir_analytic(M, D, K, P0, P1)
    if (4*M*K - D^2) < 0
        return NaN
    end
    return abs(P1 - P0) / K * (1 + exp((-D*π) / sqrt(4*M*K - D^2)))
end

function frequency_nadir(M, D, K, P0, P1)
    if (4*M*K - D^2) < 0
        return NaN
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    t = t_ω_nadir(M, D, K, P0, P1)
    frequency(t, 0, M, D, K, P0, P1)
end

function t_ω_nadir(M, D, K, P0, P1)
    α, β, g, w = calc_constants(M, D, K, P0, P1)

    if (g^2 - w^2) == 0.0
        @warn("No maximum found, for M=$M, D=$D, K=$K it is (g^2 - w^2) == 0.0.")
        return nothing
    end

    t = 0.0
    if (g^2 - w^2) < 0.0
        t = 1/w * (atan(2*w*g / (g^2 - w^2)) - α) + π/w
    else (g^2 - w^2) > 0.0
        t = 1/w * (atan(2*w*g / (g^2 - w^2)) - α)
    end

    second_deriv = deriv_RocoF(t, 0, M, D, K, P0, P1)
    if second_deriv > 0.0
        @warn("Second derivative is $second_deriv > 0.0 so this is a minimum.")
    end
    return t
end


function calc_constants(M, D, K, P0, P1)
    w = sqrt(4*M*K-D^2) / (2*M)
    g = D / (2*M)
    α = atan(w / g)
    β = (P0 - P1) / (K*sin(α))
    return α, β, g, w
end

# function t_ω_nadir(M, D, K, P0, P1)
#     α, β, g, w = calc_constants(M, D, K, P0, P1)
#
#     if (g^2 - w^2) < 0.0
#         t = 1/w * (atan(2*w*g / (g^2 - w^2)) - α) + π/w
#         return t
#     elseif (g^2 - w^2) > 0.0
#         t = 1/w * (atan(2*w*g / (g^2 - w^2)) - α)
#         return t
#     elseif (g^2 - w^2) == 0.0
#         @warn("No maximum found, for M=$M, D=$D, K=$K it is (g^2 - w^2) == 0.0.")
#         return nothing
#     end
# end

############################# analytic functions ###############################


function phase(t, failtime, M, D, K, P0, P1)
    if t < failtime
        return 0
    end
    if (4*M*K - D^2) < 0
        return NaN
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    θ = β * exp(-g * (t - failtime)) * sin(w * (t - failtime) + α) + P1 / K
    return θ
end

function frequency(t, failtime, M, D, K, P0, P1)
    if isnothing(t)
        return NaN
    end
    if t < failtime
        return 0
    end
    if (4*M*K - D^2) < 0
        return NaN
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    ω = β * exp(-g * (t - failtime)) * (cos(w * (t - failtime) + α) * w - sin(w * (t - failtime) + α) * g)
    return ω
end

function RocoF(t, failtime, M, D, K, P0, P1)
    if isnothing(t)
        return NaN
    end
    if t < failtime
        return 0
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    RocoF = β * exp(-g * (t - failtime)) * (sin(w * (t - failtime) + α) * (g^2 - w^2) - 2 * cos(w * (t - failtime) + α) * w * g)
    return RocoF
end

""" This is the time-derivative of the RoCoF"""
function deriv_RocoF(t, failtime, M, D, K, P0, P1)
    if isnothing(t)
        return NaN
    end
    if t < failtime
        return 0
    end
    α, β, g, w = calc_constants(M, D, K, P0, P1)
    deriv_RocoF = β * exp(-g * (t - failtime)) * (sin(w * (t - failtime) + α) * (3*w^2*g - g^3) + cos(w * (t - failtime) + α) * (3*w*g^2 - w^3))
    return deriv_RocoF
end

# function frequency_nadir_wrong(M, D, K, P0, P1)
#     if (4*M*K - D^2) < 0
#         return NaN
#     end
#     α, β, g, w = calc_constants(M, D, K, P0, P1)
#     return abs(P1 - P0) / K * exp(-(g*π/2) / w) * (w + sin(π/2 + α) * g / sin(α))
# end

############################# helper functions #################################
function R_squared(y, f)
    e = y .- f
    y_mean = mean(y)

    SS_res = sum(e.^2)
    SS_tot = sum((y .- y_mean).^2)

    R_2 = 1 - SS_res / SS_tot

    return R_2
end
