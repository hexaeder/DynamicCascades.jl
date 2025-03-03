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
############################### plot nadirs ####################################
################################################################################
# phase nadir
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel= L"τ1", ylabel= L"$\Delta \theta$ [rad]", xticks=1:1:12)
τ1 = range(0.0, 3, length=1200)
τ2_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for τ2 in τ2_values
    y = phase_nadir.((τ1*τ2), τ1, K, P0, P1)
    lines!(ax, τ1, y; label="τ2=$τ2", linewidth=4)
end
axislegend(ax, position=:rt, nbanks = 2)
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "phase_nadir_tau_variation.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "phase_nadir_tau_variation.png"), fig)

# phase nadir
fig = Figure(fontsize = 50)
ax = Axis(fig[1, 1]; xlabel= L"τ1", ylabel= L"$\Delta \theta$ [rad]", xticks=1:1:12)
τ1 = range(0.0, 3, length=1200)
τ2_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for τ2 in τ2_values
    y = phase_nadir.((τ1*τ2), τ1, K, P0, P1)
    lines!(ax, τ1, y; label="τ2=$τ2", linewidth=4)
end
y = phase_nadir.(τ1.^2, τ1, K, P0, P1)
lines!(ax, τ1, y; label="τ2=τ1", color=:red, linewidth=4)
y = phase_nadir.(τ1.^2*c, τ1, K, P0, P1)
lines!(ax, τ1, y; label="τ2=c*τ1, c=$c", color=:red, linestyle=:dot, linewidth=4)

axislegend(ax, position=:rt, nbanks = 2)
fig
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "phase_nadir_combined_tau.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "phase_nadir_combined_tau.png"), fig)

fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel=L"$\Delta \theta$ [rad]", xticks=1:1:12)
x = range(0.0, 3, length=1200)
M_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for M in M_values
    y = phase_nadir.(M, x, K, P0, P1)
    lines!(ax, x, y; label="Inertia I=$M", linewidth=4)
end
axislegend()
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "phase_nadir_D_variation.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "phase_nadir_D_variation.png"), fig)

# frequency nadir
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel= L"τ1", ylabel= L"$\Delta \omega$ [rad/s]")
τ1 = range(0.01, 3, length=1200)
τ2_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for τ2 in τ2_values
    y = frequency_nadir.((τ1*τ2), τ1, K, P0, P1)
    lines!(ax, τ1, y; label="τ2=$τ2", linewidth=4)
end
axislegend()
ylims!(ax, 0., 3.)
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_tau_variation.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_tau_variation.png"), fig)

fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel=L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
M_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for M in M_values
    y = frequency_nadir.(M, x, K, P0, P1)
    lines!(ax, x, y; label="Inertia I=$M", linewidth=4)
end
axislegend()
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_D_variation.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_D_variation.png"), fig)

# I/D^2=const.=1 in abh. von I.
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel="inertia I", ylabel=L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
y = frequency_nadir.(x, sqrt.(x), K, P0, P1)
lines!(ax, x, y; label="D = sqrt(I)", linewidth=4)
axislegend()
# ylims!(ax, 0., 3.)
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_I_over_D_squared_constant_I.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_I_over_D_squared_constant_I.png"), fig)

# I/D^2=const.=1 in abh. von D.
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel=L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
y = frequency_nadir.(x.^2, x, K, P0, P1)
lines!(ax, x, y; label="I = D^2", linewidth=4)
axislegend()
xlims!(ax, -0.1, 3.)
ylims!(ax, 0., 3.)
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_I_over_D_squared_constant_D.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_I_over_D_squared_constant_D.png"), fig)


# combined
fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel= L"τ1", ylabel= L"$\Delta \omega$ [rad/s]")
τ1 = range(0.00, 3, length=1200)
τ2_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for τ2 in τ2_values
    y = frequency_nadir.((τ1*τ2), τ1, K, P0, P1)
    lines!(ax, τ1, y; label="τ2=$τ2", linewidth=4)
end
y = frequency_nadir.(τ1.^2, τ1, K, P0, P1)
lines!(ax, τ1, y; label="τ2=τ1", color =:red, linewidth=4)
y = frequency_nadir.(τ1.^2*c, τ1, K, P0, P1)
lines!(ax, τ1, y; label="τ2=c*τ1, c=$c", color=:red, linestyle=:dot, linewidth=4)
axislegend()
xlims!(ax, -0.1, 3.)
ylims!(ax, 0., 3.)
fig
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "frequency_nadir_combined_tau.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "frequency_nadir_combined_tau.png"), fig)

fig = Figure(fontsize = 40)
ax = Axis(fig[1, 1]; xlabel="damping D", ylabel=L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
M_values = [0.1, 0.5, 1.0, 2.0, 3.0]
for M in M_values
    y = frequency_nadir.(M, x, K, P0, P1)
    lines!(ax, x, y; label="Inertia I=$M", linewidth=4)
end
y = frequency_nadir.(x.^2, x, K, P0, P1)
lines!(ax, x, y; label="I = D^2", color=:red, linewidth=4)
axislegend()
xlims!(ax, -0.1, 3.)
ylims!(ax, 0., 3.)
fig
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_combined_D.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/Private_MA/paper/relation_I_D", "frequency_nadir_combined_D.png"), fig)

c = 0.7
fig = Figure(fontsize = 25)
ax = Axis(fig[1, 1]; xlabel= L"Inertia I [$s^2$]", ylabel= L"$\Delta \theta$ [rad]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
for D in D_values
    y = phase_nadir.(x, D, K, P0, P1)
    lines!(ax, x, y; label="D=$D", linewidth=4)
end
y = phase_nadir.(x, sqrt.(x), K, P0, P1)
lines!(ax, x, y; label="D = sqrt(I)", color=:red, linewidth=4)
y = phase_nadir.(x, sqrt.(x/c), K, P0, P1)
lines!(ax, x, y; label="D = sqrt(I/c), c=$c", color=:red, linestyle=:dot, linewidth=4)
y = phase_nadir.(x, x, K, P0, P1)
lines!(ax, x, y; label="D = I", color=:black, linewidth=4)
y = phase_nadir.(x, x./c, K, P0, P1)
lines!(ax, x, y; label="D = I/c, c=$c", color=:black, linestyle=:dot, linewidth=4)
axislegend(ax, position=:rc, nbanks = 3)
fig
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "phase_nadir_combined_both_I.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "phase_nadir_combined_both_I.png"), fig)

fig = Figure(fontsize = 35)
ax = Axis(fig[1, 1]; xlabel= L"Inertia I [$s^2$]", ylabel= L"$\Delta \omega$ [rad/s]", xticks=1:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
for D in D_values
    y = frequency_nadir.(x, D, K, P0, P1)
    lines!(ax, x, y; label="D=$D", linewidth=4)
end
y = frequency_nadir.(x, sqrt.(x), K, P0, P1)
lines!(ax, x, y; label="D = sqrt(I)", color=:red, linewidth=4)
y = frequency_nadir.(x, sqrt.(x/c), K, P0, P1)
lines!(ax, x, y; label="D = sqrt(I/c), c=$c", color=:red, linestyle=:dot, linewidth=4)
y = frequency_nadir.(x, x, K, P0, P1)
lines!(ax, x, y; label="D = I", color=:black, linewidth=4)
y = frequency_nadir.(x, x./c, K, P0, P1)
lines!(ax, x, y; label="D = I/c, c=$c", color=:black, linestyle=:dot, linewidth=4)
axislegend()
ylims!(ax, 0., 2.)
fig
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "frequency_nadir_combined_both_I.pdf"), fig)
save(joinpath("/home/brandner/nb_data/repos/NLCp/paper_figs/relation_I_D/", "frequency_nadir_combined_both_I.png"), fig)

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
