"""
Plots conceptual model. For more context see `slack_gen_analytic_vs_numeric.jl`
"""



using DynamicCascades
using GraphMakie
using Colors
using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()
# display(inspect_solution(sol))

include(abspath(@__DIR__, "paper_plots_helpers_and_parameters.jl"))


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

############################# helper functions #################################
function R_squared(y, f)
    e = y .- f
    y_mean = mean(y)

    SS_res = sum(e.^2)
    SS_tot = sum((y .- y_mean).^2)

    R_2 = 1 - SS_res / SS_tot

    return R_2
end


################################################################################
############################### plot nadirs ####################################
################################################################################
patchsize = (32, 30)
colors = ["#5C3D99FF", "#56B4E9FF", "#D55E00FF", "#D55E00FF"]
colors = ["#1A1A1AFF", "#595959FF", "#999999FF", "#D9D9D9FF"]
colors = ["#000000FF", "#4D4D4DFF", "#8C8C8CFF", "#CFCFCFFF"]
linestyles = [:solid, :dash, :dot, :dashdot]


# phase nadir
fig = Figure(size=(800,600), fontsize = fontsize, titlefont = "TeX Gyre Heros Makie", titlefontweight = :regular)
ax = Axis(fig[1, 1]; title="Conceptual model", titlefont = "TeX Gyre Heros Makie", xlabel = L"Inertia $I$ [$s^2$]", ylabel= L"$|S^{\mathrm{max}}|$ [p.u.]", xticks=0:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
# for D in D_values
for (index, D) in enumerate(D_values)    
    # y = phase_nadir.(x, D, K, P0, P1)
    y = 2*abs.(sin.(phase_nadir.(x, D, K, P0, P1) / 2))
    lines!(ax, x, y; label="D=$D s", linewidth=4, color = colors[index], linestyle=linestyles[index])
end
axislegend(ax, position=:rb, nbanks = 2, patchsize = patchsize, labelsize=fontsize-5)
Label(fig[1, 1, TopLeft()], fontsize=fontsize+labellettersize, "e", font = :bold, padding = (-95, 0, 5, 0))
ax.xticks = [0,2,4,6,8,10,12]

save(joinpath(MA_DIR, "phase_nadir_different_dampings.pdf"), fig)
# save(joinpath(MA_DIR, "phase_nadir_different_dampings.png"), fig)
fig

# frequency nadir
fig = Figure(size=(800,600), fontsize = fontsize)
ax = Axis(fig[1, 1]; title="Conceptual model", titlefont = "TeX Gyre Heros Makie", xlabel= L"Inertia $I$ [$s^2$]", ylabel= L"$|f^{\mathrm{max}}|$ [Hz]", xticks=0:1:12)
x = range(0.0, 12, length=1200)
D_values = [0.1, 0.5, 1.0, 2.0]
for (index, D) in enumerate(D_values)    
    # y = frequency_nadir.(x, D, K, P0, P1)
    y = frequency_nadir.(x, D, K, P0, P1)/(2π)
    lines!(ax, x, y; label="D=$D s", linewidth=4, color = colors[index], linestyle=linestyles[index])
end
axislegend(patchsize = patchsize)
Label(fig[1, 1, TopLeft()], fontsize=fontsize+labellettersize, "f", font = :bold, padding = (-85, 0, 5, 0))
ax.xticks = [0,2,4,6,8,10,12]

save(joinpath(MA_DIR, "frequency_nadir_different_dampings.pdf"), fig)
# save(joinpath(MA_DIR, "frequency_nadir_different_dampings.png"), fig)
fig