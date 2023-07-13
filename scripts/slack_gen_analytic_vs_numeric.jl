using Revise
using DynamicCascades
using NetworkDynamics
using Graphs
using MetaGraphs
using Unitful
using Statistics
using GraphMakie
using Colors
using DynamicCascades: PLOT_DIR

using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()
display(inspect_solution(sol))

# parameters
M = 1.0
D = 0.5
K = 1.0
P0 = 0.0
P1 = P_perturb = 1.0
failtime = 1.0
Z = 4*M*K - D^2

main()

function main()
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
    fig = Figure(resolution=(1800,1000), fontsize=35)

    # phase | flow #################################################################
    # fig[1,1] = ax = Axis(fig; xlabel="time t in s", ylabel="phase angle θ | apparent power flow in p.u.", title="Phase | Power flow")
    fig[1,1] = ax = Axis(fig; xlabel="time t [s]", ylabel="phase angle θ", title="Phase")
    # phase of gen
    t = sol.sol.t
    y = [sol.sol.u[i][1] for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="numeric", linewidth=3)
    # # phase of slack
    # y = [sol.sol.u[i][3] for i in 1:length(sol.sol.t)]
    # lines!(ax, t, y; label="slack", linewidth=3)

    # phase analytic solution
    y = [phase(t[i], failtime) for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="analytic", linewidth=3)

    # # power flow
    # t = sol.load_S.t
    # y = [sol.load_S.saveval[i][1] for i in 1:length(sol.load_S.t)]
    # lines!(ax, t, y; label="flow", linewidth=3)
    axislegend()

    # frequency ####################################################################
    fig[1,2] = ax = Axis(fig; xlabel="time t [s]", ylabel="ang. frequency ω [Hz]", title="Frequency")
    # frequency of gen
    t = sol.sol.t
    y = [sol.sol.u[i][2] for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="numeric", linewidth=3)

    # frequency analytic solution
    y = [frequency(t[i], failtime) for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="analytic", linewidth=3)

    # # frequency of slack
    # t = sol.frequencies_load_nodes.t
    # y = [sol.frequencies_load_nodes.saveval[i][1] for i in 1:length(sol.frequencies_load_nodes.t)]
    # lines!(ax, t, y; label="slack", linewidth=3)
    axislegend()

    # RoCoF ########################################################################
    fig[2,2] = ax = Axis(fig; xlabel="time t [s]", ylabel="RoCoF [Hz/s]", title="Rate of change of angular frequency")
    t = sol.sol.t
    y = [sol.sol(t[i], Val{1})[2] for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="numeric", linewidth=3)

    # RocoF analytic solution
    y = [RocoF(t[i], failtime) for i in 1:length(sol.sol.t)]
    lines!(ax, t, y; label="analytic", linewidth=3)
    axislegend()

    parameters = "Parameters: inertia M=$M [s^2], damping D=$D [s], coupling K=$K, power perturbation P_perturb=$P_perturb [p.u.]"
    supertitle = Label(fig[0, :], parameters)
    # sideinfo = Label(fig[3, 1:2], parameters)
    fig

    save(joinpath(PLOT_DIR, "nadir_sim_M=$M,D=$D,K=$K,P_perturb=$P_perturb.pdf"), fig)
    save(joinpath(PLOT_DIR, "nadir_sim_M=$M,D=$D,K=$K,P_perturb=$P_perturb.png"), fig)
end




############################# analytic functions ################################

function calc_constants(M, D, K, P0, P1)
    w = sqrt(4*M*K-D^2) / (2*M)
    g = D / (2*M)
    β = (P0 - P1) / (K*sin(atan(w / g)))
    α = atan(w / g)
    return w, g, β, α
end

################################# phase angle ##################################
function phase(t, failtime; M=M, D=D, K=K, P0=P0, P1=P1)
    if t < failtime
        return 0
    end
    w, g, β, α = calc_constants(M, D, K, P0, P1)
    θ = β * exp(-g * (t - failtime)) * sin(w * (t - failtime) + α) + P1 / K
    return θ
end

################################# frequency ####################################
function frequency(t, failtime; M=M, D=D, K=K, P0=P0, P1=P1)
    if t < failtime
        return 0
    end
    w, g, β, α = calc_constants(M, D, K, P0, P1)
    ω = β * exp(-g * (t - failtime)) * (cos(w * (t - failtime) + α) * w - sin(w * (t - failtime) + α) * g)
    return ω
end

################################# RocoF ####################################
function RocoF(t, failtime; M=M, D=D, K=K, P0=P0, P1=P1)
    if t < failtime
        return 0
    end
    w, g, β, α = calc_constants(M, D, K, P0, P1)
    RocoF = β * exp(-g * (t - failtime)) * (sin(w * (t - failtime) + α) * (g^2 -w^2) - 2 * cos(w * (t - failtime) + α) * w * g)
    return RocoF
end
