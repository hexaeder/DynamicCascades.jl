using NetworkDynamics
using NetworkDynamicsInspector
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using OrdinaryDiffEq
# using OrdinaryDiffEqTsit5

using CairoMakie # for normal plots
CairoMakie.activate!()
using WGLMakie # for inspector
# using Bonito # for using plot pane and memorizing plots
# Bonito.set_cleanup_time!(720)

@mtkmodel DynLoad begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
    end
    @parameters begin
        τ = 1.0, [description = "dyn load time constant"]
        Pload = 1.0, [description = "Load Power"]
    end
    @equations begin
        Dt(θ) ~ 1/τ * (-Pload + Pel)
        Pinj ~ -Pel
    end
end

@mtkmodel SwingDynLoad begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        ω(t) = 0.0, [description = "Rotor frequency"]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 1, [description = "Damping"]
        Pmech = 1.0, [description = "Mechanical Power"]
        τ = 1.0, [description = "dyn load time constant"]
        Pload = 1.0, [description = "Load Power"]
        # parameters for callback
        functional = 1, [description = "If 1, the node is functional (gen + load), if 0, it is only a load"]
        ωmax = Inf, [description = "Maximum rotor frequency, used in callback"]
    end
    @equations begin
        Dt(θ) ~ ifelse(functional > 0 ,
            ω,                  # gen mode
            1/τ * (-Pload + Pel) # load mode
        )
        Dt(ω) ~ ifelse(functional > 0,
            1/M * (Pmech - Pload - D*ω + Pel), # gen mode
            0.0                                # load mode
        )
        Pinj ~ -Pel
    end
end

@mtkmodel StaticPowerLine begin
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        P(t), [description = "active power flow towards node at dst end", output=true]
        Δθ(t), [description = "voltage angle difference"]
        # srcV(t) = 1.0, [description = "voltage magnitude at src end"]
        # dstV(t) = 1.0, [description = "voltage magnitude at dst end"]
        S(t), [description = "apparent power flow towards node at dst end"]
    end
    @parameters begin
        K = 3, [description = "coupling constant"]
        # parameters for callback
        active = 1, [description = "If 1, the line is active, if 0, it is tripped"]
        rating = Inf, [description = "active power line rating"]
    end
    @equations begin
        Δθ ~ srcθ - dstθ
        P ~ active*K*sin(Δθ)
        # S ~ active * max(srcV, dstV) * K * √(srcV^2 + dstV^2 - 2*srcV*dstV*cos(Δθ))
        S ~ active*K*√(2 - 2*cos(Δθ))
    end
end

function SwingDynLoadModel(; vidx=nothing, kwargs...)
    model = SwingDynLoad(name = :swing_dyn_load; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)

    # define callback, but only if :ωmax != Inf
    get_default(vm, :ωmax) == Inf && return vm

    cond = ComponentCondition([:ω], [:ωmax]) do u, p, t
        abs(u[:ω]) - p[:ωmax]
    end
    affect = ComponentAffect([:ω], [:functional]) do u, p, ctx
        println("Vertex $(ctx.vidx) tripped at t=$(ctx.integrator.t)")
        u[:ω] = 0.0
        p[:functional] = 0 # TODO hier muss man vermutlich den index hinzufügen für vector CB
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(vm, cb)
    vm
end

function DynLoadModel(; vidx=nothing, kwargs...)
    model = DynLoad(name = :load; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)
    vm
end

function Line(; src=nothing, dst=nothing, kwargs...)
    model = StaticPowerLine(name = :line; kwargs...)
    em = EdgeModel(model, [:srcθ], [:dstθ], AntiSymmetric([:P]))
    !isnothing(src) && set_graphelement!(em, src => dst)

    # define callback, but only if :rating != Inf
    get_default(em, :rating) == Inf && return em

    cond = ComponentCondition([:S], [:rating]) do u, p, t
        u[:S] - p[:rating]
    end
    affect = ComponentAffect([], [:active]) do u, p, ctx
        println("Line $(ctx.eidx) tripped at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(em, cb)
    em
end