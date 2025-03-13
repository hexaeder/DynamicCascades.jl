using NetworkDynamics
using NetworkDynamicsInspector
using ModelingToolkit
using ModelingToolkit: D_nounits as Dt, t_nounits as t
using OrdinaryDiffEqTsit5
using WGLMakie

@mtkmodel SlackNode begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
    end
    @parameters begin
        θ_set = 0.0, [description = "voltage angle setpoint"]
    end
    @equations begin
        θ ~ θ_set
        Pinj ~ -Pel # TODO Only done for convenience (practical purposes) to have Pinj with loads?
    end
end

@mtkmodel DynLoad begin
    @variables begin
        θ(t) = 0.0, [description = "Voltage angle", output=true]
        Pel(t), [description = "Electical power flowing from network into node", input=true]
        Pinj(t), [description = "Electical power injected into network"]
    end
    @parameters begin
        τ = 1.0, [description = "dyn load time constant"]
        Pload = -1.0, [description = "Load Power"]
    end
    @equations begin
        Dt(θ) ~ 1/τ * (Pload + Pel)
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
        D = 0.1, [description = "Damping"]
        Pmech = 1.0, [description = "Mechanical Power"]
        τ = 1.0, [description = "dyn load time constant"]
        Pload = -1.0, [description = "Load Power"]
        # parameters for callback
        functional = 1, [description = "If 1, the node is functional (gen + load),if 0, it is only a load"]
        ωmax = Inf, [description = "Maximum rotor frequency, used in callback"]
    end
    @equations begin
        Dt(θ) ~ ifelse(functional > 0 ,
            ω,                  # gen mode
            1/τ * (Pload + Pel) # load mode
        )
        Dt(ω) ~ ifelse(functional > 0,
            1/M * (Pmech + Pload - D*ω + Pel), # gen mode
            0.0                                # load mode
        )
        Pinj ~ -Pel
    end
end

@mtkmodel StaticPowerLine begin
    @variables begin
        srcθ(t), [description = "voltage angle at src end", input=true]
        dstθ(t), [description = "voltage angle at dst end", input=true]
        P(t), [description = "flow towards node at dst end", output=true]
        Δθ(t), [description = "Voltage angle difference"]
    end
    @parameters begin
        K = -1.63, [description = "Coupling constant"]
        # parameters for callback
        active = 1, [description = "If 1, the line is active, if 0, it is tripped"]
        limit = Inf, [description = "Active power line limit"]
    end
    @equations begin
        Δθ ~ srcθ - dstθ
        P ~ -active*K*sin(Δθ)
    end
end

function SlackModel(; vidx=nothing, kwargs...)
    model = SlackNode(name = :slack; kwargs...)
    vm = VertexModel(model, [:Pel], [:θ])
    !isnothing(vidx) && set_graphelement!(vm, vidx)
    vm
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

    # define callback, but only if :limit != Inf
    get_default(em, :limit) == Inf && return em

    cond = ComponentCondition([:P], [:limit]) do u, p, t
        abs(u[:P]) - p[:limit]
    end
    affect = ComponentAffect([], [:active]) do u, p, ctx
        println("Line $(ctx.eidx) tripped at t=$(ctx.integrator.t)")
        p[:active] = 0
    end
    cb = ContinousComponentCallback(cond, affect)
    set_callback!(em, cb)
    em
end

#### Test functionality of node failure in 2 bus system
####
v1 = SlackModel(vidx=1)
v2 = SwingDynLoadModel(vidx=2, ωmax=0.1, Pmech=1, Pload=-0.5)
l = Line(src=1, dst=2)

# add a perturbation to the system
perturb = PresetTimeComponentCallback(1.0,
    ComponentAffect([], [:θ_set]) do u, p, ctx
        println("Jump slack voltage angle at t=$(ctx.integrator.t)")
        p[:θ_set] = 0.1
    end
)
set_callback!(v1, perturb)

nw = Network([v1,v2], l)
s0 = find_fixpoint(nw)
prob = ODEProblem(nw, uflat(s0), (0, 5), pflat(s0), callback=get_callbacks(nw));
sol = solve(prob, Tsit5())

# call for interactive inspection
# inspect(sol)

let
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1]; xlabel="time", ylabel="voltage angel θ")
    lines!(ax, sol, idxs=vidxs(1:2, :θ))
    axislegend(ax)
    ax = Axis(fig[2, 1]; xlabel="time", ylabel="rotor frequency ω")
    lines!(ax, sol, idxs=vidxs(2, :ω))
    axislegend(ax)
    ax = Axis(fig[3, 1]; xlabel="time", ylabel="injected power P")
    lines!(ax, sol, idxs=vidxs(1:2, :Pinj))
    axislegend(ax)
    fig
end

####
#### Test functionality of line failure in 2 bus system
####
v1 = SwingDynLoadModel(vidx=1, Pmech=1, Pload=0)
v2 = DynLoadModel(vidx=2, Pload=-1)
l = Line(src=1, dst=2, limit=1.1)
nw = Network([v1,v2], l)

# at time t=1 increase the load to 1.2
perturb = PresetTimeComponentCallback(1.0,
    ComponentAffect([], [:Pload]) do u, p, ctx
        println("Increase load at t=$(ctx.integrator.t)")
        p[:Pload] = -1.2
    end
)
set_callback!(v2, perturb)

s0 = find_fixpoint(nw)
@assert s0.e[1, :P] < 1.0 # fixpoint is still valid
prob = ODEProblem(nw, uflat(s0), (0, 2), pflat(s0), callback=get_callbacks(nw));
sol = solve(prob, Tsit5())


fig = Figure()
ax = Axis(fig[1, 1]; xlabel="time2", ylabel="powerflow through line")
lines!(ax, sol, idxs=eidxs(1, :P))
fig

# call for interactive inspection
# inspect(sol)