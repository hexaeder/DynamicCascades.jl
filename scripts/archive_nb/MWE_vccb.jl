function f(du, u, p, t)
    du[1] = u[2]
    du[2] = -p
    du[3] = u[4]
    du[4] = 0.0
end

function condition(out, u, t, integrator) # Event when condition(out,u,t,integrator) == 0
    # out[1] = u[1]
    # out[2] = (u[3] - 10.0)u[3]
    out[1] = -1.0
    out[2] = integrator.t - 3.0
    out[1] = -1.0
    out[2] = t - 3.0
end

function affect!(integrator, idx)
    t = integrator.t
    println("affect triggered in idx $idx at time $t \n")
end

cb = VectorContinuousCallback(condition, affect!, 2)

u0 = [50.0, 0.0, 0.0, 2.0]
tspan = (0.0, 6.0)
p = 9.8
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, Tsit5(), callback = cb, dt = 1e-3, adaptive = false)
