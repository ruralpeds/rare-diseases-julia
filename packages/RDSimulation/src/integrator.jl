# Explicit fourth-order Runge-Kutta integrator.
#
# Hand-rolled to keep RDSimulation dependency-light. For stiff systems or
# adaptive tolerances we'll graduate to OrdinaryDiffEq.jl in a separate
# extension package; the production worked-example PKU model is non-stiff
# at physiologic timescales so RK4 with a small fixed step is fine.

"""
    rk4(f!, u0, p, t0, tend, dt) -> (ts, U)

Integrate `du/dt = f!(du, u, p, t)` from `t0` to `tend` with step `dt`.
Returns the time vector `ts` (length `N+1`) and the state matrix `U`
(`length(u0) × (N+1)`).
"""
function rk4(
    f!::Function,
    u0::AbstractVector{Float64},
    p,
    t0::Float64,
    tend::Float64,
    dt::Float64,
)
    @assert tend > t0 "tend must be > t0"
    @assert dt > 0 "dt must be positive"
    n = Int(ceil((tend - t0) / dt))
    ts = collect(range(t0; step=dt, length=n + 1))
    # Snap final time to exact tend so consumers can rely on ts[end] == tend
    ts[end] = tend
    U = zeros(Float64, length(u0), n + 1)
    U[:, 1] = u0

    k1 = similar(u0)
    k2 = similar(u0)
    k3 = similar(u0)
    k4 = similar(u0)
    utmp = similar(u0)

    @inbounds for i in 1:n
        t = ts[i]
        h = ts[i + 1] - t
        u = @view U[:, i]

        f!(k1, u, p, t)
        @. utmp = u + 0.5 * h * k1
        f!(k2, utmp, p, t + 0.5 * h)
        @. utmp = u + 0.5 * h * k2
        f!(k3, utmp, p, t + 0.5 * h)
        @. utmp = u + h * k3
        f!(k4, utmp, p, t + h)

        @. U[:, i + 1] = u + (h / 6) * (k1 + 2k2 + 2k3 + k4)
    end
    return ts, U
end
