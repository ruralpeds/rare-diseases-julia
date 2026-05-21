# Lightweight Markov-style cohort scaffold.
#
# Agents.jl will be the production engine for Phase 9.3; this module
# provides a hand-rolled Markov-chain cohort runner so that worked
# examples (PKU diet vs sapropterin) can land without pulling in the
# Agents.jl dependency tree.

"""
    CohortAgent

One simulated patient. `state` is a free-form symbol (e.g. `:on_diet`,
`:on_sapropterin`, `:metabolic_crisis`). `covariates` is a small dict for
genotype/age/etc.
"""
mutable struct CohortAgent
    id::Int
    state::Symbol
    covariates::Dict{Symbol,Any}
    history::Vector{Tuple{Float64,Symbol}}
end
CohortAgent(id, state) = CohortAgent(id, state, Dict{Symbol,Any}(), [(0.0, state)])

"""
    cohort_simulate(agents, transition!; tspan, dt=1.0, rng) -> agents

Run a discrete-time Markov-style cohort simulation. `transition!` is
called once per agent per tick with signature
`transition!(agent, t, dt, rng) -> nothing` and is responsible for
sampling the next state.
"""
function cohort_simulate(
    agents::Vector{CohortAgent},
    transition!::Function;
    tspan::Tuple{<:Real,<:Real},
    dt::Real=1.0,
    rng::AbstractRNG=Random.default_rng(),
)
    t = Float64(tspan[1])
    tend = Float64(tspan[2])
    while t < tend
        step = min(Float64(dt), tend - t)
        for a in agents
            transition!(a, t, step, rng)
            push!(a.history, (t + step, a.state))
        end
        t += step
    end
    return agents
end
