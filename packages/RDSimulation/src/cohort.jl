# Patient cohort built on Agents.jl.

"""
    Patient

Agent type for cohort simulations. `state` is a free-form `Symbol`
(e.g. `:on_diet`, `:on_sapropterin`, `:metabolic_crisis`). `covariates`
holds genotype/age/etc.
"""
@agent struct Patient(NoSpaceAgent)
    state::Symbol
    covariates::Dict{Symbol,Any}
end

"""
    build_cohort_model(; n_agents, initial_state=:on_diet,
                        agent_step!,
                        rng::AbstractRNG=Random.default_rng()) -> AgentBasedModel

Construct an `Agents.StandardABM` with `n_agents` patients in the
specified `initial_state` and the given per-tick `agent_step!`
transition function.

The `agent_step!` callback receives `(agent, model)` and is expected to
mutate `agent.state` in place. Use `abmrng(model)` for the seeded RNG.
"""
function build_cohort_model(;
    n_agents::Int,
    initial_state::Symbol=:on_diet,
    agent_step!::Function,
    rng::AbstractRNG=Random.default_rng(),
)
    model = StandardABM(Patient; agent_step!=agent_step!, rng=rng)
    for _ in 1:n_agents
        add_agent!(model, initial_state, Dict{Symbol,Any}())
    end
    return model
end

"""
    run_cohort!(model, n_steps) -> model

Advance the cohort `n_steps` ticks. Returns the same model so calls
chain. Use `Agents.allagents(model)` to inspect post-run states.
"""
function run_cohort!(model, n_steps::Integer)
    run!(model, n_steps)
    return model
end
