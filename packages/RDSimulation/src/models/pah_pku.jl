# Phenylketonuria (PKU) — PAH enzyme kinetics, on top of Catalyst.
#
# State (whole-body lumped plasma, mmol/L):
#   Phe = phenylalanine
#   Tyr = tyrosine
#
# Dynamics:
#   ∅       --> Phe   at rate intake_phe
#   Phe     --> Tyr   at rate mm(Phe, Vmax_eff, Km)
#   Tyr     --> ∅     at rate k_clear_Tyr · Tyr
#
# Vmax_eff = kcat · E_total · residual_activity · bh4_factor
#
# `residual_activity ∈ [0, 1]` is genotype-derived; `bh4_factor ≥ 1` is
# the sapropterin (BH4 cofactor) potentiation factor.
#
# Parameters are order-of-magnitude PKU literature values; production
# work should replace them with Bayesian fits to longitudinal cohort data.

"""
    PAH_PKU

`Catalyst.ReactionSystem` for the simplified PAH/PKU mass-action model.
Use with `ODEProblem(PAH_PKU, u0, tspan, p)` and any `OrdinaryDiffEq`
solver.
"""
const PAH_PKU = @reaction_network PAH_PKU begin
    @parameters intake_phe kcat E_total residual_activity bh4_factor Km k_clear_Tyr
    @species Phe(t) Tyr(t)

    intake_phe,                                            ∅   --> Phe
    mm(Phe, kcat * E_total * residual_activity * bh4_factor, Km), Phe --> Tyr
    k_clear_Tyr,                                            Tyr --> ∅
end

"""
    default_pah_parameters(; kwargs...) -> Vector{Pair{Symbol,Float64}}

Default PAH parameter vector in the symbol-pair form accepted by
`ODEProblem`. Keyword arguments override individual entries.
"""
function default_pah_parameters(;
    intake_phe::Float64=2.0,
    kcat::Float64=60.0,
    E_total::Float64=1.0,
    residual_activity::Float64=1.0,
    bh4_factor::Float64=1.0,
    Km::Float64=1.0,
    k_clear_Tyr::Float64=0.3,
)
    return [
        :intake_phe        => intake_phe,
        :kcat              => kcat,
        :E_total           => E_total,
        :residual_activity => residual_activity,
        :bh4_factor        => bh4_factor,
        :Km                => Km,
        :k_clear_Tyr       => k_clear_Tyr,
    ]
end

"""
    pah_residual_activity(class::Symbol) -> Float64

Canonical PAH-variant functional classes from the BIOPKU literature,
collapsed to a single residual-activity multiplier.

| class                | residual_activity |
|----------------------|-------------------|
| `:null`              | 0.0               |
| `:classical_pku`     | 0.02              |
| `:moderate_pku`      | 0.10              |
| `:mild_pku`          | 0.25              |
| `:mild_hpa`          | 0.50              |
| `:wildtype`          | 1.0               |
"""
function pah_residual_activity(class::Symbol)
    class === :null          && return 0.0
    class === :classical_pku && return 0.02
    class === :moderate_pku  && return 0.10
    class === :mild_pku      && return 0.25
    class === :mild_hpa      && return 0.50
    class === :wildtype      && return 1.0
    throw(ArgumentError("unknown PAH variant class: $class"))
end

"""
    variant_effect(params, class; bh4_factor=1.0) -> Vector{Pair{Symbol,Float64}}

Return a parameter vector with `residual_activity` set from `class` and
`bh4_factor` set as supplied.
"""
function variant_effect(
    params::AbstractVector{<:Pair{Symbol,Float64}},
    class::Symbol;
    bh4_factor::Float64=1.0,
)
    d = Dict(params)
    d[:residual_activity] = pah_residual_activity(class)
    d[:bh4_factor] = bh4_factor
    return [k => d[k] for k in (:intake_phe, :kcat, :E_total,
                                 :residual_activity, :bh4_factor,
                                 :Km, :k_clear_Tyr)]
end

"""
    pah_pku_problem(; variant=:wildtype, bh4=1.0, u0=[:Phe=>0.1,:Tyr=>0.1],
                     tspan=(0.0, 48.0), kwargs...) -> ODEProblem

Convenience constructor that bundles `variant_effect`, default parameters,
and initial conditions into a ready-to-`solve` `ODEProblem`.

```julia
using OrdinaryDiffEq, RDSimulation
prob = pah_pku_problem(; variant=:classical_pku, bh4=1.0)
sol  = solve(prob, Tsit5())
```
"""
function pah_pku_problem(;
    variant::Symbol=:wildtype,
    bh4::Float64=1.0,
    u0=[:Phe => 0.1, :Tyr => 0.1],
    tspan::Tuple{<:Real,<:Real}=(0.0, 48.0),
    kwargs...,
)
    p = variant_effect(default_pah_parameters(; kwargs...), variant; bh4_factor=bh4)
    return ODEProblem(PAH_PKU, u0, tspan, p)
end
