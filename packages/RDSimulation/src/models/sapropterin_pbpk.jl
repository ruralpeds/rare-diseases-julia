# Sapropterin (BH4) — minimal one-compartment PBPK for oral dosing.
#
# Simplified to a single central compartment with first-order absorption
# from the gut and first-order clearance. This is the right granularity
# to feed the PAH bh4_factor parameter in pah_pku.jl with a time-varying
# exposure rather than a static multiplier.
#
# State (mg):
#   A_gut     -- mass in gut
#   A_central -- mass in plasma compartment
#
# Dynamics:
#   gut         A_gut --> A_central        rate ka
#   A_central                              cleared first-order, rate kel
#   ∅       --> A_gut at dose times        impulse (handled by callback)
#
# Parameters chosen for order-of-magnitude alignment with the human
# sapropterin label (oral, half-life ~6h). Production work should fit
# these to published PK studies via Turing.jl.

"""
    SAPROPTERIN_PBPK

`Catalyst.ReactionSystem` for a one-compartment sapropterin PBPK. The
absorption rate `ka` and elimination rate `kel` are exposed as parameters
so callers can vary them or fit them.

```julia
prob = ODEProblem(SAPROPTERIN_PBPK, [:A_gut=>10.0, :A_central=>0.0],
                  (0.0, 24.0),
                  [:ka => 1.2, :kel => 0.115])
sol  = solve(prob, Tsit5())
```
"""
const SAPROPTERIN_PBPK = @reaction_network SAPROPTERIN_PBPK begin
    @parameters ka kel
    @species A_gut(t) A_central(t)

    ka,  A_gut     --> A_central
    kel, A_central --> ∅
end

"""
    sapropterin_pbpk_problem(; dose_mg=10.0, ka=1.2, kel=0.115,
                              tspan=(0.0, 24.0)) -> ODEProblem

One-compartment PBPK with a single oral dose at t=0. Defaults give a
half-life of ~6 hours.
"""
function sapropterin_pbpk_problem(;
    dose_mg::Float64=10.0,
    ka::Float64=1.2,
    kel::Float64=0.115,
    tspan::Tuple{<:Real,<:Real}=(0.0, 24.0),
)
    u0 = [:A_gut => dose_mg, :A_central => 0.0]
    p = [:ka => ka, :kel => kel]
    return ODEProblem(SAPROPTERIN_PBPK, u0, tspan, p)
end
