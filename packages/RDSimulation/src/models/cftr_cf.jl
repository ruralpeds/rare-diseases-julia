# Cystic fibrosis (CF) — CFTR-mediated airway surface liquid (ASL) model.
#
# Reduced lumped one-state ODE: ASL volume on the apical surface of an
# airway epithelial unit. CFTR-mediated Cl⁻ and water secretion balances
# ENaC-mediated Na⁺/water absorption. Loss of CFTR drops ASL volume,
# which drives the clinical mucus-clearance failure.
#
# State (µL per epithelial unit, arbitrary scaling):
#   ASL
#
# Dynamics:
#   ∅       --> ASL  at rate secretion_base · CFTR_activity
#   ASL     --> ∅    at rate absorption_base · ASL
#
# CFTR_activity = class_factor · potentiator_factor · corrector_factor
#
# Class factors collapsed from the CFTR variant classification:
#   I   (no protein synthesis)         -> 0.00
#   II  (misfolding, e.g. F508del)     -> 0.02
#   III (gating defect, e.g. G551D)    -> 0.08
#   IV  (residual function)            -> 0.20
#   V   (reduced amount)               -> 0.30
#   VI  (turnover defect)              -> 0.40
#   wildtype                           -> 1.00
#
# Modulator multipliers:
#   ivacaftor   (potentiator)   — multiplies class III and gating
#                                 sub-population of class IV by ~3x
#   tezacaftor  (corrector)     — multiplies class II by ~2x
#   elexacaftor (next-gen corrector) — multiplies class II by ~4x
#                                       (trikafta combination effect)

"""
    CFTR_CF

`Catalyst.ReactionSystem` for the simplified CFTR / ASL volume model.
"""
const CFTR_CF = @reaction_network CFTR_CF begin
    @parameters secretion_base absorption_base cftr_activity
    @species ASL(t)

    secretion_base * cftr_activity,     ∅   --> ASL
    absorption_base,                    ASL --> ∅
end

"""
    cftr_class_factor(class::Symbol) -> Float64

Activity multiplier for a CFTR variant class. Throws on unknown class.

| class           | factor |
|-----------------|--------|
| `:I`            | 0.00   |
| `:II`           | 0.02   |
| `:III`          | 0.08   |
| `:IV`           | 0.20   |
| `:V`            | 0.30   |
| `:VI`           | 0.40   |
| `:wildtype`     | 1.00   |
"""
function cftr_class_factor(class::Symbol)
    class === :I        && return 0.00
    class === :II       && return 0.02
    class === :III      && return 0.08
    class === :IV       && return 0.20
    class === :V        && return 0.30
    class === :VI       && return 0.40
    class === :wildtype && return 1.00
    throw(ArgumentError("unknown CFTR variant class: $class"))
end

"""
    cftr_modulator_factor(; class, ivacaftor=false, tezacaftor=false,
                           elexacaftor=false) -> Float64

Composite multiplier for the CFTR-modulator combinations on the market.
Calibrated for qualitative behavior, not pharmacokinetic accuracy.

* Potentiators (`ivacaftor`) only help classes with functional protein
  on the membrane: class III and class IV.
* Correctors (`tezacaftor`, `elexacaftor`) rescue misfolded protein:
  class II.
* Trikafta = elexacaftor + tezacaftor + ivacaftor.
"""
function cftr_modulator_factor(;
    class::Symbol,
    ivacaftor::Bool=false,
    tezacaftor::Bool=false,
    elexacaftor::Bool=false,
)
    f = 1.0
    if ivacaftor && (class === :III || class === :IV)
        f *= 3.0
    end
    if class === :II
        tezacaftor  && (f *= 2.0)
        elexacaftor && (f *= 2.0)   # multiplicative on top of tezacaftor
        ivacaftor   && (f *= 1.5)   # marginal gain after correction
    end
    return f
end

"""
    cftr_cf_problem(; class, modulators=(;), secretion_base=1.0,
                     absorption_base=0.5,
                     u0=[:ASL=>0.0], tspan=(0.0, 24.0)) -> ODEProblem

Build an `ODEProblem` for the CF ASL model. `modulators` is a NamedTuple
with optional `ivacaftor`, `tezacaftor`, `elexacaftor` booleans.

```julia
# F508del/F508del on Trikafta
prob = cftr_cf_problem(; class=:II,
    modulators=(; tezacaftor=true, elexacaftor=true, ivacaftor=true))
sol  = solve(prob, Tsit5())
```
"""
function cftr_cf_problem(;
    class::Symbol,
    modulators::NamedTuple=(;),
    secretion_base::Float64=1.0,
    absorption_base::Float64=0.5,
    u0=[:ASL => 0.0],
    tspan::Tuple{<:Real,<:Real}=(0.0, 24.0),
)
    cls = cftr_class_factor(class)
    mod = cftr_modulator_factor(;
        class=class,
        ivacaftor   = get(modulators, :ivacaftor,   false),
        tezacaftor  = get(modulators, :tezacaftor,  false),
        elexacaftor = get(modulators, :elexacaftor, false),
    )
    p = [
        :secretion_base => secretion_base,
        :absorption_base => absorption_base,
        :cftr_activity => cls * mod,
    ]
    return ODEProblem(CFTR_CF, u0, tspan, p)
end
