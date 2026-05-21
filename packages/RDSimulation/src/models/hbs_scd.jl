# Sickle cell disease (SCD) — HbS polymerization, simplified.
#
# HbS (HBB:Glu6Val) polymerizes when deoxygenated, producing the long
# fibers that distort red cells. Hydroxyurea and gene-modifying therapies
# raise fetal hemoglobin (HbF), which dilutes deoxy-HbS and inhibits
# polymerization (Eaton-Hofrichter formalism).
#
# State (dimensionless):
#   polymer_fraction ∈ [0, 1]
#
# Dynamics:
#   d(polymer)/dt = k_form * (1 - polymer) * effective_deoxyHbS^n
#                 - k_diss * polymer
#
#   effective_deoxyHbS = deoxyHbS_total * (1 - hbf_inhibition * hbf_fraction)
#
# Eaton-Hofrichter has n≈10 (extreme cooperativity); we use n=4 here
# so RK-class explicit solvers don't burn time on the steep nonlinearity
# — qualitative behavior is preserved.

"""
    HBS_SCD

`Catalyst.ReactionSystem` for the polymerization of deoxy-HbS into
fiber polymer. `hbf_fraction` is the input that hydroxyurea (and
gene-therapy regimens) modulate; `hbf_inhibition` is how strongly HbF
dilutes the effective polymerizable pool.

Two species:
  * `Mono` — non-polymerized fraction (= 1 - polymer)
  * `Poly` — polymer fraction
"""
const HBS_SCD = @reaction_network HBS_SCD begin
    @parameters k_form k_diss deoxyHbS n hbf_fraction hbf_inhibition
    @species Mono(t) Poly(t)

    k_form * (deoxyHbS * (1 - hbf_inhibition * hbf_fraction))^n, Mono --> Poly
    k_diss,                                                        Poly --> Mono
end

"""
    hbs_scd_problem(; deoxyHbS=1.0, hbf_fraction=0.01,
                     hbf_inhibition=2.0,
                     k_form=1.0, k_diss=0.5, n=4.0,
                     u0=[:Mono=>1.0, :Poly=>0.0],
                     tspan=(0.0, 5.0)) -> ODEProblem

Build an `ODEProblem` for the simplified HbS polymerization model.

```julia
# Baseline SCD (low HbF)
prob_lo  = hbs_scd_problem(; hbf_fraction=0.01)
sol_lo   = solve(prob_lo, Tsit5())

# On hydroxyurea (HbF ~15%) — polymer fraction drops substantially
prob_hi  = hbs_scd_problem(; hbf_fraction=0.15)
sol_hi   = solve(prob_hi, Tsit5())
```
"""
function hbs_scd_problem(;
    deoxyHbS::Float64=1.0,
    hbf_fraction::Float64=0.01,
    hbf_inhibition::Float64=2.0,
    k_form::Float64=1.0,
    k_diss::Float64=0.5,
    n::Float64=4.0,
    u0=[:Mono => 1.0, :Poly => 0.0],
    tspan::Tuple{<:Real,<:Real}=(0.0, 5.0),
)
    p = [
        :k_form         => k_form,
        :k_diss         => k_diss,
        :deoxyHbS       => deoxyHbS,
        :n              => n,
        :hbf_fraction   => hbf_fraction,
        :hbf_inhibition => hbf_inhibition,
    ]
    return ODEProblem(HBS_SCD, u0, tspan, p)
end
