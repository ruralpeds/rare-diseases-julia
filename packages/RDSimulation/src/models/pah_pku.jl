# Phenylketonuria (PKU) — PAH enzyme kinetics, simplified.
#
# State (whole-body, lumped one-compartment plasma in mmol):
#   u[1] = Phe   (phenylalanine, blood)
#   u[2] = Tyr   (tyrosine, blood)
#
# Dynamics:
#   d[Phe]/dt = I_Phe(t) - V_PAH([Phe])
#   d[Tyr]/dt = V_PAH([Phe]) - k_clear_Tyr * [Tyr]
#
# V_PAH is Michaelis-Menten with substrate inhibition disabled; effective
# enzyme activity is scaled by `residual_activity ∈ [0, 1]` (genotype) and
# by `bh4_factor ∈ [1, +∞)` (sapropterin / BH4 cofactor supplementation).
#
# Parameters are order-of-magnitude PKU literature values and are
# explicitly approximate — production work would replace them with
# Bayesian fits to longitudinal patient cohort data.

"""
    PAHParameters

Parameters for the PAH kinetic model. Units: mmol, hour.

* `kcat`            - turnover number for wild-type PAH (1/h)
* `Km`              - Michaelis constant for Phe (mmol/L)
* `E_total`         - effective enzyme amount (arbitrary units;
                      paired with kcat to set `Vmax`)
* `residual_activity` - genotype-derived activity multiplier, ∈ [0, 1]
* `bh4_factor`      - cofactor multiplier, ∈ [1, 5]; 1.0 means no treatment
* `k_clear_Tyr`     - first-order tyrosine clearance (1/h)
* `intake_phe`      - dietary phenylalanine influx (mmol/h)
"""
Base.@kwdef mutable struct PAHParameters
    kcat::Float64              = 60.0     # 1/h
    Km::Float64                = 1.0      # mmol/L
    E_total::Float64           = 1.0      # AU
    residual_activity::Float64 = 1.0
    bh4_factor::Float64        = 1.0
    k_clear_Tyr::Float64       = 0.3      # 1/h
    intake_phe::Float64        = 2.0      # mmol/h
end

"""
    default_pah_parameters(; kwargs...) -> PAHParameters

Convenience constructor accepting overrides for any field.
"""
default_pah_parameters(; kwargs...) = PAHParameters(; kwargs...)

# Right-hand-side closure factory. Bound a parameter struct to a function
# matching the rk4 signature `f!(du, u, p, t)`.
function _pah_rhs(p::PAHParameters)
    return function (du, u, _, _t)
        phe = max(u[1], 0.0)
        tyr = max(u[2], 0.0)
        vmax = p.kcat * p.E_total * p.residual_activity * p.bh4_factor
        v = vmax * phe / (p.Km + phe)
        du[1] = p.intake_phe - v
        du[2] = v - p.k_clear_Tyr * tyr
        return nothing
    end
end

"""
    pah_residual_activity(class::Symbol) -> Float64

Lookup of canonical PAH-variant functional classes from the BIOPKU
literature, abstracted to a single residual-activity multiplier.

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
    variant_effect(p::PAHParameters, class::Symbol; bh4_factor=p.bh4_factor) -> PAHParameters

Return a copy of `p` with `residual_activity` set from `class` and
`bh4_factor` overridden if supplied.
"""
function variant_effect(p::PAHParameters, class::Symbol; bh4_factor=p.bh4_factor)
    q = deepcopy(p)
    q.residual_activity = pah_residual_activity(class)
    q.bh4_factor = bh4_factor
    return q
end

"""
    MolecularModel(PAHParameters[, name]) -> MolecularModel

Build a `MolecularModel` for the PKU/PAH ODE.
"""
function MolecularModel(p::PAHParameters, name::AbstractString="PAH-PKU")
    rhs = _pah_rhs(p)
    pdict = Dict{Symbol,Float64}(
        :kcat => p.kcat, :Km => p.Km, :E_total => p.E_total,
        :residual_activity => p.residual_activity,
        :bh4_factor => p.bh4_factor,
        :k_clear_Tyr => p.k_clear_Tyr, :intake_phe => p.intake_phe,
    )
    return MolecularModel(String(name), [:Phe, :Tyr], pdict, rhs)
end
