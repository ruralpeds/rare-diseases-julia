# Worked-example diseases

Three diseases anchor the simulation work. Each has a Catalyst reaction
network in `RDSimulation`, a Pluto notebook under `notebooks/disease/`,
and integration coverage in the test suite.

!!! warning "Not for clinical use"
    These models are intentionally simplified for research and
    education. Parameters are order-of-magnitude only; production
    work should refit them to longitudinal cohort data with
    `Turing.jl`-based Bayesian inference.

## Phenylketonuria (PKU)

**MONDO:0009861** · **ORPHA:716** · **OMIM:261600** ·
Gene: PAH (HGNC:8582) · UniProt P00439

Two-species plasma ODE (Phe, Tyr) with Michaelis-Menten PAH activity:

```julia
using OrdinaryDiffEq, RDSimulation
prob = pah_pku_problem(; variant=:classical_pku, bh4=2.0)
sol  = solve(prob, Tsit5())
```

Variant classes from BIOPKU (`:null`, `:classical_pku`,
`:moderate_pku`, `:mild_pku`, `:mild_hpa`, `:wildtype`) collapse to a
single `residual_activity` multiplier. Sapropterin (BH4) is captured by
`bh4_factor`; the standalone one-compartment PK is `SAPROPTERIN_PBPK`.

Notebook: `notebooks/disease/phenylketonuria.jl`.

## Cystic fibrosis (CF)

**MONDO:0009061** · Gene: CFTR (HGNC:1884) · UniProt P13569

Lumped ASL-volume ODE balancing CFTR-driven secretion against
ENaC-driven absorption. Variant-class factor × modulator factor
determines effective CFTR activity.

```julia
prob = cftr_cf_problem(; class=:II,
    modulators=(; tezacaftor=true, elexacaftor=true, ivacaftor=true))
sol  = solve(prob, Tsit5())
```

Modulator rules:

| Modulator   | Mechanism   | Helps classes |
|-------------|-------------|---------------|
| Ivacaftor   | Potentiator | III, IV       |
| Tezacaftor  | Corrector   | II            |
| Elexacaftor | Corrector   | II            |

Class I (no protein) gets no benefit from any modulator.

Notebook: `notebooks/disease/cystic_fibrosis.jl`.

## Sickle cell disease (SCD)

**MONDO:0011382** · Gene: HBB (HGNC:4827) · UniProt P68871

Monomer ↔ polymer interconversion with cooperative polymerization
(softened to n=4 for numerical tractability). HbF fraction (raised by
hydroxyurea or by gene-modifying therapies such as exa-cel) dilutes the
effective polymerizable pool.

```julia
prob = hbs_scd_problem(; hbf_fraction=0.15)   # on hydroxyurea
sol  = solve(prob, Tsit5())
```

Notebook: `notebooks/disease/sickle_cell.jl`.
