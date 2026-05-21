### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ banner
md"""
> **⚠ Not for clinical use.** Research and education only.

# Phenylketonuria (PKU) — worked example

**MONDO:0009861** · **ORPHA:716** · **OMIM:261600**

Gene: **PAH** (HGNC:8582) → UniProt **P00439**.

This notebook exercises the full platform stack against PKU:

1. Build a PAH/PKU reaction network with `Catalyst`.
2. Solve for several variant classes with `OrdinaryDiffEq`.
3. Overlay a sapropterin (BH4) PBPK template.
4. Run a small Agents.jl cohort and tally outcomes.
"""

# ╔═╡ deps
using Catalyst
using OrdinaryDiffEq
using Agents
using Random

using RDSimulation

# ╔═╡ variants
md"## 1. Variant class → steady-state [Phe]"

function steady_phe(class; bh4=1.0, tend=96.0)
    prob = pah_pku_problem(; variant=class, bh4=bh4, tspan=(0.0, tend))
    sol  = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
    return sol.u[end][1]  # final Phe; species[1] in PAH_PKU
end

variant_table = [
    (:wildtype,        steady_phe(:wildtype)),
    (:mild_hpa,        steady_phe(:mild_hpa)),
    (:mild_pku,        steady_phe(:mild_pku)),
    (:moderate_pku,    steady_phe(:moderate_pku)),
    (:classical_pku,   steady_phe(:classical_pku)),
    (:null,            steady_phe(:null)),
]

# ╔═╡ sapropterin
md"## 2. Sapropterin (BH4) effect on residual-activity variants"

bh4_off = steady_phe(:classical_pku; bh4=1.0)
bh4_on  = steady_phe(:classical_pku; bh4=2.0)
md"""
Classical PKU steady-state Phe:
- without sapropterin: $(round(bh4_off; digits=3)) mmol/L
- with sapropterin (bh4_factor=2): $(round(bh4_on; digits=3)) mmol/L
"""

# ╔═╡ pbpk
md"## 3. Sapropterin PBPK — single oral dose"

pbpk_sol = solve(sapropterin_pbpk_problem(; dose_mg=10.0, tspan=(0.0, 24.0)),
                 Tsit5())

# ╔═╡ cohort
md"""
## 4. Cohort outcome simulation

50 classical-PKU patients on diet, small per-tick probability of metabolic
crisis. The hook `agent_step!` is where genotype-specific PK/PD will
plug in once the variant-aware diagnostic + treatment pipelines land.
"""

function patient_step!(agent, model)
    if agent.state == :on_diet && rand(abmrng(model)) < 0.02
        agent.state = :metabolic_crisis
    end
end

cohort = build_cohort_model(;
    n_agents=50,
    initial_state=:on_diet,
    agent_step! = patient_step!,
    rng=Xoshiro(0),
)
run_cohort!(cohort, 100)

ncrisis = count(a -> a.state == :metabolic_crisis, allagents(cohort))
md"After 100 ticks, $(ncrisis) / 50 patients reached metabolic crisis."

# ╔═╡ Cell order:
# ╠═banner
# ╠═deps
# ╠═variants
# ╠═sapropterin
# ╠═pbpk
# ╠═cohort
