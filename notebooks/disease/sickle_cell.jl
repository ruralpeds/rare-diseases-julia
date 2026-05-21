### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ banner
md"""
> **⚠ Not for clinical use.** Research and education only.

# Sickle cell disease — worked example

**MONDO:0011382** · Gene: **HBB** (HGNC:4827) · UniProt **P68871**.

Single canonical variant (HbS = Glu6Val). The model tracks
deoxy-HbS polymerization with hydroxyurea / gene-therapy effects
captured through HbF fraction.
"""

# ╔═╡ deps
using Catalyst
using OrdinaryDiffEq
using RDSimulation

# ╔═╡ baseline
md"## 1. Baseline vs. hydroxyurea vs. gene therapy"

function steady_polymer(hbf)
    prob = hbs_scd_problem(; hbf_fraction=hbf, tspan=(0.0, 50.0))
    sol  = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
    return sol.u[end][2]
end

baseline   = steady_polymer(0.01)
hydroxyurea = steady_polymer(0.15)
gene_tx    = steady_polymer(0.40)

md"""
Steady-state polymer fraction:

- Baseline SCD (HbF ≈ 1%): **$(round(baseline; digits=3))**
- On hydroxyurea (HbF ≈ 15%): **$(round(hydroxyurea; digits=3))**
- Post-exa-cel / lyfgenia (HbF ≈ 40%): **$(round(gene_tx; digits=3))**
"""

# ╔═╡ sweep
md"## 2. HbF dose-response sweep"

hbf_sweep = 0.0:0.05:0.6
polymer_sweep = [steady_polymer(h) for h in hbf_sweep]

# ╔═╡ Cell order:
# ╠═banner
# ╠═deps
# ╠═baseline
# ╠═sweep
