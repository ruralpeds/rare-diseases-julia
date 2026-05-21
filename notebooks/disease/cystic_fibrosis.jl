### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ banner
md"""
> **⚠ Not for clinical use.** Research and education only.

# Cystic fibrosis — worked example

**MONDO:0009061** · Gene: **CFTR** (HGNC:1884) · UniProt **P13569**.

Tests variant-class-specific treatment logic: ivacaftor for gating
mutations (class III), Trikafta (elexacaftor/tezacaftor/ivacaftor) for
F508del (class II).
"""

# ╔═╡ deps
using Catalyst
using OrdinaryDiffEq
using RDSimulation

# ╔═╡ classes
md"## 1. CFTR variant classes and modulator response"

function steady_asl(class; mods=(;))
    prob = cftr_cf_problem(; class=class, modulators=mods, tspan=(0.0, 48.0))
    sol  = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
    return sol.u[end][1]
end

untreated_table = [
    :I        => steady_asl(:I),
    :II       => steady_asl(:II),
    :III      => steady_asl(:III),
    :IV       => steady_asl(:IV),
    :V        => steady_asl(:V),
    :wildtype => steady_asl(:wildtype),
]

# ╔═╡ treatments
md"## 2. Targeted therapies recover ASL volume by class"

ivacaftor_g551d = steady_asl(:III; mods=(; ivacaftor=true))
trikafta_f508   = steady_asl(:II;
    mods=(; tezacaftor=true, elexacaftor=true, ivacaftor=true))
class_i_modulators = steady_asl(:I; mods=(; ivacaftor=true))

md"""
- G551D (class III) on ivacaftor:                  **$(round(ivacaftor_g551d; digits=3))**
- F508del/F508del (class II) on Trikafta:          **$(round(trikafta_f508; digits=3))**
- Class I (nonsense) on any modulator (no benefit): **$(round(class_i_modulators; digits=3))**
"""

# ╔═╡ Cell order:
# ╠═banner
# ╠═deps
# ╠═classes
# ╠═treatments
