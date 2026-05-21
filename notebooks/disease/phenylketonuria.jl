### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 00000001
md"""
# Phenylketonuria (PKU) — worked example

**MONDO:0009861** · **ORPHA:716** · **OMIM:261600**

Gene: **PAH** (HGNC:8582) → UniProt **P00439** → AlphaFold structure.

Pipeline (to be filled in as phases land):

1. Pull disease record from MONDO + Orphanet + GARD.
2. Pull PAH variants from ClinVar; annotate with gnomAD AF and MANE
   transcript HGVS.
3. Map each variant onto P00439 residues; overlay pLDDT, domains, and
   active-site distance from AlphaFold.
4. Load BioModels PAH kinetic model; perturb `kcat`/`Km` for selected
   variants.
5. Run cohort simulation (Agents.jl) of diet vs sapropterin.
6. Recover known treatments via `rank_treatments(MondoId("MONDO:0009861"))`.

> Not for clinical use — research example only.
"""

# ╔═╡ Cell order:
# ╠═00000001
