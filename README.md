# rare-diseases-julia

A Julia-native research platform for rare diseases — open and government data
only, with mechanistic simulation, diagnostic reasoning, and treatment
exploration.

> **Not for clinical use.** This is a research and educational project.
> Diagnostic and treatment outputs are hypotheses, not medical advice.

## What it does

1. **Unified knowledge graph** over HPO, MONDO, Orphanet, ClinVar, gnomAD,
   UniProt, AlphaFold, Reactome, ChEMBL, FDA orphan drug data, AACT
   (ClinicalTrials.gov), FAERS, GARD, and more.
2. **Mechanistic simulation** at three scales — molecular (ODE/reaction
   networks via SciML), physiological (PBPK / QSP), and cohort
   (agent-based via `Agents.jl`).
3. **Diagnostic reasoning** — phenotype + variant scoring (Phenomizer /
   Phrank / ACMG style).
4. **Treatment exploration** — approved + trial + repurposing candidates
   with mechanism rationale and PK/PD feasibility.

See [`PLAN.md`](./PLAN.md) for the full build plan and rationale.

## Repository layout

```
packages/
  RareDiseaseCore.jl   — IDs, base types, provenance
  RDDataSources.jl     — source registry + downloaders/parsers
  RDOntology.jl        — HPO, MONDO, ORDO, OMIM, DO graph
  RDGenomics.jl        — variants, VCF, ClinVar, gnomAD
  RDProteomics.jl      — UniProt, AlphaFold, InterPro
  RDPathways.jl        — Reactome, WikiPathways, SIGNOR, BioModels
  RDPharmacology.jl    — ChEMBL, DrugBank-open, FDA OOPD, RxNorm
  RDClinical.jl        — AACT trials, FAERS, GARD, Orphanet
  RDDiagnostics.jl     — phenotype + variant differential ranking
  RDSimulation.jl      — ODE / PBPK / ABM disease models
  RDTreatment.jl       — drug-target effect and repurposing
  RDImmunology.jl      — comprehensive immune ecosystem base types
  RDApp.jl             — REST + WebSocket service, HTMX UI
data/
  SOURCES.tsv          — canonical source registry (license, sha, version)
  manifest.toml        — pinned dataset versions
notebooks/             — Pluto notebooks per worked-example disease
```

## Quick start (dev)

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test("RareDiseaseCore")'
```

## License

Code: MIT (see [`LICENSE`](./LICENSE)).
Data: each source retains its own license; see [`data/LICENSES/`](./data/LICENSES/)
and [`data/SOURCES.tsv`](./data/SOURCES.tsv).

## Ethics

See [`ETHICS.md`](./ETHICS.md) — not-for-clinical-use stance, dual-use
considerations, and data scope.

## Citation

See [`CITATION.cff`](./CITATION.cff).
