# rare-diseases-julia

A Julia-native research platform for rare diseases — open and government
data only, with mechanistic simulation, diagnostic reasoning, and
treatment exploration.

!!! warning "Not for clinical use"
    This is a research and educational project. Diagnostic and
    treatment outputs are hypotheses, not medical advice.

## What it does

1. **Unified knowledge graph** over HPO, MONDO, Orphanet, ClinVar,
   gnomAD, UniProt, AlphaFold, Reactome, ChEMBL, FDA orphan drug data,
   AACT (ClinicalTrials.gov), FAERS, GARD, and more.
2. **Mechanistic simulation** at three scales — molecular (ODE /
   reaction networks via `Catalyst` + `OrdinaryDiffEq`), physiological
   (PBPK / QSP via `ModelingToolkit`), and cohort (agent-based via
   `Agents.jl`).
3. **Diagnostic reasoning** — phenotype + variant scoring
   (Phenomizer / Phrank / ACMG style).
4. **Treatment exploration** — approved + trial + repurposing
   candidates with mechanism rationale (Guney 2016 network proximity)
   and PK/PD feasibility.

## Packages

| Package | Purpose |
|---|---|
| [`RareDiseaseCore`](packages/raredisease_core.md) | IDs, base types, provenance |
| [`RDDataSources`](packages/rd_data_sources.md) | Source registry, downloaders, cache |
| [`RDOntology`](packages/rd_ontology.md) | HPO/MONDO/Orphanet graph, IC similarity |
| [`RDGenomics`](packages/rd_genomics.md) | Variants, HGVS, ACMG evidence |
| [`RDProteomics`](packages/rd_proteomics.md) | UniProt, AlphaFold, BioStructures |
| [`RDPathways`](packages/rd_pathways.md) | Graphs.jl networks, Guney 2016 proximity |
| [`RDPharmacology`](packages/rd_pharmacology.md) | ChEMBL, DrugBank-open, FDA OOPD |
| [`RDClinical`](packages/rd_clinical.md) | AACT, FAERS, GARD, Orphanet |
| [`RDDiagnostics`](packages/rd_diagnostics.md) | Differential ranking |
| [`RDSimulation`](packages/rd_simulation.md) | Catalyst + Agents.jl |
| [`RDTreatment`](packages/rd_treatment.md) | Proximity-based ranking |
| [`RDApp`](packages/rd_app.md) | Oxygen.jl REST/WebSocket service |

## Get started

```bash
git clone https://github.com/ruralpeds/rare-diseases-julia.git
cd rare-diseases-julia
make dev
make test-core
```

See [the build plan](plan.md) for the phased roadmap.
