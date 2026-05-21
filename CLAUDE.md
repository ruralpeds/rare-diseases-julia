# CLAUDE.md — project context for Claude Code (and other agents)

This file gives Claude Code, Claude Code on the web, and any agent that
reads project memory files at startup the conventions for working on
this repository. `GEMINI.md` is a copy with the same content.

## What this repo is

A Julia-native research platform for rare diseases — open and government
data only, with mechanistic simulation, diagnostic reasoning, and
treatment exploration. See `PLAN.md` for the full design and phase plan.

**Not for clinical use.** Every diagnostic / treatment surface must
carry the not-for-clinical-use banner.

## Layout

Monorepo. Twelve in-tree Julia packages under `packages/`. Each has
its own `Project.toml`, `src/`, `test/`. The meta `Project.toml` at the
repo root is the dev environment.

```
packages/RareDiseaseCore   — IDs, base types, provenance
packages/RDDataSources     — Source trait, downloaders, cache + manifest
packages/RDOntology        — HPO/MONDO/Orphanet graph, IC similarity
packages/RDGenomics        — variants, HGVS, ACMG evidence
packages/RDProteomics      — UniProt, AlphaFold, BioStructures
packages/RDPathways        — Graphs.jl + MetaGraphsNext, Guney 2016
packages/RDPharmacology    — ChEMBL, DrugBank-open, FDA OOPD
packages/RDClinical        — AACT trials, FAERS, GARD, Orphanet
packages/RDDiagnostics     — phenotype + variant differential ranking
packages/RDSimulation      — Catalyst + OrdinaryDiffEq + Agents.jl
packages/RDTreatment       — proximity-based candidate ranking
packages/RDApp             — Oxygen.jl REST/WebSocket service
```

## Ecosystem-first policy

Use existing Julia packages instead of reinventing. See `PLAN.md` §8.6
for the full table. Highlights:

- ODE/SDE → `OrdinaryDiffEq` (`Tsit5`, etc.)
- Reaction networks → `Catalyst`
- Symbolic systems → `ModelingToolkit`
- Agent-based → `Agents.jl`
- Graphs → `Graphs.jl` + `MetaGraphsNext`
- Bio I/O → `FASTX`, `BioStructures`
- HTTP / routing → `Oxygen` on `HTTP`
- Bayesian → `Turing.jl`
- Sensitivity → `GlobalSensitivity`

Hand-rolled only where no clear ecosystem winner exists: OBO parsing,
HGVS parsing, PBPK templates, Guney 2016 proximity.

## Conventions

1. **Identifier types.** Use `RareDiseaseCore` types (`HPOId`, `MondoId`,
   `UniProtAcc`, …) instead of raw `String` for any namespaced ID.
2. **Provenance.** Every persisted record carries a `DataProvenance`
   stamp (source, version, sha256, retrieved_at).
3. **Source registry gate.** Every downloader must annotate
   `# source: <name>` matching a row in `data/SOURCES.tsv`. CI enforces.
4. **Not-for-clinical-use banner.** `RDDiagnostics` and `RDTreatment`
   outputs always carry the banner constant.
5. **Test fixtures, not network.** Use small bundled fixtures
   (`packages/RDOntology/test/fixtures/mini_hpo.obo`) so CI is offline.

## Running things

```bash
# Dev environment (all packages dev'd in)
julia --project=. scripts/dev_bootstrap.jl

# Test one package
julia --project=packages/RareDiseaseCore -e 'using Pkg; Pkg.test()'

# Co-load every package
julia --project=. -e 'include("test/runtests.jl")'

# Launch the API server (once domain logic is wired in)
julia --project=packages/RDApp -e 'using RDApp; RDApp.start_server()'
```

## Worked-example diseases

PKU, cystic fibrosis, sickle cell. Notebooks under `notebooks/disease/`.
PKU has the most plumbing today.

## When working on this repo

- Prefer editing existing files to creating new ones.
- Don't add new top-level docs unless asked; the user prefers in-code
  docstrings and `PLAN.md` updates.
- Don't add network calls to tests.
- Don't claim a phase is "done" without a runnable example.
- Keep `PLAN.md` §9 (Implementation Progress) updated when phases move.
