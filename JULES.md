# JULES.md — starter file for Jules (jules.google) on `ruralpeds/rare-diseases-julia`

## How to use this file

**If you're a human kicking off a Jules session:** open [jules.google.com](https://jules.google.com), point Jules at this repo, and either let Jules read this file from the repo automatically, or paste the **Starter task** block below into Jules' task input box.

**If you're Jules opening this repo:** this is your operational memory. Read this file first, then the strategic context, then your specific task.

## Read order for Jules

1. **JULES.md** (this file) — conventions, constraints, session scope rules
2. **CLAUDE.md** / **GEMINI.md** — identical project context written for sibling agents; useful background
3. **PLAN.md** — strategic plan and Sprint 1 history (where the repo came from)
4. **RARE_DISEASES_JULIA_BUILD_PLAN.md** — active Sprint 2 plan, twelve phases (where the repo is going)
5. The assigned GitHub issue, when present — your specific task scope

---

## Starter task — Phase 00: Foundation refresh and Parquet writer

> **For the human:** paste the entire fenced block below into Jules' task input to kick off Sprint 2 Phase 00. Jules will create a branch, execute the six steps, run the Done check, and open a PR.

```
Open RARE_DISEASES_JULIA_BUILD_PLAN.md in the repository ruralpeds/rare-diseases-julia. Locate the section "### Phase 00 — Foundation refresh and Parquet writer" under "## 7. Phased build". Execute the six numbered steps in that phase exactly as written. Do not start any work from Phase 01 or beyond in the same pull request.

Phase 00 deliverables:
1. Update PLAN.md §9 Implementation Progress table to reflect the current state of the repo and mark this commit as the Sprint 2 starting point. Do not change the historical rows; just append a Sprint 2 marker.
2. Confirm that RARE_DISEASES_JULIA_BUILD_PLAN.md exists at the repo root (it was added in PR #4). If missing, copy it from PR #4's diff.
3. In packages/RDDataSources, implement a Parquet/Arrow processed-data writer using Parquet2.jl and Arrow.jl that integrates with the existing Source trait. Expose a typed write_processed(::Source, ::DataFrame; format=:parquet) API plus a matching read_processed helper. Persist files under data/processed/<source>/<version>.parquet.
4. Add an integration test under packages/RDDataSources/test/ that round-trips a small fixture DataFrame through the writer and reads it back via both Parquet2 and Arrow readers, asserting schema and value equality.
5. Add Aqua.jl and JET.jl quality gates uniformly across all twelve existing package test suites using shared helpers under test/_qa_helpers.jl. Keep gates non-fatal for now (warn-only) so a single legacy issue doesn't block the sprint; later phases can promote to fatal.
6. Bump the meta Project.toml to add the new transitive dependencies (Parquet2, Arrow, Aqua, JET), run Pkg.resolve(), and commit the updated Manifest.toml.

Done check (must pass before opening the PR):
julia --project=. -e 'using Pkg; Pkg.instantiate(); include("test/runtests.jl")'

PR conventions:
- Branch name: agent/phase-00/foundation-refresh
- Target: main
- PR title: "Phase 00 — Foundation refresh and Parquet writer"
- PR description: include the Done check command and its full output (last 50 lines is fine), plus a one-line summary per numbered step.

Constraints (must hold):
- No new package dependencies beyond those four (Parquet2, Arrow, Aqua, JET). If you discover another dependency is needed, surface it in the PR description and stop; do not silently add it.
- No real PHI in fixtures. No network calls in tests. Banner constants preserved.
- USA-only citations if any documentation is added.
- Do not edit JULES.md, CLAUDE.md, GEMINI.md, or RARE_DISEASES_JULIA_BUILD_PLAN.md in this PR.
```

---

## What this repo is

A Julia-native research platform for rare diseases. **Open and government data only.** Mechanistic simulation at three scales (molecular reaction networks via SciML; physiological PBPK/QSP; cohort agent-based via `Agents.jl`), diagnostic reasoning (Phenomizer/Phrank/ACMG-style), and treatment exploration (approved, trial, repurposing candidates with mechanism rationale and PK/PD feasibility).

**Not for clinical use.** Every diagnostic and treatment surface must carry the not-for-clinical-use banner.

## Layout

Monorepo. Twelve in-tree Julia packages today; two more (`RDMetabolic`, `RDMultiscale`) added during Sprint 2. Each package has its own `Project.toml`, `src/`, `test/`. The meta `Project.toml` at the repo root is the dev environment.

```
packages/
  RareDiseaseCore  — IDs, base types, provenance
  RDDataSources    — Source trait, downloaders, cache, manifest
  RDOntology       — HPO/MONDO/Orphanet graph, IC similarity
  RDGenomics       — variants, HGVS, ACMG evidence
  RDProteomics     — UniProt, AlphaFold, BioStructures
  RDPathways       — Graphs.jl + MetaGraphsNext, Guney 2016
  RDPharmacology   — ChEMBL, DrugBank-open, FDA OOPD
  RDClinical       — AACT trials, FAERS, GARD, Orphanet
  RDDiagnostics    — phenotype + variant differential ranking
  RDSimulation     — Catalyst + OrdinaryDiffEq + Agents.jl
  RDTreatment      — proximity-based candidate ranking
  RDApp            — Oxygen.jl REST + WebSocket service
  RDMetabolic      — (Sprint 2 Phase 07) COBREXA + tissue reconstructions + IEM perturbations
  RDMultiscale     — (Sprint 2 Phase 10) adapter contract for future PedNeoTwin.jl
```

## Constraints (non-negotiable for every Jules session)

These apply universally:

1. **USA-only citations** for medical and clinical references — AAP, AHA, NRP, CDC, NICHD, ACOG, PubMed USA-authored journals. International sources allowed only when scientifically indispensable (Orphanet for nomenclature; UniProt / PDB / AlphaFold for protein data; Ensembl for variant annotation). Flag international sources inline when used.
2. **Open and government data sources only.** No proprietary databases. ChEMBL is OK. DrugBank only the openly-licensed subset. KEGG only academic use with explicit notation.
3. **No new dependencies without explicit approval.** If a phase requires a package not already in the meta `Project.toml`, surface it in the PR description and stop work; do not silently add it.
4. **No real PHI** in any fixture, test, or example data. Use Synthea synthetic data only.
5. **No network calls in tests.** Use bundled fixtures under `packages/<Name>/test/fixtures/`. CI must run offline.
6. **Not-for-clinical-use banner** preserved on every `RDDiagnostics` and `RDTreatment` output path. The banner constant lives in `RareDiseaseCore`.

## Ecosystem-first policy

Use the documented Julia packages instead of reinventing. Full table in `PLAN.md` §8.6. Highlights:

| Need | Use |
|---|---|
| ODE / SDE solvers | `OrdinaryDiffEq` (Tsit5 default, Rodas5P for stiff) |
| Reaction networks | `Catalyst.jl` |
| Symbolic systems | `ModelingToolkit.jl` |
| Agent-based modeling | `Agents.jl` |
| Graphs | `Graphs.jl` + `MetaGraphsNext.jl` |
| Bio I/O | `FASTX.jl`, `BioStructures.jl`, `VariantCallFormat.jl` |
| HTTP / routing | `Oxygen.jl` on `HTTP.jl` |
| Bayesian inference | `Turing.jl` |
| Sensitivity analysis | `GlobalSensitivity.jl` |
| Constraint-based metabolic | `COBREXA.jl` (new in Sprint 2) |
| CellML model import | `CellMLToolkit.jl` |
| SBML model import | `SBMLToolkit.jl` |
| Tabular persistence | `Parquet2.jl`, `Arrow.jl` (added Phase 00) |
| Quality gates | `Aqua.jl`, `JET.jl` (added Phase 00) |

Hand-rolled code only where no clear ecosystem winner exists: OBO parsing, HGVS parsing, PBPK templates, Guney 2016 proximity.

## Conventions

1. **Identifier types.** Use `RareDiseaseCore` types (`HPOId`, `MondoId`, `UniProtAcc`, etc.) instead of raw `String` for any namespaced ID.
2. **Provenance.** Every persisted record carries a `DataProvenance` stamp (source, version, sha256, retrieved_at).
3. **Source registry gate.** Every downloader is annotated with `# source: <name>` matching a row in `data/SOURCES.tsv`. CI enforces this.
4. **Banner constant.** `RDDiagnostics` and `RDTreatment` outputs always carry the not-for-clinical-use banner.
5. **Test fixtures, not network.** Bundle small fixtures (`packages/RDOntology/test/fixtures/mini_hpo.obo` is the canonical example) so CI is offline.

## How to scope a Jules session

**One phase per session.** For fan-out phases (Phase 03 across pathway sources; Phases 08 + 09 across disease modules), **one fan-out item per session**.

Per-session checklist:

1. Read your assigned phase's numbered steps in `RARE_DISEASES_JULIA_BUILD_PLAN.md`.
2. Create a feature branch named `agent/phase-NN/<short-slug>` (e.g., `agent/phase-00/foundation-refresh`, `agent/phase-08/glycogen-storage-Ia`).
3. Execute the numbered steps in order. Do not reorder, skip, or add steps without surfacing it in the PR description.
4. Run the phase's Done check command and capture its output.
5. Open a PR against `main`. Include the Done check command and its full output (last 50 lines suffices) in the PR description, plus a one-line summary per numbered step.
6. **Stop.** Do not start the next phase in the same PR.

## Files Jules MUST NOT edit

- `audit-log/` (when present)
- `dhf/risk/` (when present)
- `policies/rulesets/` (when present)
- `.github/workflows/copilot-task-guardrails.yml` (when present)
- `JULES.md`, `CLAUDE.md`, `GEMINI.md` — request human review for any change to operational docs
- `RARE_DISEASES_JULIA_BUILD_PLAN.md` — the plan is authoritative; if a step is wrong, surface it in the PR description rather than silently editing the plan

## Files Jules SHOULD update during work

- `PLAN.md` §9 Implementation Progress table — mark your phase's row as `done` when the Done check passes
- `CHANGELOG.md` — append your phase's deliverables when one exists (it will be created in Phase 11)
- `docs/` pages corresponding to any new module or package you add

## When you're stuck

If a step in `RARE_DISEASES_JULIA_BUILD_PLAN.md` is ambiguous, do the most conservative thing that satisfies the Done check, add a `# TODO(jules)` comment in the code with your interpretation, and surface the ambiguity explicitly in the PR description so a human can correct course before merge.

If a CellML or SBML model import fails structural simplification or numerical integration (relevant to Phases 06, 07, 08, 09): tag the model `:needs_modification` in the model registry and move on. Do not try to fix upstream biology in this repo — that's outside scope.

If a Done check passes but you suspect the implementation is incomplete or fragile: still open the PR, but call out residual concerns in the PR description rather than expanding scope unilaterally.

If a step requires a new dependency you weren't authorized to add: stop, push what you have, and open a draft PR with a `needs-approval` label explaining the proposed dependency and why.

## Running things locally (parity check)

```bash
# Dev environment (all packages dev'd in)
julia --project=. scripts/dev_bootstrap.jl

# Test one package
julia --project=packages/RareDiseaseCore -e 'using Pkg; Pkg.test()'

# Co-load every package
julia --project=. -e 'include("test/runtests.jl")'

# Launch the API server (after Sprint 1 domain logic is wired)
julia --project=packages/RDApp -e 'using RDApp; RDApp.start_server()'
```

## Active sprint at a glance

Sprint 2 (`RARE_DISEASES_JULIA_BUILD_PLAN.md`), twelve phases:

```
00 — Foundation refresh + Parquet writer                ← START HERE
01 — VCF + ClinVar + gnomAD ingestion                   (depends on 00)
02 — UniProt + AlphaFold + InterPro                     (depends on 01)
03 — Reactome + WikiPathways + SIGNOR  [fan-out × 3]    (depends on 02)
04 — ChEMBL + FDA OOPD + AACT + FAERS                   (depends on 03)
05 — Variant-aware diagnostics + PK/PD feasibility      (depends on 04)
06 — CellML + SBML deep integration                     (depends on 00)
07 — RDMetabolic.jl (NEW): COBREXA + IEM perturbations  (depends on 06)
08 — Pediatric IEM modules         [fan-out × 8]        (depends on 07)
09 — Pediatric channelopathy + LSD [fan-out × 7]        (depends on 08)
10 — RDMultiscale.jl (NEW): adapter contract            (depends on 07)
11 — Docs refresh + v0.1.0 release                      (depends on 10)
```

---

*End of JULES.md. When you finish a phase, append a one-line entry to PLAN.md §9 with the phase number, deliverable, and `done`.*
