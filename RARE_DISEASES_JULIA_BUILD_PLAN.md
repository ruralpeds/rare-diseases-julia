# rare-diseases-julia — Comprehensive Build Plan

**Repo:** ruralpeds/rare-diseases-julia
**Target:** Sprint 2 — close Sprint 1 open gaps + layer in v2 systems biology + add a pediatric rare disease module library
**Stack:** Julia (monorepo; 12 in-tree packages with a new 13th `RDMetabolic` + a new 14th `RDMultiscale`)
**Buildable by:** GitHub Copilot agents in a devcontainer pinned to `ghcr.io/ruralpeds/sprint-base:v1`, following `ruralpeds/.github` sprint standard v1 and this repo's existing AGENTS conventions (`CLAUDE.md`, `GEMINI.md`).

---

## 0. Why this design

The repo's original `PLAN.md` (Sprint 1) established a twelve-package monorepo covering rare-disease data, ontologies, genomics, proteomics, pathways, pharmacology, clinical, diagnostics, simulation, treatment, and an Oxygen.jl REST surface. As of the current commit, most of the API surfaces and several worked-example diseases (PKU, CF, sickle cell) are in place, but a clearly-bounded set of "not started" rows in §9 of `PLAN.md` represent the remaining real-data ingestion work, and the simulation engine — while functional — does not yet exercise the deeper Julia systems-biology stack (COBREXA, fuller CellML/SBML coverage) nor a broad pediatric rare-disease module library.

This sprint is built around three principles:

First, **close before extending.** Every "not started" row in `PLAN.md` §9 maps directly onto a phase in this sprint. The variant-aware diagnostic ranker cannot be honest until ClinVar and gnomAD are wired; the treatment ranker cannot account for PK/PD feasibility until pharmacology data is fully ingested. We close those gaps first.

Second, **lean into the Julia ecosystem.** Per the existing repo's ecosystem-first policy in `CLAUDE.md`, this sprint adds COBREXA.jl for constraint-based metabolic modeling and CellMLToolkit.jl for importing the hundreds of curated, peer-reviewed mechanistic models from the CellML and BioModels repositories. The goal is to make every pediatric rare disease module a thin domain-specific composition over well-validated upstream models rather than reinventing biology in this repo.

Third, **scale across diseases via fan-out.** Phases 08 and 09 each fan out across a fixed set of pediatric rare diseases (priority IEMs, channelopathies, lysosomal storage disorders). The fan-out items are deliberately tractable for an autonomous Copilot agent in one work-day each, with shared templates so the agents are not improvising structural decisions per disease.

The sprint deliberately stops short of adopting the full multi-organ digital twin from the broader PedNeoTwin v2 plan. Phase 10 provides a stub adapter layer (`RDMultiscale`) so a downstream `PedNeoTwin.jl` repo can consume this repo's molecular trajectories without this repo growing organ-level physiology code.

---

## 1. Repository layout

After this sprint completes, the repository tree looks like:

```
rare-diseases-julia/
├── packages/
│   ├── RareDiseaseCore/        (existing — no changes this sprint)
│   ├── RDDataSources/          (Phase 00: Parquet writer; Phase 03–05: new sources)
│   ├── RDOntology/             (existing — no changes this sprint)
│   ├── RDGenomics/             (Phase 01: VCF + ClinVar + gnomAD)
│   ├── RDProteomics/           (Phase 02: UniProt + AlphaFold + InterPro)
│   ├── RDPathways/             (Phase 03: Reactome + WikiPathways + SIGNOR)
│   ├── RDPharmacology/         (Phase 04: ChEMBL + DrugBank-open + FDA OOPD)
│   ├── RDClinical/             (Phase 04: AACT + FAERS)
│   ├── RDDiagnostics/          (Phase 05: variant-aware rank_diagnoses)
│   ├── RDTreatment/            (Phase 05: PK/PD feasibility)
│   ├── RDSimulation/           (Phase 06: CellML/SBML deep; Phase 08–09: disease modules)
│   ├── RDMetabolic/            (Phase 07: NEW — COBREXA wrapper + tissue reconstructions + IEM perturbations)
│   ├── RDMultiscale/           (Phase 10: NEW — adapter contract for downstream physiology)
│   └── RDApp/                  (Phase 11: route additions for new modules)
├── data/
│   ├── SOURCES.tsv             (Phase 01–05 additions)
│   ├── LICENSES/               (Phase 01–05 additions)
│   └── manifest.toml           (pinned dataset versions)
├── notebooks/disease/
│   ├── phenylketonuria.jl      (existing)
│   ├── cystic_fibrosis.jl      (existing)
│   ├── sickle_cell.jl          (existing)
│   ├── glycogen_storage_Ia.jl  (Phase 08)
│   ├── pompe.jl                (Phase 08)
│   ├── otc_deficiency.jl       (Phase 08)
│   ├── msud.jl                 (Phase 08)
│   ├── mcad_deficiency.jl      (Phase 08)
│   ├── mito_complex_I.jl       (Phase 08)
│   ├── galactosemia.jl         (Phase 08)
│   ├── long_qt_1.jl            (Phase 09)
│   ├── long_qt_2.jl            (Phase 09)
│   ├── dravet.jl               (Phase 09)
│   ├── gaucher.jl              (Phase 09)
│   ├── niemann_pick_C.jl       (Phase 09)
│   ├── mps_I.jl                (Phase 09)
│   └── krabbe.jl               (Phase 09)
├── docs/                       (Phase 11: expanded Documenter site)
├── scripts/
│   ├── dev_bootstrap.jl        (existing)
│   ├── refresh_data.jl         (existing — extended Phase 01–05)
│   ├── demo.jl                 (existing — extended Phase 11)
│   └── module_smoke_test.jl    (Phase 08: NEW)
├── test/
│   ├── runtests.jl             (existing — extended)
│   └── disease_modules_test.jl (Phase 08–09: NEW)
├── PLAN.md                                     (existing strategic plan)
├── RARE_DISEASES_JULIA_BUILD_PLAN.md            (this file — Sprint 2 plan)
└── ... (existing governance files unchanged)
```

Two new packages introduced this sprint: `RDMetabolic` (constraint-based metabolic modeling) and `RDMultiscale` (adapter contract for downstream consumers). Everything else extends existing packages.

---

## 2. Design Principles for This Sprint

Carry over the principles already codified in `PLAN.md` §0 and `CLAUDE.md`, then add three sprint-specific ones:

The first sprint-specific principle is **upstream-model-first.** Whenever a CellML or SBML model exists for a phenomenon (e.g., O'Hara-Rudy ventricular action potential for LQT modeling, Bergman minimal model for glucose-insulin), we import and validate rather than re-deriving. The new `RDSimulation` curated model registry tags each imported model with `:pediatric_validated`, `:adult_only`, or `:needs_modification` so downstream module authors can pick safely.

The second is **disease module template uniformity.** Phases 08 and 09 follow an identical template: each disease module is a Julia file at `packages/RDSimulation/src/diseases/<slug>.jl` exporting a `build_<slug>_model(; ...)` function that returns an `RDDiseaseModel` struct, with a paired Pluto notebook at `notebooks/disease/<slug>.jl` and a paired test file at `packages/RDSimulation/test/diseases/test_<slug>.jl`. This uniformity is what makes Copilot agent fan-out actually tractable.

The third is **provenance for every imported model.** Every CellML/SBML import gets a `DataProvenance` stamp identical in form to the data-source provenance the repo already uses. A model whose provenance cannot be established does not enter the registry.

---

## 3. Dependency Graph for This Sprint

Phases 01–04 are formally serial in the plan (each depends on the previous merged) but their work is largely independent and can run in close succession. Phase 05 (reasoning closure) genuinely needs Phases 01–04 merged because variant-aware diagnostics requires the full data layer.

Phase 06 (CellML/SBML deep integration) is independent of Phases 01–05 in content but depends on Phase 00 (Parquet writer) for the model registry's persistence.

Phase 07 (COBREXA) depends on Phase 06 because the IEM perturbation library cross-references imported metabolic SBML models.

Phases 08 and 09 (disease modules, fan-out) depend on Phases 06 and 07 because every module composes imported models with metabolic perturbations.

Phase 10 (multiscale adapter) is independent of Phases 08–09 in content but depends on Phase 07 because the adapter contract references the metabolic trajectory types.

Phase 11 (docs + release) depends on everything.

---

## 4. Key Julia Packages Added or Promoted This Sprint

The existing `CLAUDE.md` already documents the ecosystem-first stack. This sprint adds or more deeply uses:

`COBREXA.jl` becomes the constraint-based metabolic modeling backbone in the new `RDMetabolic` package. Used for flux balance analysis (FBA), flux variability analysis (FVA), parsimonious FBA (pFBA), and a dynamic-FBA scaffold against the time-course glucose/lactate/ammonia variables that the existing PKU and CF models already track.

`CellMLToolkit.jl` becomes the bridge to the CellML model repository, exposing roughly 1,500 published mechanistic models. The new `RDSimulation/src/imports/cellml.jl` module exposes a typed import path with a curated subset registry.

`SBMLToolkit.jl` (already partially in the repo for the BioModels loader) is promoted to first-class status with a deeper curated subset registry covering the BioModels Database entries most relevant to pediatric rare disease.

`VariantCallFormat.jl` and `GeneticVariation.jl` enter `RDGenomics` for VCF parsing.

`Parquet2.jl` and `Arrow.jl` enter `RDDataSources` for the processed-data writer.

`JET.jl` and `Aqua.jl` are enforced uniformly across all packages in CI by Phase 00.

---

## 5. Risks specific to this sprint

The CellML and SBML import registry is the highest-leverage element of the sprint but also the most fragile: upstream model heterogeneity means some imports will fail structural simplification or numerical integration with default solvers. Phase 06 mitigates this with explicit registry tagging (`:needs_modification`) and per-model smoke tests; agents must be empowered to tag a model `:needs_modification` and skip it rather than try to fix upstream biology themselves.

The fan-out in Phases 08 and 09 has fifteen disease modules total. The seeder enforces a maximum of eight simultaneous tasks per the org standard, so the fan-out runs in two waves. The dependency between Phase 08 and Phase 09 is on the merged Phase 08 work, not on every Phase-8 module individually, so partial Phase-8 completion does not block Phase 09 from starting.

The `RDMetabolic` package introduces COBREXA.jl as a heavy new dependency. Phase 07 must lock COBREXA.jl version in the meta `Project.toml` to avoid resolver-induced version churn breaking other packages.

The variant-aware diagnostic ranker in Phase 05 requires meaningful test fixtures with realistic VCFs. Curating those fixtures from public Synthea data plus deliberately hand-crafted Mendelian cases is non-trivial; Phase 01 includes the fixture authoring.

---

## 6. Out of scope for this sprint

Multi-organ physiological simulation (the v2 PedNeoTwin organ-level twin) remains out of scope. Phase 10 provides the adapter contract only; the actual organ-level twin lives in a future `ruralpeds/PedNeoTwin.jl` repo.

LLM-backed natural-language interfaces over the diagnostic and treatment surfaces remain out of scope.

EHR integration (FHIR R4 client, HL7 v2 listener) remains out of scope and is reserved for a downstream clinical-integration sprint.

Real-time data assimilation against bedside monitor streams remains out of scope.

---

## 7. Phased build

### Phase 00 — Foundation refresh and Parquet writer

1. Update `PLAN.md` §9 Implementation Progress table to reflect any drift since the last commit; mark this Sprint 2 starting state.
2. Add `RARE_DISEASES_JULIA_BUILD_PLAN.md` (this file) to the repo root if not already present from sprint bootstrap.
3. In `packages/RDDataSources`, add a Parquet/Arrow processed-data writer with the existing `Source` trait, using `Parquet2.jl` for the canonical processed-data format and `Arrow.jl` as an alternative reader path.
4. Add an integration test that round-trips a small fixture through the writer and reads it back with both Parquet2 and Arrow.
5. Add Aqua.jl + JET.jl quality gates uniformly to all twelve existing package test suites with shared helpers under `test/_qa_helpers.jl`.
6. Bump the meta `Project.toml` to add the new transitive dependencies, run `Pkg.resolve()`, commit the updated `Manifest.toml`.

Done check: `julia --project=. -e 'using Pkg; Pkg.instantiate(); include("test/runtests.jl")'` exits 0 with all twelve existing package test suites green and a new Parquet round-trip test reported.

### Phase 01 — Genomics closure: VCF, ClinVar, gnomAD

1. In `packages/RDGenomics`, add a VCF reader using `VariantCallFormat.jl` with a typed `Variant` struct that interoperates with the existing HGVS parser and `ACMGEvidence` accumulator.
2. Add a ClinVar ingestion path: download the weekly ClinVar VCF and `variant_summary.txt.gz` via the existing `RDDataSources` HTTP downloader, parse into a typed `ClinVarRecord` Parquet table, and add a sha256 + manifest entry.
3. Add a gnomAD GraphQL client returning typed `AlleleFrequency` records with a local SQLite/Parquet cache. Default cache policy `:cache_first`.
4. Wire ClinVar and gnomAD into the ACMG evidence accumulator so PS1, PM1, PM2, PM5, BA1, BS1, BS2 can be auto-populated from the data layer where applicable.
5. Add Synthea-generated and small public fixtures under `packages/RDGenomics/test/fixtures/` so CI does not depend on the network.

Depends on: Phase 00 merged.
Done check: `julia --project=packages/RDGenomics -e 'using Pkg; Pkg.test()'` passes including new ClinVar/gnomAD ACMG integration test.

### Phase 02 — Proteomics closure: UniProt, AlphaFold, InterPro

1. In `packages/RDProteomics`, add a UniProt Swiss-Prot ingestion path: stream the XML from the existing `RDDataSources` downloader and write a typed Parquet table keyed by `UniProtAcc`.
2. Add an on-demand AlphaFold structure fetcher (per-accession PDB pull) with a content-addressed local cache, sha256 verified.
3. Add InterPro domain annotation ingestion via the InterPro REST API with the existing throttle/retry plumbing.
4. Expose `lookup_protein(::UniProtAcc)` returning a unified `ProteinRecord` combining UniProt metadata, AlphaFold structure handle, and InterPro domains.
5. Add a structural-impact helper using `BioStructures.jl` that reports solvent accessibility and secondary-structure context for a residue position — feeds into Phase 05's variant-aware ranker.

Depends on: Phase 01 merged.
Done check: `julia --project=packages/RDProteomics -e 'using Pkg; Pkg.test()'` passes including a round-trip test that pulls UniProt + AlphaFold + InterPro for a fixture accession and reports a structural-impact score for a fixture residue.

### Phase 03 — Pathways closure: Reactome, WikiPathways, SIGNOR

1. In `packages/RDPathways`, add a Reactome ingestion path consuming Reactome's SBML and BioPAX exports; convert into the existing `MetaGraphsNext`-backed `PathwayNetwork`.
2. Add a WikiPathways ingestion path consuming WikiPathways' GPML exports.
3. Add a SIGNOR ingestion path consuming SIGNOR's tab-separated signed-edge export, preserving the activation/inhibition sign on each edge.
4. Expose a unified `load_pathway_corpus(; sources=[:reactome, :wikipathways, :signor])` constructor that merges all three into a single `PathwayNetwork` with provenance preserved per edge.
5. Validate the merged network's connectivity and reachability metrics in a new integration test.

Depends on: Phase 02 merged.
Fan out across: Reactome, WikiPathways, SIGNOR.
Done check: `julia --project=packages/RDPathways -e 'using Pkg; Pkg.test()'` passes including the merged-corpus connectivity test.

### Phase 04 — Pharmacology and Clinical data closure

1. In `packages/RDPharmacology`, add a ChEMBL SQLite ingestion path that loads the latest ChEMBL release into a local SQLite file, plus a typed `drugs_for_target(::UniProtAcc)` query returning ranked candidate compounds with activity thresholds.
2. Add a DrugBank-open subset loader (the openly-licensed structures and identifiers only).
3. Add an FDA OOPD (Orphan Drug Designation) loader from the FDA's public CSV export, joining by RxCUI / ChEMBL ID where possible.
4. In `packages/RDClinical`, add an AACT (ClinicalTrials.gov AACT relational snapshot) loader using the daily PostgreSQL dump or its CSV mirror; ingest the `studies`, `conditions`, `interventions`, and `eligibilities` tables.
5. Add a FAERS quarterly file loader for adverse event signals, with per-drug aggregation helpers.
6. Add small fixture databases for both `RDPharmacology` and `RDClinical` so CI remains offline.

Depends on: Phase 03 merged.
Done check: `julia --project=. -e 'include("packages/RDPharmacology/test/runtests.jl"); include("packages/RDClinical/test/runtests.jl")'` exits 0 with new ChEMBL, OOPD, AACT, and FAERS smoke tests reported.

### Phase 05 — Reasoning closure: variant-aware diagnostics and PK/PD feasibility

1. In `packages/RDDiagnostics`, extend `rank_diagnoses` to consume a VCF (or already-parsed `Vector{Variant}`) and produce an Exomiser-style integrated score combining HPO semantic similarity (already in place) with ACMG-evidence-weighted variant pathogenicity (newly available from Phase 01).
2. Add `inheritance_pattern_fit(::Disease, ::Vector{Variant}, ::Pedigree)` to weight candidate diseases by Mendelian inheritance pattern fit; require at minimum a `MinimalPedigree` (proband-only) interface.
3. In `packages/RDTreatment`, add a PK/PD feasibility scorer that consumes the existing one-compartment PBPK template plus pharmacology data (Phase 04) to estimate whether a candidate drug can reach its predicted target site at therapeutic concentrations; surface this as a tiered evidence score.
4. Update the not-for-clinical-use banner constant to remain visible in every variant-aware and PK/PD-aware output path.
5. Add an end-to-end integration test: a fixture proband with HPO terms and a VCF resolves to the expected ClinVar-curated diagnosis in the top-five ranked candidates, with at least one ranked treatment surfacing a PK/PD feasibility score.

Depends on: Phase 04 merged.
Done check: `julia --project=. -e 'include("test/integration/variant_aware_and_pkpd_test.jl")'` exits 0 with all assertions passing.

### Phase 06 — CellML and SBML deep integration

1. In `packages/RDSimulation`, formalize the imports subdirectory with `src/imports/cellml.jl` and `src/imports/sbml.jl`. The existing BioModels SBML loader is refactored into this layout without behavior change.
2. Add a curated `ModelRegistry` in `packages/RDSimulation/data/model_registry.toml` enumerating roughly 30 CellML and 20 SBML imports relevant to pediatric rare disease, each tagged with `:pediatric_validated`, `:adult_only`, `:needs_modification`, plus provenance fields (source, version, sha256, original publication PMID).
3. Implement `import_cellml(::AbstractString)` and `import_sbml(::AbstractString)` that return a `RegisteredModel` struct exposing the underlying `ODESystem` (via ModelingToolkit), the provenance, and the registry tags.
4. Add smoke tests that import every entry in the registry, perform a structural simplification, and run a single time integration over a default time span without numerical failure.
5. Pin upstream CellMLToolkit.jl and SBMLToolkit.jl versions in the meta `Project.toml` with explicit `[compat]` bounds.

Depends on: Phase 00 merged.
Done check: `julia --project=packages/RDSimulation -e 'using RDSimulation; for m in RDSimulation.list_registry(); RDSimulation.smoke_test(m); end'` exits 0 with all 50 registered models passing their smoke test.

### Phase 07 — Constraint-based metabolic modeling: RDMetabolic.jl

1. Create the new `packages/RDMetabolic` package with its own `Project.toml`, `src/RDMetabolic.jl`, `test/runtests.jl`, and dev-in entry in the meta `Project.toml`.
2. Wrap `COBREXA.jl` with a typed `MetabolicModel` struct holding a flux balance model plus annotations (tissue tag, organism, source, sha256). Expose `flux_balance(::MetabolicModel)`, `flux_variability(::MetabolicModel)`, and `parsimonious_flux_balance(::MetabolicModel)`.
3. Add four tissue-specific reconstructions derived from Recon3D: `liver`, `brain`, `muscle`, `neonatal_hepatic`. Each is a `MetabolicModel` instance with curated exchange-flux boundary conditions documented inline.
4. Implement an `IEMPerturbation` library: for each priority IEM in Phase 08, a function `apply_iem!(::MetabolicModel, ::Symbol; severity::Float64)` that bounds the relevant enzymatic flux. Initial library covers GSD Ia (G6PC), Pompe (GAA), OTC (OTC), MSUD (BCKDHA/B/E2/E3), MCAD (ACADM), Mitochondrial complex I (NDUFS subunits), Galactosemia (GALT).
5. Implement a `dynamic_fba` scaffold that integrates an FBA snapshot against the existing time-course glucose/lactate/ammonia state variables in `RDSimulation`.
6. Add validation tests confirming that an unperturbed `neonatal_hepatic` model produces plausible glucose output flux and that each `IEMPerturbation` reduces its target flux as expected.

Depends on: Phase 06 merged.
Done check: `julia --project=packages/RDMetabolic -e 'using Pkg; Pkg.test()'` passes with all tissue reconstructions and all seven Phase-7 IEM perturbations validated.

### Phase 08 — Pediatric IEM module library

1. Author one disease module per fan-out item below, each following the standing template: a Julia source file at `packages/RDSimulation/src/diseases/<slug>.jl` exporting `build_<slug>_model(; severity, age_days, weight_kg)` returning an `RDDiseaseModel`; a Pluto notebook at `notebooks/disease/<slug>.jl` walking through model build, perturbation, and a 48-hour simulation; a test file at `packages/RDSimulation/test/diseases/test_<slug>.jl` asserting clinically recognizable trajectory features (e.g., GSD Ia fasting hypoglycemia trajectory, OTC hyperammonemia onset, MSUD elevated branched-chain amino acids).
2. For each module, compose: an `RDMetabolic.IEMPerturbation` from Phase 07; any directly relevant Catalyst.jl reaction network for affected signaling; and any pediatric-validated CellML/SBML import from the Phase 06 registry where applicable.
3. For each module, register the disease against the existing `MONDO` and `OMIM` identifier maps in `RDOntology` so the diagnostic ranker can surface it.
4. Each module's notebook concludes with a treatment exploration cell calling `rank_treatments` against the disease and surfacing approved + repurposing candidates with PK/PD feasibility (from Phase 05).
5. Add a collective `test/disease_modules_test.jl` harness that loads all eight Phase-8 disease modules, runs each smoke test, and reports aggregate timings.

Depends on: Phase 07 merged.
Fan out across: glycogen_storage_Ia, pompe, otc_deficiency, msud, mcad_deficiency, mito_complex_I, phenylketonuria_v2, galactosemia.
Done check: `julia --project=. -e 'include("test/disease_modules_test.jl")'` exits 0 with all eight Phase-8 disease modules' smoke tests passing.

### Phase 09 — Pediatric channelopathy and lysosomal storage module library

1. Author one disease module per fan-out item below, following the same standing template as Phase 08, with two adaptations: channelopathy modules import the relevant pediatric-validated cardiac or cortical CellML model from the Phase 06 registry (ten Tusscher or O'Hara-Rudy for cardiac LQT modules; an appropriate cortical mean-field model for Dravet) and apply the variant-specific kinetic-parameter perturbation; LSD modules use a Catalyst.jl reaction network for the affected hydrolase pathway plus a slow-accumulation state variable.
2. For each channelopathy module, expose a `simulate_action_potential` helper returning a typed `APRecord` with computed QT interval (for cardiac) or a seizure-burden surrogate (for cortical).
3. For each LSD module, expose an `accumulation_trajectory` helper returning the projected substrate accumulation over 90 simulated days with and without enzyme-replacement-therapy parameters from Phase 04 pharmacology data.
4. Register each module's disease against `MONDO` and `OMIM` maps; ensure the diagnostic ranker surfaces them given representative HPO phenotype sets.
5. Extend `test/disease_modules_test.jl` to load all seven Phase-9 modules in addition to the eight from Phase 08.

Depends on: Phase 08 merged.
Fan out across: long_qt_1, long_qt_2, dravet, gaucher, niemann_pick_C, mps_I, krabbe.
Done check: `julia --project=. -e 'include("test/disease_modules_test.jl")'` exits 0 with all fifteen Phase-8 + Phase-9 disease modules' smoke tests passing.

### Phase 10 — Multi-scale adapter contract: RDMultiscale.jl

1. Create the new `packages/RDMultiscale` package with its own `Project.toml`, `src/RDMultiscale.jl`, `test/runtests.jl`, and dev-in entry in the meta `Project.toml`. This package is intentionally minimal in this sprint.
2. Define abstract types `MolecularTrajectory`, `PhysiologyParameterPatch`, `MultiscaleAdapter` plus an interface protocol that downstream consumers (notably a future `PedNeoTwin.jl`) implement to receive molecular-network outputs and emit organ-level parameter perturbations.
3. Implement a reference `RDPhysiologyAdapter` that, given a `MolecularTrajectory` from any Phase-8 or Phase-9 disease module, emits a `PhysiologyParameterPatch` documenting the intended effect on organ-level variables (e.g., GSD Ia → reduced hepatic glucose output rate; LQT1 → prolonged QT-equivalent parameter). The adapter does not itself simulate physiology; it documents the contract.
4. Publish an `interface.md` in `packages/RDMultiscale/docs/` specifying the contract for downstream physiology consumers, including the unit conventions and provenance requirements.
5. Add an integration test demonstrating round-trip: a Phase-8 disease module produces a `MolecularTrajectory`, the adapter emits a `PhysiologyParameterPatch`, and a mock physiology consumer in the test suite consumes the patch and reports a documented expected delta.

Depends on: Phase 07 merged.
Done check: `julia --project=packages/RDMultiscale -e 'using Pkg; Pkg.test()'` passes with the round-trip integration test reported.

### Phase 11 — Docs site refresh, worked-diseases gallery, v0.1.0 release

1. Expand the Documenter site with one page per new disease module from Phases 08–09 (fifteen new pages), each linking to its Pluto notebook and to its model registry entries.
2. Add a new "Systems Biology" docs section covering RDMetabolic, the CellML/SBML registry, and the multi-scale adapter contract.
3. Update the home page worked-diseases gallery to include screenshots and short rationales for all new modules.
4. Refresh `PLAN.md` §9 Implementation Progress table to mark every Sprint-2 row as done.
5. Bump version to `0.1.0` in the meta `Project.toml` and in each subpackage's `Project.toml`.
6. Create a `CHANGELOG.md` entry describing Sprint 2 deliverables; tag the release `v0.1.0`.
7. Verify Documenter build via the existing `Docs.yml` workflow and that the deployed site renders correctly.

Depends on: Phase 10 merged.
Done check: `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'` exits 0 and the deployed Documenter site lists every Sprint-2 disease module under the worked-diseases gallery.

---

## 8. Acceptance criteria

The sprint is considered complete when all of the following hold simultaneously on `main`:

The repository CI passes with `julia --project=. -e 'using Pkg; Pkg.instantiate(); include("test/runtests.jl")'` exiting 0 on Ubuntu, macOS, and Windows runners with Julia LTS and Julia current.

Every package in `packages/` has Aqua.jl and JET.jl gates green, including the two new packages `RDMetabolic` and `RDMultiscale`.

All fifteen disease modules from Phases 08 and 09 pass their smoke tests under `test/disease_modules_test.jl`, and the variant-aware diagnostic ranker from Phase 05 returns the expected diagnosis in the top five candidates for at least twelve of the fifteen modules' fixture probands.

The Documenter site deploys and lists every new module page; the worked-diseases gallery on the home page shows the full set; `CHANGELOG.md` records `v0.1.0`; the `v0.1.0` tag exists on `main`.

`PLAN.md` §9 reflects all Sprint-2 rows as done; no Sprint-1 rows that were "not started" at the start of this sprint remain "not started"; any new "not started" rows added during this sprint reflect forward work intentionally deferred to a future sprint.

---

*End of Sprint 2 build plan. Maintain in `ruralpeds/rare-diseases-julia/RARE_DISEASES_JULIA_BUILD_PLAN.md`.*
