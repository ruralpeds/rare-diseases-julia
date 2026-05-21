# Rare Diseases Julia — Build Plan

A long-form, opinionated build plan for a Julia-native platform that ingests open
clinical, genomic, proteomic, and pharmacologic data on rare diseases, and uses it
to (a) reason about diagnostic features and (b) simulate disease mechanisms,
progression, and candidate treatments.

The plan is structured in phases. Each phase has goals, data sources, Julia
package choices, deliverables, and exit criteria. Stick to the phase boundaries —
they exist so that early phases produce a usable artifact even if later phases
slip.

---

## 0. Guiding Principles

1. **Open data only.** Every primary source must be redistributable or freely
   re-derivable. Track licensing in `data/SOURCES.tsv` with one row per source
   (name, URL, license, version, retrieved-at, sha256). No source enters the
   ingestion pipeline without a row here.
2. **Julia first, FFI when needed.** Use native Julia for data, modeling, and
   simulation. Wrap C/C++/Rust tools (HTSlib, Open Babel, RDKit, PLINK2) via
   `Libdl`/`CxxWrap`/`PythonCall` only when no Julia equivalent exists.
3. **Reproducibility.** Every dataset version is content-addressed (sha256) and
   pinned in `Project.toml` + `data/manifest.toml`. Every model run emits a
   `RunManifest` with code git-sha, data hashes, RNG seed, and solver tolerances.
4. **Clinical safety stance.** This is a research and education tool. Diagnostic
   and treatment outputs are explicitly labeled "not for clinical use," gated
   behind a banner in any UI, and every recommendation surface cites primary
   evidence.
5. **Small core, fat ecosystem.** The monorepo hosts a small set of stable core
   packages and a wider set of disease-specific or model-specific add-ons that
   can move independently.

---

## 1. Repository & Package Layout

Monorepo, multiple registered Julia packages, one top-level `Project.toml` for
the development environment.

```
rare-diseases-julia/
├── Project.toml                # dev meta-environment
├── Manifest.toml
├── PLAN.md
├── README.md
├── LICENSE                     # MIT for code
├── data/
│   ├── SOURCES.tsv             # canonical source registry
│   ├── manifest.toml           # pinned versions + sha256
│   ├── raw/                    # gitignored, content-addressed
│   ├── interim/                # gitignored
│   └── processed/              # parquet/arrow, gitignored
├── packages/
│   ├── RareDiseaseCore.jl      # types, IDs, ontology graph
│   ├── RDDataSources.jl        # downloaders + parsers per source
│   ├── RDOntology.jl           # HPO, MONDO, ORDO, OMIM, MeSH bridges
│   ├── RDGenomics.jl           # variants, VCF, ClinVar, gnomAD
│   ├── RDProteomics.jl         # UniProt, AlphaFold, InterPro
│   ├── RDPathways.jl           # Reactome, KEGG-ish (open), WikiPathways
│   ├── RDPharmacology.jl       # DrugBank-open, ChEMBL, RxNorm, FDA OOPD
│   ├── RDClinical.jl           # GARD, Orphanet, ClinicalTrials.gov, FDA AERS
│   ├── RDDiagnostics.jl        # phenotype matching, differential ranking
│   ├── RDSimulation.jl         # ODE/SDE/ABM disease models
│   ├── RDTreatment.jl          # PK/PD, drug-target effect models
│   └── RDApp.jl                # Genie/Oxygen web UI + REST API
├── notebooks/                  # Pluto + Jupyter, narrative analyses
├── scripts/                    # one-off ETL and benchmarking
├── test/                       # cross-package integration tests
└── .github/workflows/          # CI (test, docs, data-build smoke)
```

Each `packages/*` directory is its own Julia package with `Project.toml`,
`src/`, `test/`, and `docs/`. They are added to the meta-environment via
`dev ./packages/<Name>`.

---

## 2. Data Sources

All sources below are open or government-released. Rows here go into
`data/SOURCES.tsv` verbatim. Group by domain.

### 2.1 Disease & Phenotype Ontologies
- **HPO** (Human Phenotype Ontology) — phenotypes, disease-phenotype links.
  `obo` + annotation TSVs. https://hpo.jax.org
- **MONDO** — unified disease ontology, cross-refs to OMIM/Orphanet/DOID.
- **Orphanet / ORDO** — rare-disease-specific, with epidemiology and
  inheritance. Orphadata XML dumps.
- **OMIM** — open subset via mim2gene; full OMIM requires registration, treat
  as optional.
- **DO (Disease Ontology)** — alternative cross-refs.
- **MeSH** — for literature linking.

### 2.2 Government / Clinical Registries
- **GARD** (NIH Genetic and Rare Diseases Information Center) — disease
  summaries.
- **ClinicalTrials.gov** — AACT relational dump (open) for trials by condition.
- **FDA Orphan Drug Designations** (OOPD) — designated/approved orphan drugs.
- **FDA FAERS** — adverse-event reports, quarterly ASCII dumps.
- **FDA Drugs@FDA + Orange Book** — approval and exclusivity data.
- **DailyMed (NLM)** — structured product labels (SPL XML).
- **RxNorm / RxClass (NLM)** — drug normalization.

### 2.3 Genetics & Genomics
- **ClinVar** — variant-disease assertions, monthly VCF + XML.
- **gnomAD** — population allele frequencies (v4 genomes/exomes).
- **dbSNP, dbVar** — variant identifiers and structural variants.
- **GTEx** — tissue expression (open summary stats).
- **GENCODE / Ensembl / RefSeq** — gene models, GTF/GFF3.
- **HGNC** — gene symbols and aliases.
- **PanelApp (Genomics England, AU)** — gene panels per condition.
- **MANE Select** — canonical transcripts.

### 2.4 Proteins & Structures
- **UniProt / SwissProt** — canonical proteins, XML + FASTA.
- **AlphaFold DB** — predicted structures by UniProt accession.
- **PDB** — experimental structures (mmCIF).
- **InterPro / Pfam** — domains.
- **STRING** (open subset) — protein-protein interactions.

### 2.5 Pathways & Networks
- **Reactome** — pathways, BioPAX/SBML/SBGN exports.
- **WikiPathways** — community pathways, GPML.
- **SIGNOR** — signed signaling network.
- **BioModels** — curated SBML models, many disease-relevant.

### 2.6 Chemistry & Pharmacology
- **DrugBank-open** (limited fields) + **ChEMBL** for bioactivity.
- **PubChem** — compounds, bioassays.
- **BindingDB** — affinities.
- **PharmGKB** (open subset) — pharmacogenomics.

### 2.7 Literature
- **PubMed / PMC** — abstracts (PubMed Baseline) + OA full text (PMC OAS).
- **Europe PMC** — alternative API.
- **bioRxiv / medRxiv** — preprints via API.

Versioning: every source download is stored under
`data/raw/<source>/<ISO8601-date>/` and registered in `data/manifest.toml`.

---

## 3. Phase Roadmap

Twelve phases. Each is roughly one to three weeks of focused work; the back
half can be parallelized across contributors.

### Phase 1 — Bootstrap (Week 1)

Goal: a minimal repo that other phases can build on.

- Create monorepo layout from §1.
- Set up CI (`.github/workflows/CI.yml`) running `Pkg.test` on Julia 1.10 LTS
  and current stable.
- Add Documenter.jl for each package, deploy to GitHub Pages.
- Add pre-commit (JuliaFormatter, end-of-file fixer).
- Write `RareDiseaseCore.jl` with:
  - Strongly-typed identifiers: `HPOId`, `MondoId`, `OrphaId`, `OmimId`,
    `HgncId`, `EnsemblGeneId`, `UniProtAcc`, `RxCui`, `ClinVarId`, `PubMedId`,
    `Pmid`. Each is a `struct` wrapping a validated string.
  - `Disease`, `Phenotype`, `Gene`, `Variant`, `Protein`, `Drug`, `Trial`
    abstract types and concrete records.
  - `DataProvenance` struct attached to every record (source, version, sha256,
    retrieved_at).

Exit: `]test RareDiseaseCore` green; identifiers parse/round-trip.

### Phase 2 — Data Layer (Weeks 2–3)

Goal: deterministic, resumable ingestion.

- `RDDataSources.jl` exposes a `Source` trait with:
  - `manifest(::Source)` — declared URLs, expected sha, license, citation.
  - `fetch!(::Source; cache_dir)` — idempotent HTTP/FTP download with
    resumable transfers (`Downloads.jl` + range requests).
  - `parse(::Source)` — streaming parser producing typed records.
- Storage: write parsed records to **Parquet** (`Parquet2.jl`) partitioned by
  source and date; expose as `DataFrame` via `DataFrames.jl` and as Arrow
  streams via `Arrow.jl` for cross-language use.
- Concurrency: `OhMyThreads.jl` for parallel parses; per-source rate limits
  enforced via a token-bucket in `RDDataSources`.
- Smoke test in CI: download a 1-MB sample of each source and parse it.

Exit: `make data-small` produces a complete small fixture under
`data/processed/small/` in under five minutes.

### Phase 3 — Ontology Layer (Weeks 3–4)

Goal: a single in-memory graph that unifies disease and phenotype concepts.

- `RDOntology.jl`:
  - Parse OBO/OWL with a hand-rolled OBO parser plus `EzXML.jl` for OWL.
  - Build a `MetaGraphsNext.jl` graph (`Graphs.jl` ecosystem) with typed
    nodes (`OntologyTerm`) and edges (`is_a`, `part_of`, `has_phenotype`,
    `xref`).
  - Equivalence resolution across MONDO/Orphanet/OMIM/DOID using
    `oboInOwl:hasDbXref` edges, with conflict logging.
  - Indices: term-by-id, term-by-synonym (case-folded), ancestor/descendant
    bitsets via `BitSet` per term for O(1) "is-a" reachability checks.
- Provide a `phenotype_similarity(a, b; method=:resnik)` function using
  Resnik / Lin / Jiang–Conrath information content, with IC computed from
  disease-phenotype annotation frequencies.

Exit: HPO + MONDO loaded in <5 s; `phenotype_similarity` matches reference
values from published HPO benchmarks within 1e-6.

### Phase 4 — Genomics Layer (Weeks 4–6)

Goal: variant- and gene-level reasoning.

- `RDGenomics.jl`:
  - VCF I/O via `GeneticVariation.jl` / `VariantCallFormat.jl` (or wrap HTSlib
    via `XAM.jl` companion if needed).
  - Coordinate models on GRCh37 and GRCh38; liftover via `Chain` files
    (`Liftover.jl` or our own minimal port).
  - ClinVar ingestion: variant → disease assertions with review status.
  - gnomAD AF lookup keyed by `(chrom, pos, ref, alt, assembly)`.
  - Variant annotation: gene, transcript, HGVS.c/HGVS.p via Ensembl GTF and
    `BioSequences.jl` translation. Don't try to reimplement VEP; instead
    provide a thin Julia API that runs VEP via container if available, else
    falls back to in-Julia annotation for SNVs and small indels in coding
    exons.
  - ACMG-style classification skeleton: PVS1/PM2/PP3 evidence codes wired to
    ClinVar + gnomAD + in silico predictors (CADD, REVEL — distributed as
    pre-computed score files).

Exit: given a ClinVar variant ID, return disease(s), gene, HGVS, AF, and a
provenance trail.

### Phase 5 — Proteomics & Structure Layer (Weeks 6–7)

Goal: link genes/variants to protein consequences and 3-D structure.

- `RDProteomics.jl`:
  - UniProt XML parser → `Protein` records with features (domains, active
    sites, PTMs, variants).
  - AlphaFold structure fetcher (PDB or mmCIF) by `UniProtAcc`.
  - `BioStructures.jl` for parsing; expose residue-level access.
  - Map gene-coordinate variants to UniProt residue positions via Ensembl
    canonical transcript + MANE.
  - Pre-computed per-residue features: pLDDT bin, secondary structure
    (DSSP via `BioStructures` if available), domain membership, distance to
    nearest active site.

Exit: given a missense variant, return UniProt position, AlphaFold pLDDT at
that residue, containing domain, and distance to nearest annotated site.

### Phase 6 — Pathways & Networks Layer (Week 7)

Goal: place genes/proteins inside mechanistic networks.

- `RDPathways.jl`:
  - Reactome BioPAX/SBML parser → typed `Pathway`, `Reaction`, `Complex`.
  - WikiPathways GPML parser.
  - SIGNOR signed graph; merge with STRING for unsigned PPI.
  - BioModels SBML loader (use Phase 9's SBMLToolkit hook).
- Network queries: shortest path, neighborhood, centrality
  (`Graphs.jl` + `GraphsFlows.jl`), with edge weights from confidence scores.

Exit: query "what pathways contain gene X?" returns Reactome IDs with
evidence and one-line descriptions.

### Phase 7 — Pharmacology Layer (Week 8)

Goal: drugs, targets, and orphan-drug context.

- `RDPharmacology.jl`:
  - ChEMBL SQLite snapshot loader; bioactivity queries by target UniProt.
  - DrugBank-open + FDA Orange Book + OOPD for approval/orphan status.
  - RxNorm graph for ingredient → brand → product mapping.
  - DailyMed SPL parser for structured indications/contraindications/dosing.
  - Drug-target affinities exposed as `(drug, target, pKi, source)` rows.

Exit: given a UniProt accession, list drugs with affinity ≤ 1 μM, their
approval status, and any orphan-disease indications.

### Phase 8 — Clinical Layer (Week 9)

Goal: link diseases to trials, registries, and outcomes data.

- `RDClinical.jl`:
  - AACT (ClinicalTrials.gov relational dump) via Postgres or as Parquet.
    Provide a `trials_for(disease)` returning trial records joined to
    interventions and outcomes.
  - FAERS quarterly loader; signal-detection helpers (PRR, ROR) with
    multiple-testing caveats clearly documented.
  - GARD and Orphanet disease summaries with epidemiology
    (prevalence/incidence) attached to `Disease`.

Exit: a disease record exposes `.summary`, `.epidemiology`, `.trials`,
`.adverse_event_signals`, each with provenance.

### Phase 9 — Simulation Engine (Weeks 10–13)

Goal: actually simulate disease mechanisms and progression.

Three layers of model, each appropriate for different questions:

#### 9.1 Molecular & cellular dynamics (SciML core)
- `RDSimulation.jl` built on `DifferentialEquations.jl`,
  `ModelingToolkit.jl`, `Catalyst.jl`, and `SBMLToolkit.jl`.
- Provide an importer that turns a BioModels SBML model into a
  `ReactionSystem`, with a registry mapping pathway IDs to canonical models.
- For diseases where no curated model exists, build a "default-mechanism"
  generator: from a small set of disease-relevant pathways and known
  loss-/gain-of-function variants, instantiate a toy mass-action model with
  rate parameters drawn from BRENDA / SABIO-RK with explicit uncertainty.
- Sensitivity analysis (`GlobalSensitivity.jl`), parameter estimation
  (`SciMLSensitivity.jl` + `Optimization.jl`), and Bayesian inference
  (`Turing.jl`) wired in as first-class workflows.

#### 9.2 Physiological / PBPK & QSP
- Whole-body PBPK templates (Jones-Rowland-style 13-compartment) in
  `ModelingToolkit.jl`, parameterized by species, weight, organ blood flows
  (use open IT'IS tissue parameter sets).
- QSP scaffolds for common disease axes (immune, metabolic, neuro) that
  diseases can plug into.

#### 9.3 Cohort / progression
- Agent-based progression models via `Agents.jl`. Each agent is a
  patient with a genotype draw, environmental covariates, and a per-tick
  state transition driven by either a Markov model (data-fit) or a coupled
  ODE excerpt from §9.1.
- Survival / time-to-event analysis (`Survival.jl`) for calibration against
  registry data when available.

All models emit a `SimulationResult` with: trajectories, parameter
posteriors (if Bayesian), code+data hashes, and a citation list of every
parameter source.

Exit: at least one end-to-end worked example for a tractable disease (e.g.
**phenylketonuria**: PAH enzyme kinetics → phenylalanine accumulation →
cohort outcomes under diet/sapropterin treatment) reproduces published
qualitative behavior.

### Phase 10 — Diagnostic Reasoning (Weeks 13–14)

Goal: phenotype- and variant-driven differential diagnosis.

- `RDDiagnostics.jl`:
  - Input: a `PatientCase` with `Set{HPOId}` observed/excluded phenotypes,
    optional `Vector{Variant}`, demographics, family history.
  - Phenotype-only ranking: implement Phenomizer / Phrank / Exomiser-style
    semantic similarity scoring with information-content weights. Provide
    benchmarks on HPO's published patient cohort.
  - Variant-aware ranking: combine ACMG evidence + gene-phenotype
    likelihood ratios (Exomiser's `hiPHIVE`-like score, reimplemented from
    open papers).
  - Output: ranked `DifferentialDiagnosis` with per-candidate evidence
    breakdown and explicit "not for clinical use" tag.

Exit: on a held-out simulated patient set, top-10 recall ≥ 0.7 (matching
published baselines).

### Phase 11 — Treatment Reasoning (Weeks 14–15)

Goal: candidate-treatment ranking grounded in mechanism and evidence.

- `RDTreatment.jl`:
  - For a given disease, gather: approved drugs (OOPD/Orange Book), drugs
    in trials (AACT), drugs hitting disease-associated targets
    (Pharmacology + Pathways layers), and drug-repurposing candidates
    (network-proximity score against disease module à la Guney 2016 —
    method is openly published).
  - PK/PD layer: hook into PBPK models from 9.2 to estimate steady-state
    target exposure vs IC50.
  - Adverse-event prior from FAERS.
  - Output: `TreatmentCandidate` records with mechanism rationale, evidence
    tier (approved / trial / preclinical / repurposing hypothesis), and the
    same "not for clinical use" tag.

Exit: for the Phase 9 worked example disease, recover known approved
treatments in the top results and produce at least one plausible
repurposing hypothesis with literature support.

### Phase 12 — Interface, Docs, Release (Weeks 15–17)

- `RDApp.jl`: a Genie.jl or Oxygen.jl service exposing:
  - REST: `/disease/{mondo}`, `/variant/{clinvar}`, `/diagnose`, `/simulate`.
  - WebSocket streaming for long simulations.
  - A small HTMX-driven UI; visualizations via `Makie.jl` server-side render
    plus client-side `vega-lite` for interactive charts.
- Pluto notebooks for each worked disease example.
- Documenter sites for each package + a top-level "book" with the data
  dictionary, methods, and limitations.
- Register stable packages (`RareDiseaseCore`, `RDOntology`, `RDGenomics`,
  `RDProteomics`, `RDSimulation`) in the Julia General registry.
- Citable release via Zenodo, DOI in README.

Exit: a tagged `v0.1.0` of each registered package, a deployed docs site,
and a reproducible Docker image that loads the small fixture and runs all
worked examples in under ten minutes.

---

## 4. Cross-Cutting Concerns

### 4.1 Performance
- Prefer `StructArrays.jl` for record collections to keep columns contiguous.
- Memory-map large Arrow/Parquet tables; never load gnomAD or ClinVar fully
  into RAM.
- Use `LoopVectorization.jl` / `Tullio.jl` for hot kernels (e.g. pairwise
  HPO similarity matrices).
- Benchmark each package's hot paths with `BenchmarkTools.jl`; commit
  baseline JSON so regressions are visible in CI.

### 4.2 Testing
- Unit tests per package.
- Integration tests in top-level `test/` that exercise multi-package flows
  (e.g. variant → gene → protein → AlphaFold residue → pathway → drug).
- Property tests via `Supposition.jl` for parsers (round-trip, malformed
  input handling).
- Golden-file tests for ontology graphs (term counts, known ancestor
  relationships).

### 4.3 Documentation & Citation
- Every public function has a docstring with at least one runnable example.
- Every data-derived output includes a `cite()` method that returns BibTeX
  for the sources it consumed.
- Maintain `CITATION.cff` at repo root.

### 4.4 Ethics, Safety, Licensing
- Top-level `ETHICS.md` covering: not-for-clinical-use stance, dual-use
  considerations for variant interpretation, and the absence of patient-level
  data in the repository.
- `LICENSE` MIT for code; `data/LICENSES/` mirrors each source's license text.
- Refuse to ingest sources that prohibit redistribution; downloaders for
  such sources fetch on the user's machine and stay out of the cached
  artifact.

### 4.5 Governance
- `CONTRIBUTING.md` with the package layout rules and the
  "no source without a SOURCES row" gate.
- CODEOWNERS per `packages/*`.
- Quarterly data refresh cycle; a GitHub Action opens a PR with updated
  manifests and a diff summary.

---

## 5. Key Julia Package Choices (Summary Table)

| Concern | Primary | Backup / Notes |
|---|---|---|
| HTTP & download | `Downloads`, `HTTP.jl` | `Curl_jll` |
| Tabular | `DataFrames.jl`, `Arrow.jl`, `Parquet2.jl` | `CSV.jl` |
| Storage | Arrow + Parquet on disk, optional DuckDB via `DuckDB.jl` | SQLite for ChEMBL |
| Graphs | `Graphs.jl`, `MetaGraphsNext.jl` | `GraphIO.jl` |
| Biosequences | `BioSequences.jl`, `FASTX.jl`, `XAM.jl`, `GFF3.jl` | — |
| Variants | `GeneticVariation.jl`, `VariantCallFormat.jl` | HTSlib via FFI |
| Structures | `BioStructures.jl` | — |
| ODE/SDE | `DifferentialEquations.jl`, `ModelingToolkit.jl` | — |
| Reaction networks | `Catalyst.jl`, `SBMLToolkit.jl` | — |
| ABM | `Agents.jl` | — |
| Bayesian | `Turing.jl`, `MCMCChains.jl` | `Pigeons.jl` |
| Optimization | `Optimization.jl`, `Optim.jl` | `JuMP.jl` |
| Sensitivity | `GlobalSensitivity.jl`, `SciMLSensitivity.jl` | — |
| Visualization | `Makie.jl`, `AlgebraOfGraphics.jl` | `Plots.jl` |
| Notebooks | `Pluto.jl` | Jupyter via `IJulia.jl` |
| Web | `Genie.jl` or `Oxygen.jl` | `HTTP.jl` directly |
| Concurrency | `OhMyThreads.jl` | `Dagger.jl` for distributed |
| Interop | `PythonCall.jl`, `RCall.jl`, `CxxWrap.jl` | only where Julia gap exists |

---

## 6. Worked-Example Diseases (Pick Three)

To keep the simulation work honest, commit early to a small set of diseases
with sufficient open mechanism + data:

1. **Phenylketonuria (PKU, MONDO:0009861)** — well-characterized PAH
   enzymology, established dietary and sapropterin therapy, abundant
   pharmacokinetic data. Good first ODE/PBPK example.
2. **Cystic fibrosis (MONDO:0009061)** — CFTR variant landscape, ivacaftor /
   elexacaftor mechanism, registry data via CF Foundation public reports.
   Tests variant-class-specific treatment logic.
3. **Sickle cell disease (MONDO:0011382)** — well-defined single variant,
   hemoglobin polymerization biophysics, hydroxyurea PK/PD, growing
   gene-therapy literature. Tests cohort + intervention modeling.

Each gets its own `notebooks/disease/<slug>.jl` Pluto notebook and a
corresponding integration test.

---

## 7. Milestones & Definition of Done

| Milestone | When | Definition of Done |
|---|---|---|
| M1: Core+Data | end of Phase 2 | small fixture builds, CI green |
| M2: Knowledge graph | end of Phase 6 | disease → genes → proteins → pathways traversal works |
| M3: Diagnosis | end of Phase 10 | differential ranking benchmark hits target |
| M4: Simulation | end of Phase 11 | one disease end-to-end (mechanism → cohort → treatment) |
| M5: v0.1.0 | end of Phase 12 | packages registered, docs deployed, Docker image published |

---

## 8. Risks & Mitigations

- **Source schema drift.** Mitigation: per-source parser tests pinned to a
  committed sample; CI nightly job pulls latest and reports diffs.
- **License creep.** Mitigation: `data/SOURCES.tsv` gate in CI; PR blocked if
  a new download appears without a row.
- **Mechanistic models without parameters.** Mitigation: be explicit about
  parameter uncertainty; default to Bayesian fits with wide priors and
  surface posterior intervals, not point predictions.
- **Scope creep into clinical decision support.** Mitigation: §0 principle 4;
  every diagnostic/treatment output carries the not-for-clinical-use tag and
  the UI banner is non-dismissable.
- **Julia ecosystem gaps (e.g. VEP-equivalent).** Mitigation: define clear
  FFI boundaries early; don't reimplement mature tools.

---

## 8.6 Ecosystem-First Policy

Per maintainer direction: use existing Julia ecosystem packages instead of
reimplementing whenever a maintained option exists. Concrete bindings:

| Concern | Package | Notes |
|---|---|---|
| ODE / SDE integration | `OrdinaryDiffEq` | `Tsit5` default; `DifferentialEquations` umbrella when needed |
| Symbolic models | `ModelingToolkit` | underpins Catalyst |
| Reaction networks / mass-action | `Catalyst` | `@reaction_network` |
| SBML import | `SBML` + `SBMLToolkit` | for BioModels |
| Agent-based models | `Agents` | cohort simulations |
| Graphs (algorithms) | `Graphs` | shortest path, centrality, BFS |
| Typed property graphs | `MetaGraphsNext` | replaces hand-rolled adjacency |
| Bio sequences / FASTA / FASTQ | `FASTX` | |
| Protein structures | `BioStructures` | PDB + mmCIF |
| Bayesian inference | `Turing` | parameter fits |
| Sensitivity analysis | `GlobalSensitivity` | Sobol, Morris |
| Optimization | `Optimization` | unifies Optim/NLopt/etc. |
| HTTP server / routing | `Oxygen` over `HTTP` | FastAPI-style |
| HTTP client / downloads | stdlib `Downloads` | |
| TOML manifests | stdlib `TOML` | |

Hand-rolled retained where no clear ecosystem winner exists:

| Concern | Reason |
|---|---|
| OBO ontology parsing | no maintained pure-Julia OBO parser |
| HGVS nomenclature | no Julia HGVS package (Python `hgvs` via PythonCall is the alternative) |
| PBPK / QSP templates | no Julia library; will build on `ModelingToolkit` + `Catalyst` |
| Guney 2016 network proximity | disease-specific scoring layered on `Graphs.jl` |

## 9. Implementation Progress

Status as of the current commit. Phase numbering follows §3.

| Phase | Component | Status |
|---|---|---|
| 1 | Monorepo + CI + licensing + ETHICS | done |
| 1 | `RareDiseaseCore` identifier types + provenance | done |
| 2 | `RDDataSources` `Source` trait + registry | done |
| 2 | Reference `HPOSource` manifest | done |
| 2 | HTTP downloader + content-addressed cache + sha256 verify + manifest writer | done |
| 2 | Parquet/Arrow processed-data writer | not started |
| 3 | OBO parser (`HPO`/`MONDO`/`Orphanet` shaped) | done |
| 3 | `OntologyGraph` with `ancestors`/`descendants`/`is_a`/`resolve_xref` | done |
| 3 | Information content + Resnik / Lin / Jiang-Conrath similarity | done |
| 4 | HGVS `.c` / `.p` parser | done |
| 4 | ACMG evidence accumulator + Richards 2015 combining rules | done |
| 4 | VCF reader + ClinVar/gnomAD lookups | not started |
| 5 | `RDProteomics` API + BioStructures pass-through | done |
| 5 | UniProt + AlphaFold ingestion | not started |
| 6 | Graphs.jl + MetaGraphsNext-backed `PathwayNetwork` | done |
| 6 | BFS / shortest-path / neighborhood queries | done |
| 6 | Guney 2016 closest-distance + z-scored network proximity | done |
| 6 | Reactome / SIGNOR / WikiPathways loaders | not started |
| 7 | `RDPharmacology` API surface | done |
| 7 | ChEMBL SQLite ingestion + drugs-for-target query | not started |
| 8 | `RDClinical` API surface | done |
| 8 | AACT / FAERS loaders | not started |
| 9 | Simulation types + `RunManifest` | done |
| 9 | Catalyst PAH/PKU reaction network + OrdinaryDiffEq solve | done |
| 9 | Sapropterin one-compartment PBPK | done |
| 9 | CFTR/CF ASL-volume model with variant classes + modulators | done |
| 9 | HbS/SCD polymerization with HbF inhibition | done |
| 9 | `Agents.jl` cohort scaffold | done |
| 9 | SBML/SBMLToolkit BioModels loader | done |
| 10 | Phenotype-only `rank_diagnoses` | done |
| 10 | Variant-aware ranking | not started (needs Phase 4 data) |
| 11 | `TreatmentCandidate` types + evidence tiers | done |
| 11 | Network-proximity-driven `rank_treatments` | done |
| 11 | PK/PD feasibility integration | not started |
| 12 | `RDApp` route table | done |
| 12 | Oxygen.jl route handlers (returning stubs + banner) | done |
| 12 | Domain logic wired into `/diagnose`, `/simulate`, `/treatments` | done |
| 12 | WebSocket streaming for long simulations | not started |
| 12 | Documenter site + Docs.yml workflow | done |
| 12 | Pluto notebooks for PKU + CF + SCD with runnable code | done |

## 10. Immediate Next Actions (First Two Weeks)

1. Land this `PLAN.md` on `main`.
2. Scaffold `packages/RareDiseaseCore.jl` with identifier types and tests.
3. Scaffold `packages/RDDataSources.jl` with the `Source` trait and one
   reference implementation (HPO `.obo`).
4. Set up CI on Julia 1.10 + nightly.
5. Create `data/SOURCES.tsv` with all sources from §2 pre-registered (URLs,
   licenses, citations) even before parsers exist.
6. Open tracking issues for Phases 3–12.
