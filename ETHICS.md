# Ethics & Scope

## Not for clinical use

This project is research and educational software. It does **not** provide
medical advice, diagnosis, or treatment. Any output describing candidate
diagnoses or therapies is a hypothesis derived from open data and published
methods, with known and substantial false-positive and false-negative rates.

All user-facing surfaces (REST responses, web UI, notebooks, function
docstrings on `RDDiagnostics` and `RDTreatment`) must carry a
non-dismissable banner stating this.

## Data scope

- We ingest only **open** and **government-released** data.
- We do **not** host patient-level data of any kind.
- We do **not** redistribute data whose license forbids it. Such sources may
  be supported via on-machine downloaders, but their bytes never appear in
  the published artifact or cached image.

## Variant interpretation — dual use

Variant pathogenicity prediction has obvious benefits but can be misused
(e.g. eugenic screening, discrimination by insurers). To reduce harm:

- All variant calls return calibration / confidence intervals, not point
  classifications.
- We do not provide pre-implantation or pre-natal screening workflows.
- We follow the ACMG / AMP guideline structure (evidence codes) rather than
  emitting opaque scores.

## Drug repurposing — responsibility

Network-proximity and target-based repurposing produces many hypotheses.
We surface them with explicit evidence tiers and link every hypothesis to
its supporting publications. Hypotheses without literature support are
labeled as such.

## Reporting concerns

Open an issue with the label `ethics` or email the maintainers privately
(see `CITATION.cff`).
