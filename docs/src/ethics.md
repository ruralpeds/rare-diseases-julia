# Ethics & Scope

See [`ETHICS.md`](https://github.com/ruralpeds/rare-diseases-julia/blob/main/ETHICS.md)
in the repository root for the canonical statement.

Key principles:

1. **Not for clinical use.** No diagnostic or treatment output from
   this platform is validated for clinical decision-making.
2. **Open data only.** No closed sources, no patient-level data.
3. **Variant interpretation is dual-use.** We surface calibration /
   confidence intervals, not opaque pathogenicity scores; we do not
   provide pre-natal or pre-implantation screening workflows.
4. **Repurposing hypotheses are clearly labeled as such.**

Every `RDDiagnostics` and `RDTreatment` API surface carries the
not-for-clinical-use banner constant in its output.
