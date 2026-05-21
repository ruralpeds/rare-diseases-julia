"""
    RDDiagnostics

Phenotype- and variant-driven differential-diagnosis ranking using
Phenomizer / Phrank / Exomiser-style information-content scores and ACMG
evidence integration.

Phase 10 of the build plan.

!!! warning "Not for clinical use"
    All outputs of this module are research hypotheses. They are not
    validated for clinical decision making and must not be used for
    diagnosis or treatment of any patient.
"""
module RDDiagnostics

using RareDiseaseCore

export
    PatientCase, DifferentialCandidate, DifferentialDiagnosis,
    rank_diagnoses

"""
    PatientCase

Input to the ranker. `phenotypes_present` and `phenotypes_absent` are HPO
IDs; `variants` is an optional vector. `demographics` and `family_history`
are free-form for now.
"""
Base.@kwdef struct PatientCase
    phenotypes_present::Vector{HPOId} = HPOId[]
    phenotypes_absent::Vector{HPOId}  = HPOId[]
    variants::Vector{Variant}         = Variant[]
    demographics::Dict{Symbol,Any}    = Dict{Symbol,Any}()
    family_history::Vector{String}    = String[]
end

"""
    DifferentialCandidate

One ranked candidate disease with its score breakdown.
"""
struct DifferentialCandidate
    disease::MondoId
    score::Float64
    phenotype_score::Float64
    variant_score::Float64
    evidence::Vector{String}
end

struct DifferentialDiagnosis
    candidates::Vector{DifferentialCandidate}   # sorted desc by score
    warning::String
end

const NOT_FOR_CLINICAL_USE =
    "NOT FOR CLINICAL USE. Research output only — do not act on this for any patient."

rank_diagnoses(::PatientCase; topk::Int=20) =
    error("rank_diagnoses not yet implemented (Phase 10)")

end # module
