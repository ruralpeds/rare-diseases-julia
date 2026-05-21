"""
    RDTreatment

Candidate-treatment ranking: approved drugs, trial-stage drugs, and
repurposing hypotheses (network-proximity to a disease module, à la
Guney 2016). PK/PD feasibility hooks into `RDSimulation`'s PBPK layer;
adverse-event priors hook into `RDClinical`'s FAERS signals.

Phase 11 of the build plan.

!!! warning "Not for clinical use"
    Treatment candidates produced here are research hypotheses, not
    medical advice. Do not act on them for any patient.
"""
module RDTreatment

using RareDiseaseCore

export
    EvidenceTier, APPROVED, IN_TRIAL, PRECLINICAL, REPURPOSING_HYPOTHESIS,
    TreatmentCandidate,
    rank_treatments

@enum EvidenceTier APPROVED IN_TRIAL PRECLINICAL REPURPOSING_HYPOTHESIS

"""
    TreatmentCandidate

One ranked treatment for a disease. `evidence_tier` is the strongest tier
of evidence supporting the candidate; `mechanism_rationale` is a short
human-readable explanation linking targets, pathways, and disease genes.
"""
struct TreatmentCandidate
    drug_name::String
    chembl::Union{Nothing,ChemblId}
    rxcui::Union{Nothing,RxCui}
    evidence_tier::EvidenceTier
    targets::Vector{UniProtAcc}
    score::Float64
    mechanism_rationale::String
    citations::Vector{String}
end

const NOT_FOR_CLINICAL_USE =
    "NOT FOR CLINICAL USE. Research output only — do not act on this for any patient."

rank_treatments(::MondoId; topk::Int=20) =
    error("rank_treatments not yet implemented (Phase 11)")

end # module
