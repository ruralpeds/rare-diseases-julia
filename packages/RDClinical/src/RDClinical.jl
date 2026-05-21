"""
    RDClinical

ClinicalTrials.gov (AACT relational dump), FDA FAERS, GARD, Orphanet.
Links diseases to active and historical trials, adverse-event signals,
and disease summaries with epidemiology.

Phase 8 of the build plan.

# source: ClinicalTrials-AACT
# source: FDA-FAERS
# source: GARD
# source: Orphanet-ORDO
"""
module RDClinical

using RareDiseaseCore

export
    AdverseEventSignal,
    trials_for, adverse_event_signals, disease_summary

"""
    AdverseEventSignal

Disproportionality-analysis signal from FAERS. `prr` = proportional reporting
ratio; `ror` = reporting odds ratio. Both are observational signals only.
"""
struct AdverseEventSignal
    drug::String
    event::String
    n_reports::Int
    prr::Float64
    ror::Float64
end

trials_for(::MondoId)              = error("trials_for not yet implemented (Phase 8)")
adverse_event_signals(::String)    = error("adverse_event_signals not yet implemented (Phase 8)")
disease_summary(::MondoId)         = error("disease_summary not yet implemented (Phase 8)")

end # module
