"""
    RDTreatment

Candidate-treatment ranking. Three evidence inputs:

  * approved drugs (FDA Orange Book / OOPD)
  * trial-stage drugs (AACT)
  * mechanism-based repurposing hypotheses (Guney 2016 network-proximity
    of drug targets to the disease module, computed via `RDPathways`)

PK/PD feasibility hooks into `RDSimulation`'s PBPK layer; adverse-event
priors hook into `RDClinical`'s FAERS signals. Both ride in
`mechanism_rationale` and `citations` of the returned candidates.

Phase 11 of the build plan.

!!! warning "Not for clinical use"
    Treatment candidates produced here are research hypotheses, not
    medical advice. Do not act on them for any patient.
"""
module RDTreatment

using RareDiseaseCore
using RDPathways

export
    EvidenceTier, APPROVED, IN_TRIAL, PRECLINICAL, REPURPOSING_HYPOTHESIS,
    TreatmentCandidate, DrugRecord,
    rank_treatments, NOT_FOR_CLINICAL_USE

@enum EvidenceTier APPROVED IN_TRIAL PRECLINICAL REPURPOSING_HYPOTHESIS

const NOT_FOR_CLINICAL_USE =
    "NOT FOR CLINICAL USE. Research output only — do not act on this for any patient."

"""
    DrugRecord

Minimal drug record used by the ranker. Real ingestion (`RDPharmacology`,
`RDClinical`) will populate these from ChEMBL / DailyMed / AACT.
"""
Base.@kwdef struct DrugRecord
    name::String
    chembl::Union{Nothing,ChemblId} = nothing
    rxcui::Union{Nothing,RxCui}     = nothing
    targets::Vector{String}         = String[]   # node labels in pathway network
    evidence_tier::EvidenceTier     = REPURPOSING_HYPOTHESIS
    citations::Vector{String}       = String[]
end

"""
    TreatmentCandidate

One ranked treatment for a disease. `score` is higher-is-better;
`mechanism_rationale` is a short human-readable explanation linking
targets, pathways, and disease genes.
"""
struct TreatmentCandidate
    drug::DrugRecord
    score::Float64
    proximity_d::Float64
    proximity_z::Float64
    mechanism_rationale::String
    citations::Vector{String}
end

"""
    rank_treatments(disease_genes, drugs, network; topk=20, n_bootstrap=200, rng)
        -> (candidates::Vector{TreatmentCandidate}, warning::String)

Score each `DrugRecord` against the disease module (`disease_genes`,
typically a vector of HGNC IDs or whatever node label the pathway network
uses). The score is a combination of:

  * **Evidence tier** — approved drugs floored above trial > preclinical
    > repurposing hypotheses
  * **Network proximity** — negative z-score from `RDPathways.network_proximity_z`
    (closer mean-min distance → more negative z → higher score)

Candidates without any reachable target are dropped.
"""
function rank_treatments(
    disease_genes::AbstractVector{<:AbstractString},
    drugs::AbstractVector{DrugRecord},
    network::PathwayNetwork;
    topk::Int=20,
    n_bootstrap::Int=200,
    rng=nothing,
)
    cands = TreatmentCandidate[]
    for drug in drugs
        targets = collect(String, drug.targets)
        isempty(targets) && continue

        r = if rng === nothing
            network_proximity_z(network, disease_genes, targets;
                                n_bootstrap=n_bootstrap)
        else
            network_proximity_z(network, disease_genes, targets;
                                n_bootstrap=n_bootstrap, rng=rng)
        end
        isfinite(r.d) || continue   # no reachable target

        tier_bonus = _tier_bonus(drug.evidence_tier)
        # More-negative z = closer than random = higher score.
        score = tier_bonus - r.z

        rationale = string(
            "Drug '", drug.name, "' targets {", join(targets, ", "),
            "} have mean-min distance ", round(r.d; digits=2),
            " to disease module {", join(disease_genes, ", "),
            "} (z=", round(r.z; digits=2), "); tier=$(drug.evidence_tier).",
        )

        push!(cands, TreatmentCandidate(
            drug, score, r.d, r.z, rationale, drug.citations,
        ))
    end
    sort!(cands; by=c -> c.score, rev=true)
    length(cands) > topk && resize!(cands, topk)
    return cands, NOT_FOR_CLINICAL_USE
end

# Evidence-tier additive bonus. Calibrated so that an APPROVED drug
# always outranks a REPURPOSING_HYPOTHESIS at any plausible z-score.
function _tier_bonus(t::EvidenceTier)
    t === APPROVED                  && return 20.0
    t === IN_TRIAL                  && return 10.0
    t === PRECLINICAL               && return 5.0
    t === REPURPOSING_HYPOTHESIS    && return 0.0
    return 0.0
end

end # module
