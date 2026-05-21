"""
    RDDiagnostics

Phenotype- and variant-driven differential-diagnosis ranking using
Phenomizer / Phrank / Exomiser-style information-content scores and ACMG
evidence integration.

Phase 10 of the build plan. Phenotype-only ranking is implemented now;
variant-aware ranking lands when `RDGenomics` Phase 4 work completes.

!!! warning "Not for clinical use"
    All outputs of this module are research hypotheses. They are not
    validated for clinical decision making and must not be used for
    diagnosis or treatment of any patient.
"""
module RDDiagnostics

using RareDiseaseCore
using RDOntology

export
    PatientCase, DifferentialCandidate, DifferentialDiagnosis,
    rank_diagnoses, NOT_FOR_CLINICAL_USE

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
    disease::String         # ontology id; MONDO/OMIM/ORPHA depending on source
    name::String
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

"""
    rank_diagnoses(case, g; disease_annotations, topk=20) -> DifferentialDiagnosis

Rank candidate diseases by best-pairs Resnik similarity between the
patient's phenotypes and each disease's annotated phenotype set.

* `g` is an `OntologyGraph` with `information_content!` already applied
  using the same `disease_annotations` passed here.
* `disease_annotations` is `disease_id => (name, phenotype_ids)`.
"""
function rank_diagnoses(
    case::PatientCase,
    g::OntologyGraph;
    disease_annotations::AbstractDict,
    topk::Int=20,
)
    query = [string(p) for p in case.phenotypes_present]
    excluded = Set(string(p) for p in case.phenotypes_absent)

    cands = DifferentialCandidate[]
    for (did, payload) in disease_annotations
        name, phs = _unpack_annotation(payload)
        phs_kept = [p for p in phs if !(string(p) in excluded)]
        isempty(phs_kept) && continue

        pscore = isempty(query) ? 0.0 :
            phenotype_similarity(g, query, phs_kept; method=:resnik)

        # Variant scoring is a no-op until Phase 4 supplies real evidence.
        vscore = 0.0

        push!(cands, DifferentialCandidate(
            String(did), name,
            pscore + vscore, pscore, vscore,
            ["phenotype Resnik vs $(length(phs_kept)) annotations"],
        ))
    end

    sort!(cands; by=c -> c.score, rev=true)
    if length(cands) > topk
        resize!(cands, topk)
    end
    return DifferentialDiagnosis(cands, NOT_FOR_CLINICAL_USE)
end

# Accept either `[phenotype_ids...]` or `(name, [phenotype_ids...])`.
function _unpack_annotation(x::Tuple)
    length(x) == 2 || throw(ArgumentError("expected (name, phenotypes) tuple"))
    return (String(x[1]), collect(x[2]))
end
function _unpack_annotation(x::AbstractVector)
    return ("", collect(x))
end
function _unpack_annotation(x)
    return ("", collect(x))
end

end # module
