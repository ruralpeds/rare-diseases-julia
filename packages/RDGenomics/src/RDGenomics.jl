"""
    RDGenomics

Variant- and gene-level reasoning. Reads VCF, joins ClinVar assertions and
gnomAD allele frequencies, annotates to HGVS, applies ACMG/AMP-style evidence
codes.

Phase 4 of the build plan.

# source: ClinVar
# source: gnomAD
# source: GENCODE
# source: HGNC
# source: MANE
"""
module RDGenomics

using RareDiseaseCore

include("hgvs.jl")

export
    GenomeAssembly, GRCh37, GRCh38,
    HgvsCoding, HgvsProtein,
    parse_hgvs_c, parse_hgvs_p,
    annotate_variant, lookup_clinvar, lookup_gnomad_af,
    acmg_classify, AcmgEvidence

@enum GenomeAssembly GRCh37 GRCh38

"""
    AcmgEvidence

Evidence codes per Richards et al. 2015 (ACMG/AMP). Codes are stored as
`Symbol`s — `:PVS1`, `:PS1`..`:PS4`, `:PM1`..`:PM6`, `:PP1`..`:PP5`,
`:BA1`, `:BS1`..`:BS4`, `:BP1`..`:BP7`.
"""
struct AcmgEvidence
    codes::Vector{Symbol}
    rationales::Dict{Symbol,String}
end
AcmgEvidence() = AcmgEvidence(Symbol[], Dict{Symbol,String}())

function Base.push!(e::AcmgEvidence, code::Symbol, rationale::AbstractString="")
    code in e.codes || push!(e.codes, code)
    isempty(rationale) || (e.rationales[code] = String(rationale))
    return e
end

annotate_variant(::Variant)                  = error("annotate_variant not yet implemented (Phase 4)")
lookup_clinvar(::Variant)                    = error("lookup_clinvar not yet implemented (Phase 4)")
lookup_gnomad_af(::Variant)                  = error("lookup_gnomad_af not yet implemented (Phase 4)")

"""
    acmg_classify(evidence::AcmgEvidence) ->
        Symbol  # :pathogenic | :likely_pathogenic | :uncertain |
                # :likely_benign | :benign

Apply the ACMG/AMP combining rules from Richards 2015 Table 5. Inputs that
satisfy no rule return `:uncertain`.

This is a faithful reading of the published rules; it is **not** a substitute
for expert curation, and the per-code evidence is itself a research output.
"""
function acmg_classify(e::AcmgEvidence)
    has(c) = c in e.codes
    n(prefix) = count(c -> startswith(String(c), prefix), e.codes)

    # Pathogenic rules (any one is sufficient)
    pathogenic =
        (has(:PVS1) && (n("PS") ≥ 1 || n("PM") ≥ 2 || (n("PM") ≥ 1 && n("PP") ≥ 1) || n("PP") ≥ 2)) ||
        (n("PS") ≥ 2) ||
        (n("PS") ≥ 1 && (n("PM") ≥ 3 || (n("PM") ≥ 2 && n("PP") ≥ 2) || (n("PM") ≥ 1 && n("PP") ≥ 4)))

    likely_pathogenic =
        (has(:PVS1) && n("PM") ≥ 1) ||
        (n("PS") ≥ 1 && (1 ≤ n("PM") ≤ 2)) ||
        (n("PS") ≥ 1 && n("PP") ≥ 2) ||
        (n("PM") ≥ 3) ||
        (n("PM") ≥ 2 && n("PP") ≥ 2) ||
        (n("PM") ≥ 1 && n("PP") ≥ 4)

    benign = has(:BA1) || n("BS") ≥ 2
    likely_benign = (n("BS") ≥ 1 && n("BP") ≥ 1) || n("BP") ≥ 2

    pathogenic        && return :pathogenic
    likely_pathogenic && return :likely_pathogenic
    benign            && return :benign
    likely_benign     && return :likely_benign
    return :uncertain
end

end # module
