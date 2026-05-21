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

export
    GenomeAssembly, GRCh37, GRCh38,
    annotate_variant, lookup_clinvar, lookup_gnomad_af,
    acmg_classify, AcmgEvidence

@enum GenomeAssembly GRCh37 GRCh38

"""
    AcmgEvidence

Evidence codes per Richards et al. 2015 (ACMG/AMP).
"""
struct AcmgEvidence
    codes::Vector{Symbol}      # e.g. [:PVS1, :PM2, :PP3]
    rationales::Dict{Symbol,String}
end
AcmgEvidence() = AcmgEvidence(Symbol[], Dict{Symbol,String}())

annotate_variant(::Variant)                  = error("annotate_variant not yet implemented (Phase 4)")
lookup_clinvar(::Variant)                    = error("lookup_clinvar not yet implemented (Phase 4)")
lookup_gnomad_af(::Variant)                  = error("lookup_gnomad_af not yet implemented (Phase 4)")
acmg_classify(::Variant, ::AcmgEvidence)     = error("acmg_classify not yet implemented (Phase 4)")

end # module
