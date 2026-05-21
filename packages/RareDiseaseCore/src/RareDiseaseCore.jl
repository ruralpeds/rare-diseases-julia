"""
    RareDiseaseCore

Core identifier types, base records, and provenance utilities shared across
every package in the rare-diseases-julia monorepo.

Identifiers are validated, namespaced wrappers around strings/integers so that
the type system prevents mixing, e.g., an `HPOId` with a `MondoId`.
"""
module RareDiseaseCore

using Dates
using SHA

export
    # Identifier types
    HPOId, MondoId, OrphaId, OmimId, DoId, MeshId,
    HgncId, EnsemblGeneId, EnsemblTranscriptId, RefSeqId, UniProtAcc,
    ClinVarId, DbSnpId, PubMedId, PmcId, RxCui, ChemblId, PubChemCid,
    ReactomeId, PfamId, InterProId,
    # Provenance
    DataProvenance,
    # Base records
    Disease, Phenotype, Gene, Variant, Protein, Drug, Trial, Pathway,
    # Helpers
    cite, sha256_file

# ---------------------------------------------------------------------------
# Identifier types
# ---------------------------------------------------------------------------

abstract type AbstractId end

Base.show(io::IO, id::AbstractId) = print(io, idstring(id))
Base.string(id::AbstractId) = idstring(id)
Base.:(==)(a::T, b::T) where {T<:AbstractId} = idstring(a) == idstring(b)
Base.hash(id::AbstractId, h::UInt) = hash((typeof(id), idstring(id)), h)

# Each identifier validates its format at construction. We reject empties and
# enforce the documented namespace prefix when applicable.

"""Construct an identifier with prefix `P` and numeric body of width ≥ `W`."""
function _validate_prefixed(s::AbstractString, prefix::AbstractString, minwidth::Int)
    isempty(s) && throw(ArgumentError("empty id"))
    startswith(s, prefix) || throw(ArgumentError("id '$s' must start with '$prefix'"))
    body = SubString(s, lastindex(prefix) + 1)
    all(isdigit, body) || throw(ArgumentError("id '$s' body must be digits"))
    length(body) >= minwidth ||
        throw(ArgumentError("id '$s' body must have at least $minwidth digits"))
    return String(s)
end

# Phenotype / disease ontologies
struct HPOId <: AbstractId
    value::String
    HPOId(s::AbstractString) = new(_validate_prefixed(s, "HP:", 7))
end
struct MondoId <: AbstractId
    value::String
    MondoId(s::AbstractString) = new(_validate_prefixed(s, "MONDO:", 7))
end
struct OrphaId <: AbstractId
    value::String
    OrphaId(s::AbstractString) = new(_validate_prefixed(s, "ORPHA:", 1))
end
struct OmimId <: AbstractId
    value::String
    OmimId(s::AbstractString) = new(_validate_prefixed(s, "OMIM:", 6))
end
struct DoId <: AbstractId
    value::String
    DoId(s::AbstractString) = new(_validate_prefixed(s, "DOID:", 1))
end
struct MeshId <: AbstractId
    value::String
    function MeshId(s::AbstractString)
        isempty(s) && throw(ArgumentError("empty id"))
        # MeSH descriptors are e.g. "D012173", "C535309"
        m = match(r"^[A-Z]\d{6,7}$", s)
        m === nothing && throw(ArgumentError("invalid MeSH id '$s'"))
        new(String(s))
    end
end

# Genes / transcripts
struct HgncId <: AbstractId
    value::String
    HgncId(s::AbstractString) = new(_validate_prefixed(s, "HGNC:", 1))
end
struct EnsemblGeneId <: AbstractId
    value::String
    function EnsemblGeneId(s::AbstractString)
        m = match(r"^ENS[A-Z]{0,3}G\d{11}(\.\d+)?$", s)
        m === nothing && throw(ArgumentError("invalid Ensembl gene id '$s'"))
        new(String(s))
    end
end
struct EnsemblTranscriptId <: AbstractId
    value::String
    function EnsemblTranscriptId(s::AbstractString)
        m = match(r"^ENS[A-Z]{0,3}T\d{11}(\.\d+)?$", s)
        m === nothing && throw(ArgumentError("invalid Ensembl transcript id '$s'"))
        new(String(s))
    end
end
struct RefSeqId <: AbstractId
    value::String
    function RefSeqId(s::AbstractString)
        m = match(r"^[NX][MRP]_\d+(\.\d+)?$", s)
        m === nothing && throw(ArgumentError("invalid RefSeq id '$s'"))
        new(String(s))
    end
end

# Proteins
struct UniProtAcc <: AbstractId
    value::String
    function UniProtAcc(s::AbstractString)
        # Canonical UniProt accession regex from UniProt's docs.
        m = match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", s)
        m === nothing && throw(ArgumentError("invalid UniProt accession '$s'"))
        new(String(s))
    end
end

# Variants & sequence DBs
struct ClinVarId <: AbstractId
    value::String
    function ClinVarId(s::AbstractString)
        all(isdigit, s) || throw(ArgumentError("ClinVar id must be digits, got '$s'"))
        new(String(s))
    end
end
struct DbSnpId <: AbstractId
    value::String
    DbSnpId(s::AbstractString) = new(_validate_prefixed(s, "rs", 1))
end

# Literature
struct PubMedId <: AbstractId
    value::String
    function PubMedId(s::AbstractString)
        all(isdigit, s) || throw(ArgumentError("PMID must be digits, got '$s'"))
        new(String(s))
    end
end
struct PmcId <: AbstractId
    value::String
    PmcId(s::AbstractString) = new(_validate_prefixed(s, "PMC", 1))
end

# Drugs & chemistry
struct RxCui <: AbstractId
    value::String
    function RxCui(s::AbstractString)
        all(isdigit, s) || throw(ArgumentError("RxCUI must be digits, got '$s'"))
        new(String(s))
    end
end
struct ChemblId <: AbstractId
    value::String
    ChemblId(s::AbstractString) = new(_validate_prefixed(s, "CHEMBL", 1))
end
struct PubChemCid <: AbstractId
    value::String
    function PubChemCid(s::AbstractString)
        all(isdigit, s) || throw(ArgumentError("PubChem CID must be digits, got '$s'"))
        new(String(s))
    end
end

# Pathways & domains
struct ReactomeId <: AbstractId
    value::String
    ReactomeId(s::AbstractString) = new(_validate_prefixed(s, "R-", 1))
end
struct PfamId <: AbstractId
    value::String
    PfamId(s::AbstractString) = new(_validate_prefixed(s, "PF", 5))
end
struct InterProId <: AbstractId
    value::String
    InterProId(s::AbstractString) = new(_validate_prefixed(s, "IPR", 6))
end

"""Return the canonical string representation of an identifier."""
idstring(id::AbstractId) = id.value

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

"""
    DataProvenance(source, version, sha256, retrieved_at; url=nothing, citation=nothing)

Provenance metadata attached to every persisted record. `source` should
match a name in `data/SOURCES.tsv`.
"""
struct DataProvenance
    source::String
    version::String
    sha256::String
    retrieved_at::DateTime
    url::Union{Nothing,String}
    citation::Union{Nothing,String}
end
DataProvenance(source, version, sha256, retrieved_at; url=nothing, citation=nothing) =
    DataProvenance(source, version, sha256, retrieved_at, url, citation)

# ---------------------------------------------------------------------------
# Base records
#
# These structs are intentionally light. Domain packages may extend them with
# their own fields (we keep them mutable=false and field-stable so column
# stores stay efficient).
# ---------------------------------------------------------------------------

struct Phenotype
    id::HPOId
    name::String
    definition::Union{Nothing,String}
    provenance::DataProvenance
end

struct Disease
    id::MondoId
    name::String
    synonyms::Vector{String}
    xrefs::Vector{String}              # cross-ref IDs as namespaced strings
    inheritance::Vector{String}        # e.g. ["autosomal_recessive"]
    prevalence::Union{Nothing,Float64} # per 100,000 if known
    provenance::DataProvenance
end

struct Gene
    hgnc::HgncId
    symbol::String
    ensembl::Union{Nothing,EnsemblGeneId}
    chromosome::Union{Nothing,String}
    provenance::DataProvenance
end

struct Variant
    chrom::String
    pos::Int
    ref::String
    alt::String
    assembly::String              # "GRCh37" | "GRCh38"
    clinvar::Union{Nothing,ClinVarId}
    rsid::Union{Nothing,DbSnpId}
    hgvs_c::Union{Nothing,String}
    hgvs_p::Union{Nothing,String}
    provenance::DataProvenance
end

struct Protein
    accession::UniProtAcc
    name::String
    length::Int
    gene_symbol::Union{Nothing,String}
    provenance::DataProvenance
end

struct Drug
    rxcui::Union{Nothing,RxCui}
    chembl::Union{Nothing,ChemblId}
    name::String
    is_approved::Bool
    orphan_designations::Vector{String}  # disease names or MONDO strings
    provenance::DataProvenance
end

struct Trial
    nct_id::String
    phase::Union{Nothing,String}
    status::String
    conditions::Vector{String}
    interventions::Vector{String}
    provenance::DataProvenance
end

struct Pathway
    id::ReactomeId
    name::String
    species::String
    gene_members::Vector{HgncId}
    provenance::DataProvenance
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
    cite(p::DataProvenance) -> String

Return a one-line, human-readable citation derived from provenance.
"""
function cite(p::DataProvenance)
    base = p.citation === nothing ? p.source : p.citation
    return "$base ($(p.source) v$(p.version), retrieved $(Date(p.retrieved_at)))"
end

"""
    sha256_file(path) -> String

Return the lowercase hex SHA-256 of a file on disk. Used by `RDDataSources`
to verify downloads and stamp `DataProvenance`.
"""
function sha256_file(path::AbstractString)
    ctx = SHA.SHA2_256_CTX()
    open(path, "r") do io
        buf = Vector{UInt8}(undef, 1 << 20)
        while !eof(io)
            n = readbytes!(io, buf)
            SHA.update!(ctx, view(buf, 1:n))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

end # module
