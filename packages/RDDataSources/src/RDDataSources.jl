"""
    RDDataSources

Source registry, downloader, and parser dispatch for every open data source
the platform ingests.

Each concrete source subtypes `AbstractSource` and implements:

  * `manifest(::AbstractSource) -> SourceManifest`
  * `fetch!(::AbstractSource; cache_dir) -> Vector{String}`  (paths)
  * `parse(::AbstractSource, paths) -> Any`                  (typed records)

Convention: every concrete source file carries a top-of-file comment
`# source: <name>` where `<name>` matches a row in `data/SOURCES.tsv`.
The `sources-gate` CI job enforces this.
"""
module RDDataSources

using Dates
using RareDiseaseCore

export
    AbstractSource, SourceManifest,
    manifest, fetch!, parse_source,
    register_source!, registered_sources,
    HPOSource

abstract type AbstractSource end

"""
    SourceManifest

Declared metadata for a source: where to get it, what license it carries,
how to cite it.
"""
Base.@kwdef struct SourceManifest
    name::String
    urls::Vector{String}
    license::String
    citation::String
    expected_sha256::Union{Nothing,Vector{String}} = nothing
    notes::String = ""
end

manifest(s::AbstractSource) =
    error("manifest(::$(typeof(s))) not implemented")
fetch!(s::AbstractSource; cache_dir::AbstractString) =
    error("fetch!(::$(typeof(s))) not implemented")
parse_source(s::AbstractSource, paths) =
    error("parse_source(::$(typeof(s)), ...) not implemented")

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

const _REGISTRY = Dict{String,Type{<:AbstractSource}}()

function register_source!(::Type{T}, name::AbstractString) where {T<:AbstractSource}
    _REGISTRY[String(name)] = T
    return T
end

registered_sources() = sort!(collect(keys(_REGISTRY)))

# ---------------------------------------------------------------------------
# Reference source: HPO
# source: HPO
# ---------------------------------------------------------------------------

"""
    HPOSource(version="latest")

Downloader/parser for the Human Phenotype Ontology (OBO release + annotations).
The actual fetching/parsing implementation lands in Phase 2; this stub exposes
the `manifest()` contract so the registry and CI gate work today.
"""
struct HPOSource <: AbstractSource
    version::String
end
HPOSource() = HPOSource("latest")

function manifest(s::HPOSource)
    v = s.version
    base = "https://github.com/obophenotype/human-phenotype-ontology/releases"
    obo_url = v == "latest" ?
        "$base/latest/download/hp.obo" :
        "$base/download/v$(v)/hp.obo"
    SourceManifest(
        name="HPO",
        urls=[obo_url],
        license="CC-BY-4.0",
        citation="Köhler S, et al. NAR 2021. PMID:33264411",
        notes="Human Phenotype Ontology",
    )
end

register_source!(HPOSource, "HPO")

end # module
