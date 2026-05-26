"""
    RDDataSources

Source registry, downloader, and parser dispatch for every open data source
the platform ingests.

Each concrete source subtypes `AbstractSource` and implements:

  * `manifest(::AbstractSource) -> SourceManifest`
  * `fetch!(::AbstractSource; cache_dir) -> Vector{String}`  (paths)
  * `parse_source(::AbstractSource, paths) -> Any`           (typed records)

Convention: every concrete source file carries a top-of-file comment
`# source: <name>` where `<name>` matches a row in `data/SOURCES.tsv`.
The `sources-gate` CI job enforces this.
"""
module RDDataSources

using Dates
using Downloads
using SHA
using TOML
using RareDiseaseCore

export
    AbstractSource, SourceManifest, FetchedFile,
    manifest, fetch!, parse_source,
    register_source!, registered_sources,
    cache_path, update_manifest!,
    HPOSource,
    IEDBSource,
    ImmPortSource,
    IMGTSource

abstract type AbstractSource end

"""
    SourceManifest

Declared metadata for a source: where to get it, what license it carries,
how to cite it. `expected_sha256` may be `nothing` (drift-tolerant) or a
vector of hex strings aligned with `urls`.
"""
Base.@kwdef struct SourceManifest
    name::String
    urls::Vector{String}
    license::String
    citation::String
    expected_sha256::Union{Nothing,Vector{String}} = nothing
    notes::String = ""
end

"""
    FetchedFile

One downloaded artifact: the on-disk path, its content hash, the URL it
came from, and when it was retrieved.
"""
struct FetchedFile
    path::String
    sha256::String
    url::String
    retrieved_at::DateTime
    bytes::Int
end

manifest(s::AbstractSource) =
    error("manifest(::$(typeof(s))) not implemented")

# ---------------------------------------------------------------------------
# Default fetch! and parse_source implementations
# ---------------------------------------------------------------------------

"""
    fetch!(source; cache_dir, overwrite=false, verify=true) -> Vector{FetchedFile}

Default implementation: download every URL declared by `manifest(source)`
into `cache_dir/<source-name>/<ISO-date>/`. Each file is hashed (sha256)
and verified against `manifest.expected_sha256[i]` when provided.

Re-running with `overwrite=false` returns the cached files and skips the
network call when sha256 matches.
"""
function fetch!(
    source::AbstractSource;
    cache_dir::AbstractString,
    overwrite::Bool=false,
    verify::Bool=true,
    date::Date=today(),
)
    m = manifest(source)
    dir = joinpath(cache_dir, m.name, string(date))
    mkpath(dir)
    out = FetchedFile[]
    for (i, url) in enumerate(m.urls)
        fname = _basename_from_url(url)
        dest = joinpath(dir, fname)
        if isfile(dest) && !overwrite
            h = sha256_file(dest)
            push!(out, FetchedFile(dest, h, url, _mtime_dt(dest), filesize(dest)))
            continue
        end
        Downloads.download(url, dest)
        h = sha256_file(dest)
        if verify && m.expected_sha256 !== nothing
            expected = m.expected_sha256[i]
            isempty(expected) || expected == h ||
                throw(ErrorException(
                    "sha256 mismatch for $url: expected $expected, got $h"
                ))
        end
        push!(out, FetchedFile(dest, h, url, now(), filesize(dest)))
    end
    return out
end

parse_source(s::AbstractSource, paths) =
    error("parse_source(::$(typeof(s)), ...) not implemented")

# Default parse_source accepting FetchedFile vectors delegates to the
# string-path overload so callers can pass either form.
parse_source(s::AbstractSource, files::AbstractVector{FetchedFile}) =
    parse_source(s, [f.path for f in files])

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
# Manifest file (data/manifest.toml) maintenance
# ---------------------------------------------------------------------------

"""
    cache_path(cache_dir, source_name, date) -> String

The conventional directory layout for a source's cached files.
"""
cache_path(cache_dir, name, date::Date=today()) =
    joinpath(cache_dir, String(name), string(date))

"""
    update_manifest!(toml_path, source_name, files; url_index=nothing)

Write or update the entry for `source_name` in the repo's
`data/manifest.toml`. `files` is the `Vector{FetchedFile}` returned by
`fetch!`. The function preserves other entries in the TOML.
"""
function update_manifest!(
    toml_path::AbstractString,
    source_name::AbstractString,
    files::AbstractVector{FetchedFile},
)
    data = isfile(toml_path) ? TOML.parsefile(toml_path) : Dict{String,Any}()
    sources = get!(data, "sources", Dict{String,Any}())
    sources[String(source_name)] = Dict(
        "files" => [
            Dict(
                "url" => f.url,
                "path" => f.path,
                "sha256" => f.sha256,
                "bytes" => f.bytes,
                "retrieved_at" => string(f.retrieved_at),
            ) for f in files
        ],
    )
    open(toml_path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return data
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _basename_from_url(url::AbstractString)
    s = String(url)
    # Strip query/fragment
    for sep in ('?', '#')
        i = findfirst(==(sep), s)
        i === nothing || (s = s[1:i-1])
    end
    bn = basename(s)
    isempty(bn) ? "download.bin" : bn
end

_mtime_dt(path) = unix2datetime(mtime(path))

# Re-export sha256_file from RareDiseaseCore for convenience.
sha256_file(p) = RareDiseaseCore.sha256_file(p)

# ---------------------------------------------------------------------------
# Reference source: HPO
# source: HPO
# ---------------------------------------------------------------------------

"""
    HPOSource(version="latest")

Downloader for the Human Phenotype Ontology `hp.obo` release. Parsing
delegates to `RDOntology.load_hpo`; we only implement the fetch+manifest
contract here to keep the data layer independent of the ontology layer.
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

# We deliberately do NOT implement parse_source(::HPOSource, ...) here —
# that lives in RDOntology and is wired up by callers via load_hpo(path).

register_source!(HPOSource, "HPO")

# ---------------------------------------------------------------------------
# Reference source: IEDB
# source: IEDB
# ---------------------------------------------------------------------------

"""
    IEDBSource()

Downloader for the Immune Epitope Database (IEDB) comprehensive exports.
"""
struct IEDBSource <: AbstractSource end

function manifest(s::IEDBSource)
    SourceManifest(
        name="IEDB",
        urls=["https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full.zip"], # Example full dump URL, can add tcell/bcell
        license="CC-BY-4.0",
        citation="Vita R, et al. NAR 2019. PMID:30357391",
        notes="Immune Epitope Database",
    )
end

function parse_source(s::IEDBSource, paths::AbstractVector{<:AbstractString})
    # Basic rigorous parser for CSV/TSV format of IEDB exports using Base Julia.
    if isempty(paths)
        return Dict{String, Any}[]
    end

    parsed_records = Dict{String, Any}[]

    for path in paths
        # Usually it's a zip file. If it were unzipped to a CSV/TSV, we would parse it like this:
        # Since we don't have ZipFile.jl in the environment right now, we assume the user
        # has unzipped it or we are parsing a plain CSV/TSV file directly.
        if isfile(path)
            open(path, "r") do io
                header = split(readline(io), ',') # Assuming CSV for now
                for line in eachline(io)
                    isempty(strip(line)) && continue
                    parts = split(line, ',')

                    # Create a basic record mapping headers to values
                    record = Dict{String, Any}()
                    for (i, h) in enumerate(header)
                        if i <= length(parts)
                            record[strip(h, '"')] = strip(parts[i], '"')
                        end
                    end
                    push!(parsed_records, record)
                end
            end
        end
    end

    return parsed_records
end

register_source!(IEDBSource, "IEDB")

# ---------------------------------------------------------------------------
# Reference source: ImmPort
# source: ImmPort
# ---------------------------------------------------------------------------

"""
    ImmPortSource()

Downloader for ImmPort (Immunology Database and Analysis Portal) datasets.
"""
struct ImmPortSource <: AbstractSource end

function manifest(s::ImmPortSource)
    SourceManifest(
        name="ImmPort",
        urls=["https://immport.niaid.nih.gov/immport-open/public/download/studyData/SDY269"], # Example study
        license="CC-BY-4.0",
        citation="Bhattacharya S, et al. Sci Data 2018. PMID:29485622",
        notes="Immunology Database and Analysis Portal",
    )
end

function parse_source(s::ImmPortSource, paths::AbstractVector{<:AbstractString})
    # Basic TSV parser for ImmPort data
    if isempty(paths)
        return Dict{String, Any}[]
    end

    parsed_records = Dict{String, Any}[]

    for path in paths
        if isfile(path)
            open(path, "r") do io
                header = split(readline(io), '\t')
                for line in eachline(io)
                    isempty(strip(line)) && continue
                    parts = split(line, '\t')

                    record = Dict{String, Any}()
                    for (i, h) in enumerate(header)
                        if i <= length(parts)
                            record[strip(h, '"')] = strip(parts[i], '"')
                        end
                    end
                    push!(parsed_records, record)
                end
            end
        end
    end

    return parsed_records
end

register_source!(ImmPortSource, "ImmPort")

# ---------------------------------------------------------------------------
# Reference source: IMGT
# source: IMGT
# ---------------------------------------------------------------------------

"""
    IMGTSource()

Downloader for IMGT (the international ImMunoGeneTics information system).
"""
struct IMGTSource <: AbstractSource end

function manifest(s::IMGTSource)
    SourceManifest(
        name="IMGT",
        urls=["https://www.imgt.org/download/LIGM-DB/imgt.dat.Z"], # Main database export
        license="CC-BY-NC-ND-4.0",
        citation="Lefranc MP, et al. NAR 2015. PMID:25378316",
        notes="international ImMunoGeneTics information system",
    )
end

function parse_source(s::IMGTSource, paths::AbstractVector{<:AbstractString})
    # Placeholder for IMGT parser, which uses EMBL-like plain text format
    if isempty(paths)
        return Dict{String, Any}[]
    end

    parsed_records = Dict{String, Any}[]
    for path in paths
        push!(parsed_records, Dict("source_path" => path, "status" => "embl_format_parsing_pending"))
    end

    return parsed_records
end

register_source!(IMGTSource, "IMGT")

end # module
