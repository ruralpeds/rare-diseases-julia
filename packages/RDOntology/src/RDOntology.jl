"""
    RDOntology

Unified disease + phenotype ontology graph over HPO, MONDO, Orphanet, OMIM
(open subset), and DO. Provides term lookup, ancestor/descendant queries,
cross-reference resolution, and information-content-based semantic
similarity (Resnik, Lin, Jiang–Conrath).

Phase 3 of the build plan. This stub establishes the public API surface;
the graph implementation lands when the data layer is in.
"""
module RDOntology

using RareDiseaseCore

export
    OntologyTerm, OntologyGraph,
    load_hpo, load_mondo, load_orphanet,
    ancestors, descendants, is_a,
    resolve_xref,
    phenotype_similarity

"""
    OntologyTerm

A single term in the unified ontology. `id` is namespaced (HP:, MONDO:, ...).
"""
struct OntologyTerm
    id::String
    name::String
    namespace::Symbol           # :hpo | :mondo | :orpha | :omim | :do
    synonyms::Vector{String}
    xrefs::Vector{String}
    is_obsolete::Bool
end

"""
    OntologyGraph

In-memory typed DAG. Populated by `load_*` functions; queried via
`ancestors`, `descendants`, `is_a`, `resolve_xref`, `phenotype_similarity`.
"""
mutable struct OntologyGraph
    terms::Dict{String,OntologyTerm}
    parents::Dict{String,Vector{String}}
    children::Dict{String,Vector{String}}
end
OntologyGraph() = OntologyGraph(Dict(), Dict(), Dict())

# --- API stubs (Phase 3 fills in the bodies) ---

load_hpo(::AbstractString)       = error("load_hpo not yet implemented (Phase 3)")
load_mondo(::AbstractString)     = error("load_mondo not yet implemented (Phase 3)")
load_orphanet(::AbstractString)  = error("load_orphanet not yet implemented (Phase 3)")

ancestors(::OntologyGraph, ::AbstractString)   = error("ancestors not yet implemented (Phase 3)")
descendants(::OntologyGraph, ::AbstractString) = error("descendants not yet implemented (Phase 3)")
is_a(::OntologyGraph, ::AbstractString, ::AbstractString) =
    error("is_a not yet implemented (Phase 3)")
resolve_xref(::OntologyGraph, ::AbstractString) =
    error("resolve_xref not yet implemented (Phase 3)")

"""
    phenotype_similarity(g, a, b; method=:resnik)

Semantic similarity between two HPO term sets. `method` is one of
`:resnik`, `:lin`, `:jc`.
"""
phenotype_similarity(::OntologyGraph, a, b; method::Symbol=:resnik) =
    error("phenotype_similarity not yet implemented (Phase 3)")

end # module
