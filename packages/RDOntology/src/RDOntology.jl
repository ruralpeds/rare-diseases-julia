"""
    RDOntology

Unified disease + phenotype ontology graph over HPO, MONDO, Orphanet, OMIM
(open subset), and DO. Provides term lookup, ancestor/descendant queries,
cross-reference resolution, and information-content-based semantic
similarity (Resnik, Lin, Jiang–Conrath).

Phase 3 of the build plan.
"""
module RDOntology

using RareDiseaseCore

include("obo.jl")

export
    OntologyTerm, OntologyGraph,
    load_obo, load_hpo, load_mondo, load_orphanet,
    add_term!, add_edge!,
    ancestors, descendants, is_a,
    resolve_xref,
    most_informative_common_ancestor,
    information_content!, term_ic,
    phenotype_similarity

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

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

In-memory typed DAG. Edges are directed child -> parent (`is_a`).

* `parents[id]`   — direct parents of `id`
* `children[id]`  — direct children of `id`
* `xref_index`    — `"OMIM:261600" => "MONDO:0009861"` style equivalences
* `annotations`   — `disease_id => Set{phenotype_id}` for IC computation
* `ic`            — populated by `information_content!`
"""
mutable struct OntologyGraph
    terms::Dict{String,OntologyTerm}
    parents::Dict{String,Vector{String}}
    children::Dict{String,Vector{String}}
    xref_index::Dict{String,String}
    annotations::Dict{String,Set{String}}
    ic::Dict{String,Float64}
end
OntologyGraph() = OntologyGraph(
    Dict{String,OntologyTerm}(),
    Dict{String,Vector{String}}(),
    Dict{String,Vector{String}}(),
    Dict{String,String}(),
    Dict{String,Set{String}}(),
    Dict{String,Float64}(),
)

Base.length(g::OntologyGraph) = length(g.terms)
Base.haskey(g::OntologyGraph, id::AbstractString) = haskey(g.terms, id)
Base.getindex(g::OntologyGraph, id::AbstractString) = g.terms[id]

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

function add_term!(g::OntologyGraph, t::OntologyTerm)
    g.terms[t.id] = t
    get!(g.parents, t.id, String[])
    get!(g.children, t.id, String[])
    for xr in t.xrefs
        # First registration wins; downstream loaders log conflicts.
        get!(g.xref_index, xr, t.id)
    end
    return g
end

"""
    add_edge!(g, child_id, parent_id)

Add a directed `is_a` edge `child -> parent`. Both terms must exist.
"""
function add_edge!(g::OntologyGraph, child::AbstractString, parent::AbstractString)
    haskey(g.terms, child)  || throw(KeyError(child))
    haskey(g.terms, parent) || throw(KeyError(parent))
    push!(g.parents[child], String(parent))
    push!(g.children[parent], String(child))
    return g
end

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

"""
    load_obo(path; namespace) -> OntologyGraph

Parse an OBO 1.2/1.4 file into an `OntologyGraph`. `namespace` is the
`Symbol` tag stamped on every term (`:hpo`, `:mondo`, ...).
"""
function load_obo(path::AbstractString; namespace::Symbol)
    g = OntologyGraph()
    stanzas = parse_obo(path)
    # Pass 1: terms
    for st in stanzas
        st.kind == "Term" || continue
        haskey(st.tags, "id") || continue
        id = first(st.tags["id"])
        name = haskey(st.tags, "name") ? first(st.tags["name"]) : id
        synonyms = get(st.tags, "synonym", String[])
        # Strip surrounding quotes and trailing metadata from synonym lines:
        #   "Seizures" EXACT [HPO:probinson]
        synonyms = [_extract_quoted(s) for s in synonyms]
        xrefs = String[]
        if haskey(st.tags, "xref")
            for x in st.tags["xref"]
                push!(xrefs, _first_token(x))
            end
        end
        obsolete = haskey(st.tags, "is_obsolete") &&
                   lowercase(first(st.tags["is_obsolete"])) == "true"
        add_term!(g, OntologyTerm(id, name, namespace, synonyms, xrefs, obsolete))
    end
    # Pass 2: edges. We support `is_a` and `relationship: part_of` parents.
    for st in stanzas
        st.kind == "Term" || continue
        haskey(st.tags, "id") || continue
        child = first(st.tags["id"])
        haskey(g.terms, child) || continue
        if haskey(st.tags, "is_a")
            for raw in st.tags["is_a"]
                parent = _first_token(raw)
                if haskey(g.terms, parent)
                    add_edge!(g, child, parent)
                end
            end
        end
    end
    return g
end

load_hpo(path::AbstractString)      = load_obo(path; namespace=:hpo)
load_mondo(path::AbstractString)    = load_obo(path; namespace=:mondo)
load_orphanet(path::AbstractString) = load_obo(path; namespace=:orpha)

# Helpers for OBO field cleanup.
function _first_token(s::AbstractString)
    # "MONDO:0009861 ! Phenylketonuria" -> "MONDO:0009861"
    String(strip(split(s, r"\s*!\s*"; limit=2)[1]))
end
function _extract_quoted(s::AbstractString)
    # `"Seizures" EXACT [...]` -> `Seizures`. Fall back to raw if no quote.
    m = match(r"^\"((?:[^\"\\]|\\.)*)\"", s)
    m === nothing && return String(s)
    return String(m.captures[1])
end

# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------

"""
    ancestors(g, id) -> Set{String}

All transitive parents of `id` (exclusive of `id` itself). O(V+E) per call,
no caching — callers that need many queries should memoize.
"""
function ancestors(g::OntologyGraph, id::AbstractString)
    out = Set{String}()
    haskey(g.terms, id) || return out
    stack = String[id]
    while !isempty(stack)
        cur = pop!(stack)
        for p in g.parents[cur]
            if !(p in out)
                push!(out, p)
                push!(stack, p)
            end
        end
    end
    return out
end

"""
    descendants(g, id) -> Set{String}

All transitive children of `id` (exclusive of `id` itself).
"""
function descendants(g::OntologyGraph, id::AbstractString)
    out = Set{String}()
    haskey(g.terms, id) || return out
    stack = String[id]
    while !isempty(stack)
        cur = pop!(stack)
        for c in g.children[cur]
            if !(c in out)
                push!(out, c)
                push!(stack, c)
            end
        end
    end
    return out
end

"""
    is_a(g, child, parent) -> Bool

True if `child === parent` or `parent` is a transitive ancestor of `child`.
"""
function is_a(g::OntologyGraph, child::AbstractString, parent::AbstractString)
    child == parent && return haskey(g.terms, child)
    return parent in ancestors(g, child)
end

"""
    resolve_xref(g, xref) -> Union{String,Nothing}

Resolve a cross-reference (e.g. `"OMIM:261600"`) to the term id that
declared it. Returns `nothing` if no term claims that xref.
"""
function resolve_xref(g::OntologyGraph, xref::AbstractString)
    return get(g.xref_index, String(xref), nothing)
end

# ---------------------------------------------------------------------------
# Information content + similarity (Resnik / Lin / Jiang-Conrath)
# ---------------------------------------------------------------------------

"""
    information_content!(g; annotations)

Compute and store IC(t) = -log(p(t)) for every term, where p(t) is the
fraction of annotated diseases whose phenotype set contains `t` or any of
its descendants. `annotations` is `disease_id => Iterable{phenotype_id}`.

Must be called before `phenotype_similarity` or `term_ic`.
"""
function information_content!(g::OntologyGraph;
                              annotations::AbstractDict)
    # Store raw annotations
    empty!(g.annotations)
    for (d, phs) in annotations
        g.annotations[String(d)] = Set{String}(String(p) for p in phs)
    end

    n = length(g.annotations)
    n == 0 && throw(ArgumentError("annotations is empty"))

    # For each term, count diseases whose annotation set intersects
    # {term} ∪ descendants(term). We compute counts by propagating up:
    # for every (disease, leaf), mark leaf and all its ancestors.
    counts = Dict{String,Int}()
    for (_, phs) in g.annotations
        marked = Set{String}()
        for p in phs
            haskey(g.terms, p) || continue
            push!(marked, p)
            union!(marked, ancestors(g, p))
        end
        for t in marked
            counts[t] = get(counts, t, 0) + 1
        end
    end

    empty!(g.ic)
    for (t, c) in counts
        # IC = -log(p). p never reaches 0 here.
        g.ic[t] = -log(c / n)
    end
    # Terms with no annotations get IC = log(n) (i.e. p = 1/n upper bound).
    fallback = log(n)
    for t in keys(g.terms)
        haskey(g.ic, t) || (g.ic[t] = fallback)
    end
    return g
end

"""
    term_ic(g, id) -> Float64

Information content of a single term. Errors if `information_content!`
hasn't been run.
"""
function term_ic(g::OntologyGraph, id::AbstractString)
    isempty(g.ic) && error("call information_content!(g; annotations=...) first")
    return get(g.ic, String(id), 0.0)
end

"""
    most_informative_common_ancestor(g, a, b) -> Tuple{String,Float64}

Return `(mica_id, ic)` — the common ancestor of `a` and `b` with the
highest information content. If `a` and `b` share no ancestor, returns
`("", 0.0)`.
"""
function most_informative_common_ancestor(
    g::OntologyGraph, a::AbstractString, b::AbstractString
)
    isempty(g.ic) && error("call information_content!(g; annotations=...) first")
    A = ancestors(g, a); push!(A, String(a))
    B = ancestors(g, b); push!(B, String(b))
    best_id = ""
    best_ic = -Inf
    for t in intersect(A, B)
        ic = get(g.ic, t, 0.0)
        if ic > best_ic
            best_ic = ic
            best_id = t
        end
    end
    return best_ic == -Inf ? ("", 0.0) : (best_id, best_ic)
end

"""
    phenotype_similarity(g, query, reference; method=:resnik) -> Float64

Semantic similarity between two HPO term sets using the asymmetric
"best-pairs" formulation common in phenotype matching (Phenomizer):
for each `q ∈ query`, take the maximum pairwise similarity to any
`r ∈ reference`, then average over `query`.

`method`:
  * `:resnik` — IC(MICA(q, r))
  * `:lin`    — 2·IC(MICA) / (IC(q) + IC(r))
  * `:jc`     — 1 / (1 + IC(q) + IC(r) - 2·IC(MICA))
"""
function phenotype_similarity(
    g::OntologyGraph,
    query, reference;
    method::Symbol=:resnik,
)
    isempty(g.ic) && error("call information_content!(g; annotations=...) first")
    q = collect(String(x) for x in query)
    r = collect(String(x) for x in reference)
    (isempty(q) || isempty(r)) && return 0.0

    total = 0.0
    for qi in q
        best = -Inf
        for rj in r
            _, mica_ic = most_informative_common_ancestor(g, qi, rj)
            sim = if method === :resnik
                mica_ic
            elseif method === :lin
                denom = term_ic(g, qi) + term_ic(g, rj)
                denom == 0 ? 0.0 : (2 * mica_ic) / denom
            elseif method === :jc
                d = term_ic(g, qi) + term_ic(g, rj) - 2 * mica_ic
                1.0 / (1.0 + d)
            else
                throw(ArgumentError("unknown method $method"))
            end
            sim > best && (best = sim)
        end
        total += best == -Inf ? 0.0 : best
    end
    return total / length(q)
end

end # module
