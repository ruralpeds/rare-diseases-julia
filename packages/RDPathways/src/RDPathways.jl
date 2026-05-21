"""
    RDPathways

Reactome, WikiPathways, SIGNOR, BioModels. Places genes/proteins in
mechanistic networks; supports shortest-path, neighborhood, centrality,
network-proximity (Guney 2016) queries, and signed-edge logic.

Phase 6 of the build plan. Built on `Graphs.jl` + `MetaGraphsNext.jl` so
that all general-purpose graph algorithms (Dijkstra, centrality, BFS,
connected components, …) come from the Julia ecosystem; this package
contributes only the disease-specific scoring on top.

# source: Reactome
# source: WikiPathways
# source: SIGNOR
# source: BioModels
# source: STRING
"""
module RDPathways

import Random
using Graphs
using MetaGraphsNext
using RareDiseaseCore

export
    PathwayNetwork, NodeData, EdgeData,
    add_node!, add_edge_undirected!, add_edge_directed!,
    has_node, has_edge,
    neighbors_of, neighborhood, shortest_path_length, shortest_path,
    closest_distance, network_proximity_z,
    load_reactome, load_wikipathways, load_signor,
    pathways_for_gene

"""
    NodeData

Per-vertex metadata stored in the MetaGraph.
"""
mutable struct NodeData
    pathways::Set{String}
end
NodeData() = NodeData(Set{String}())

"""
    EdgeData

Per-edge metadata: `weight` (default 1.0) and `sign` (+1 activation,
−1 inhibition, 0 unsigned).
"""
mutable struct EdgeData
    weight::Float64
    sign::Int
end
EdgeData(; weight::Float64=1.0, sign::Int=0) = EdgeData(weight, sign)

"""
    PathwayNetwork

Thin wrapper around `MetaGraphsNext.MetaGraph` keyed by string labels
(e.g. `"HGNC:8582"`, `"UniProt:P00439"`). The underlying graph is a
`SimpleDiGraph`; undirected edges are simulated by adding both directions.
"""
struct PathwayNetwork
    g::Any  # MetaGraph with closure-typed weight_function; parameterizing
            # buys nothing at the API boundary so we keep this field abstract.
end

function PathwayNetwork()
    g = MetaGraph(
        SimpleDiGraph();
        label_type=String,
        vertex_data_type=NodeData,
        edge_data_type=EdgeData,
        weight_function=ed -> ed.weight,
        default_weight=1.0,
    )
    return PathwayNetwork(g)
end

Base.length(n::PathwayNetwork) = nv(n.g)
has_node(n::PathwayNetwork, id) = haskey(n.g, String(id))
has_edge(n::PathwayNetwork, a, b) = haskey(n.g, String(a), String(b))

function add_node!(n::PathwayNetwork, id::AbstractString;
                   pathways::Vector{<:AbstractString}=String[])
    s = String(id)
    if !haskey(n.g, s)
        n.g[s] = NodeData()
    end
    if !isempty(pathways)
        union!(n.g[s].pathways, String.(pathways))
    end
    return n
end

function add_edge_undirected!(n::PathwayNetwork, a, b;
                              weight::Float64=1.0, sign::Int=0)
    add_node!(n, a); add_node!(n, b)
    sa, sb = String(a), String(b)
    n.g[sa, sb] = EdgeData(weight, sign)
    n.g[sb, sa] = EdgeData(weight, sign)
    return n
end

function add_edge_directed!(n::PathwayNetwork, src, dst;
                            weight::Float64=1.0, sign::Int=0)
    add_node!(n, src); add_node!(n, dst)
    n.g[String(src), String(dst)] = EdgeData(weight, sign)
    return n
end

# ---------------------------------------------------------------------------
# Queries — delegated to Graphs.jl
# ---------------------------------------------------------------------------

"""
    neighbors_of(n, id) -> Vector{String}

Out-neighbors of `id` as labels.
"""
function neighbors_of(n::PathwayNetwork, id)
    s = String(id)
    haskey(n.g, s) || return String[]
    return [label_for(n.g, v) for v in outneighbors(n.g, code_for(n.g, s))]
end

"""
    neighborhood(n, id; radius=1) -> Set{String}

All labels within `radius` hops of `id`, exclusive of `id` itself.
Uses `Graphs.neighborhood`.
"""
function neighborhood(n::PathwayNetwork, id; radius::Int=1)
    s = String(id)
    haskey(n.g, s) || return Set{String}()
    code = code_for(n.g, s)
    codes = Graphs.neighborhood(n.g.graph, code, radius)
    out = Set{String}()
    for c in codes
        c == code && continue
        push!(out, label_for(n.g, c))
    end
    return out
end

"""
    shortest_path_length(n, src, dst) -> Int

Unweighted hop count via `Graphs.gdistances`. `typemax(Int)` if
unreachable; `0` if `src == dst`.
"""
function shortest_path_length(n::PathwayNetwork, src, dst)
    a, b = String(src), String(dst)
    (haskey(n.g, a) && haskey(n.g, b)) || return typemax(Int)
    a == b && return 0
    cs = code_for(n.g, a)
    ct = code_for(n.g, b)
    dists = gdistances(n.g.graph, cs)
    d = dists[ct]
    return d == typemax(eltype(dists)) ? typemax(Int) : Int(d)
end

"""
    shortest_path(n, src, dst) -> Vector{String}

Reconstructed path via `Graphs.a_star`. Empty if unreachable.
"""
function shortest_path(n::PathwayNetwork, src, dst)
    a, b = String(src), String(dst)
    (haskey(n.g, a) && haskey(n.g, b)) || return String[]
    a == b && return [a]
    cs = code_for(n.g, a)
    ct = code_for(n.g, b)
    edges = a_star(n.g.graph, cs, ct)
    isempty(edges) && return String[]
    path = String[label_for(n.g, Graphs.src(first(edges)))]
    for e in edges
        push!(path, label_for(n.g, Graphs.dst(e)))
    end
    return path
end

# ---------------------------------------------------------------------------
# Guney 2016 closest-distance and z-scored proximity
# ---------------------------------------------------------------------------

"""
    closest_distance(n, sources, targets) -> Float64

Mean over `sources` of `min_{t ∈ targets} d(s, t)`. `Inf` if any source
is unreachable from every target, or if either set is empty.
"""
function closest_distance(
    n::PathwayNetwork,
    sources::AbstractVector{<:AbstractString},
    targets::AbstractVector{<:AbstractString},
)
    (isempty(sources) || isempty(targets)) && return Inf
    total = 0.0
    count = 0
    for s in sources
        best = typemax(Int)
        for t in targets
            d = shortest_path_length(n, s, t)
            d < best && (best = d)
        end
        best == typemax(Int) && return Inf
        total += best
        count += 1
    end
    return total / count
end

"""
    network_proximity_z(n, sources, targets; n_bootstrap=200, rng) -> NamedTuple

Z-score of `closest_distance` against random source sets sampled from the
network. Returns `(d, μ, σ, z)`. Degree-stratified sampling can replace
the uniform sampler later without changing the call site.
"""
function network_proximity_z(
    n::PathwayNetwork,
    sources::AbstractVector{<:AbstractString},
    targets::AbstractVector{<:AbstractString};
    n_bootstrap::Int=200,
    rng=Random.default_rng(),
)
    d = closest_distance(n, sources, targets)
    isinf(d) && return (d=d, μ=Inf, σ=0.0, z=Inf)
    all_labels = String[label_for(n.g, v) for v in vertices(n.g)]
    isempty(all_labels) && return (d=d, μ=NaN, σ=NaN, z=NaN)
    ds = Float64[]
    for _ in 1:n_bootstrap
        rs = rand(rng, all_labels, length(sources))
        push!(ds, closest_distance(n, rs, targets))
    end
    finite = filter(isfinite, ds)
    isempty(finite) && return (d=d, μ=NaN, σ=NaN, z=NaN)
    μ = sum(finite) / length(finite)
    σ² = sum((x - μ)^2 for x in finite) / max(length(finite) - 1, 1)
    σ = sqrt(σ²)
    z = σ == 0 ? 0.0 : (d - μ) / σ
    return (d=d, μ=μ, σ=σ, z=z)
end

# Loader stubs — SBML/BioPAX/GPML parsers will land via SBMLToolkit + EzXML.

load_reactome(::AbstractString)        = error("load_reactome not yet implemented")
load_wikipathways(::AbstractString)    = error("load_wikipathways not yet implemented")
load_signor(::AbstractString)          = error("load_signor not yet implemented")

"""
    pathways_for_gene(n, id) -> Set{String}
"""
pathways_for_gene(n::PathwayNetwork, id) =
    haskey(n.g, String(id)) ? n.g[String(id)].pathways : Set{String}()

end # module
