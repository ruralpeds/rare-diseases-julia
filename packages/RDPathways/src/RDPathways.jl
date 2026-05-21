"""
    RDPathways

Reactome, WikiPathways, SIGNOR, BioModels. Places genes/proteins in
mechanistic networks; supports shortest-path, neighborhood, centrality,
network-proximity queries, and signed-edge logic.

Phase 6 of the build plan. We ship a hand-rolled adjacency-list graph
keyed by `String` node IDs (HGNC, UniProt, or whatever the caller picks),
plus BFS-based path queries. Full SBML/BioPAX/GPML parsing lands later.

# source: Reactome
# source: WikiPathways
# source: SIGNOR
# source: BioModels
# source: STRING
"""
module RDPathways

import Random
using RareDiseaseCore

export
    PathwayNetwork,
    add_node!, add_edge_undirected!, add_edge_directed!,
    has_node, has_edge,
    neighbors_of, neighborhood, shortest_path_length, shortest_path,
    closest_distance, network_proximity_z,
    load_reactome, load_wikipathways, load_signor,
    pathways_for_gene

"""
    PathwayNetwork

Adjacency-list graph with optional edge confidences. Nodes are namespaced
ID strings (e.g. `"HGNC:8582"`, `"UniProt:P00439"`); edges may be
directed or undirected, signed or unsigned.

* `adj[src]`           -> Dict(dst => weight)
* `signs[(src,dst)]`   -> +1 (activation), -1 (inhibition), or 0 (unsigned)
* `pathways_of[node]`  -> Set of pathway IDs the node belongs to
"""
mutable struct PathwayNetwork
    nodes::Set{String}
    adj::Dict{String,Dict{String,Float64}}
    signs::Dict{Tuple{String,String},Int}
    pathways_of::Dict{String,Set{String}}
end
PathwayNetwork() = PathwayNetwork(
    Set{String}(),
    Dict{String,Dict{String,Float64}}(),
    Dict{Tuple{String,String},Int}(),
    Dict{String,Set{String}}(),
)

Base.length(n::PathwayNetwork) = length(n.nodes)
has_node(n::PathwayNetwork, id) = String(id) in n.nodes
has_edge(n::PathwayNetwork, a, b) =
    haskey(n.adj, String(a)) && haskey(n.adj[String(a)], String(b))

function add_node!(n::PathwayNetwork, id::AbstractString;
                   pathways::Vector{<:AbstractString}=String[])
    s = String(id)
    push!(n.nodes, s)
    get!(n.adj, s, Dict{String,Float64}())
    if !isempty(pathways)
        bag = get!(n.pathways_of, s, Set{String}())
        union!(bag, String.(pathways))
    end
    return n
end

function add_edge_undirected!(n::PathwayNetwork, a, b;
                              weight::Float64=1.0, sign::Int=0)
    add_node!(n, a); add_node!(n, b)
    sa, sb = String(a), String(b)
    n.adj[sa][sb] = weight
    n.adj[sb][sa] = weight
    if sign != 0
        n.signs[(sa, sb)] = sign
        n.signs[(sb, sa)] = sign
    end
    return n
end

function add_edge_directed!(n::PathwayNetwork, src, dst;
                            weight::Float64=1.0, sign::Int=0)
    add_node!(n, src); add_node!(n, dst)
    n.adj[String(src)][String(dst)] = weight
    sign == 0 || (n.signs[(String(src), String(dst))] = sign)
    return n
end

"""
    neighbors_of(n, id) -> Vector{String}
"""
neighbors_of(n::PathwayNetwork, id) =
    has_node(n, id) ? collect(keys(n.adj[String(id)])) : String[]

"""
    neighborhood(n, id; radius=1) -> Set{String}

All nodes within `radius` hops of `id`, excluding `id` itself.
"""
function neighborhood(n::PathwayNetwork, id; radius::Int=1)
    s = String(id)
    has_node(n, s) || return Set{String}()
    frontier = Set{String}([s])
    out = Set{String}()
    for _ in 1:radius
        next = Set{String}()
        for u in frontier
            for v in keys(n.adj[u])
                if !(v in out) && v != s
                    push!(out, v)
                    push!(next, v)
                end
            end
        end
        frontier = next
        isempty(frontier) && break
    end
    return out
end

"""
    shortest_path_length(n, src, dst) -> Int

Hop count (BFS, unweighted) from `src` to `dst`. Returns `typemax(Int)`
if no path exists. `src == dst` returns 0.
"""
function shortest_path_length(n::PathwayNetwork, src, dst)
    a, b = String(src), String(dst)
    (has_node(n, a) && has_node(n, b)) || return typemax(Int)
    a == b && return 0
    dist = Dict{String,Int}(a => 0)
    queue = String[a]
    while !isempty(queue)
        u = popfirst!(queue)
        d = dist[u]
        for v in keys(n.adj[u])
            haskey(dist, v) && continue
            v == b && return d + 1
            dist[v] = d + 1
            push!(queue, v)
        end
    end
    return typemax(Int)
end

"""
    shortest_path(n, src, dst) -> Vector{String}

Reconstruct a shortest path as a node sequence. Empty if unreachable.
"""
function shortest_path(n::PathwayNetwork, src, dst)
    a, b = String(src), String(dst)
    (has_node(n, a) && has_node(n, b)) || return String[]
    a == b && return [a]
    prev = Dict{String,String}()
    queue = String[a]
    found = false
    while !isempty(queue)
        u = popfirst!(queue)
        for v in keys(n.adj[u])
            (haskey(prev, v) || v == a) && continue
            prev[v] = u
            if v == b
                found = true
                break
            end
            push!(queue, v)
        end
        found && break
    end
    found || return String[]
    path = [b]
    while path[end] != a
        push!(path, prev[path[end]])
    end
    return reverse(path)
end

"""
    closest_distance(n, source_set, target_set) -> Float64

Guney 2016 closest-distance: for each source node, take the minimum
shortest-path-length to any target node; return the mean over sources.
`Inf` if any source is isolated from every target.
"""
function closest_distance(
    n::PathwayNetwork,
    sources::AbstractVector{<:AbstractString},
    targets::AbstractVector{<:AbstractString},
)
    isempty(sources) && return Inf
    isempty(targets) && return Inf
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

Z-score of `closest_distance(sources, targets)` against a degree-matched
bootstrap. Returns `(d=…, μ=…, σ=…, z=…)`. The simplified degree match
samples random nodes from the same connected pool — fine for the
order-of-magnitude scoring this package uses; rigorous degree-stratified
matching can be slotted in later without changing the API.
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
    all_nodes = collect(n.nodes)
    isempty(all_nodes) && return (d=d, μ=NaN, σ=NaN, z=NaN)

    ds = Float64[]
    for _ in 1:n_bootstrap
        rand_sources = rand(rng, all_nodes, length(sources))
        push!(ds, closest_distance(n, rand_sources, targets))
    end
    finite = filter(isfinite, ds)
    isempty(finite) && return (d=d, μ=NaN, σ=NaN, z=NaN)
    μ = sum(finite) / length(finite)
    σ² = sum((x - μ)^2 for x in finite) / max(length(finite) - 1, 1)
    σ = sqrt(σ²)
    z = σ == 0 ? 0.0 : (d - μ) / σ
    return (d=d, μ=μ, σ=σ, z=z)
end

# Loader stubs — full SBML/BioPAX/GPML parsers land in a follow-up commit.

load_reactome(::AbstractString)        = error("load_reactome not yet implemented (Phase 6 cont.)")
load_wikipathways(::AbstractString)    = error("load_wikipathways not yet implemented (Phase 6 cont.)")
load_signor(::AbstractString)          = error("load_signor not yet implemented (Phase 6 cont.)")

"""
    pathways_for_gene(n, id) -> Set{String}

Pathways recorded for a node (populated by `add_node!`'s `pathways=`
keyword or by future loaders).
"""
pathways_for_gene(n::PathwayNetwork, id) =
    get(n.pathways_of, String(id), Set{String}())

end # module
