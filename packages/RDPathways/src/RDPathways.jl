"""
    RDPathways

Reactome, WikiPathways, SIGNOR, BioModels. Places genes/proteins in
mechanistic networks; supports shortest-path, neighborhood, centrality,
and signed-edge queries.

Phase 6 of the build plan.

# source: Reactome
# source: WikiPathways
# source: SIGNOR
# source: BioModels
# source: STRING
"""
module RDPathways

using RareDiseaseCore

export
    PathwayNetwork,
    load_reactome, load_wikipathways, load_signor,
    pathways_for_gene, shortest_path, neighborhood

mutable struct PathwayNetwork
    pathways::Dict{String,Pathway}
    # adjacency lookups populated by loaders
end
PathwayNetwork() = PathwayNetwork(Dict())

load_reactome(::AbstractString)        = error("load_reactome not yet implemented (Phase 6)")
load_wikipathways(::AbstractString)    = error("load_wikipathways not yet implemented (Phase 6)")
load_signor(::AbstractString)          = error("load_signor not yet implemented (Phase 6)")
pathways_for_gene(::PathwayNetwork, ::HgncId) =
    error("pathways_for_gene not yet implemented (Phase 6)")
shortest_path(::PathwayNetwork, ::HgncId, ::HgncId) =
    error("shortest_path not yet implemented (Phase 6)")
neighborhood(::PathwayNetwork, ::HgncId; radius::Int=1) =
    error("neighborhood not yet implemented (Phase 6)")

end # module
