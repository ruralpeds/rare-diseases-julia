"""
    RDApp

REST + WebSocket service and HTMX-driven UI over the rest of the platform.
Built on `Genie.jl` or `Oxygen.jl` (decision in Phase 12). Every response
that surfaces diagnostic or treatment output carries the not-for-clinical-use
banner non-dismissibly.
"""
module RDApp

using RareDiseaseCore

export start_server, routes

const NOT_FOR_CLINICAL_USE_BANNER =
    "This service is for research and education only. Not for clinical use."

"""
    routes() -> Vector{NamedTuple}

The planned route table. Useful as a single source of truth for docs and
for client SDK generation.
"""
function routes()
    return [
        (method="GET",  path="/health"),
        (method="GET",  path="/disease/:mondo"),
        (method="GET",  path="/variant/:clinvar"),
        (method="GET",  path="/protein/:uniprot"),
        (method="GET",  path="/pathway/:reactome"),
        (method="POST", path="/diagnose"),
        (method="POST", path="/simulate"),
        (method="POST", path="/treatments"),
        (method="GET",  path="/sources"),
    ]
end

start_server(; host="127.0.0.1", port=8000) =
    error("start_server not yet implemented (Phase 12)")

end # module
