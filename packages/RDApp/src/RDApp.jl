"""
    RDApp

REST + WebSocket service over the rest of the platform, built on
`Oxygen.jl` (FastAPI-style routing on top of `HTTP.jl`). Every response
that surfaces diagnostic or treatment output carries the
not-for-clinical-use banner non-dismissibly.
"""
module RDApp

using HTTP
using JSON3
using Oxygen
using RareDiseaseCore

export
    routes, register_routes!, start_server,
    NOT_FOR_CLINICAL_USE_BANNER

const NOT_FOR_CLINICAL_USE_BANNER =
    "This service is for research and education only. Not for clinical use."

"""
    routes() -> Vector{NamedTuple}

Documented route table; serves as a single source of truth for docs and
for client SDK generation.
"""
function routes()
    return [
        (method="GET",  path="/health",
         desc="Liveness probe; returns {status:'ok'}."),
        (method="GET",  path="/sources",
         desc="Registered open data sources from SOURCES.tsv."),
        (method="GET",  path="/disease/{mondo}",
         desc="Disease record by MONDO id."),
        (method="GET",  path="/variant/{clinvar}",
         desc="Variant record by ClinVar id."),
        (method="GET",  path="/protein/{uniprot}",
         desc="Protein record by UniProt accession."),
        (method="GET",  path="/pathway/{reactome}",
         desc="Pathway record by Reactome id."),
        (method="POST", path="/diagnose",
         desc="Phenotype + variant differential diagnosis."),
        (method="POST", path="/simulate",
         desc="Run a registered disease model."),
        (method="POST", path="/treatments",
         desc="Rank candidate treatments for a disease."),
    ]
end

"""
    register_routes!(; oxygen=Oxygen) -> nothing

Register every route in `routes()` against the running Oxygen app. Each
handler returns a JSON payload prefixed with the not-for-clinical-use
banner. Domain logic is delegated to the relevant `RD*` package and is
mostly placeholder until the rest of the platform fills in.
"""
function register_routes!(; oxygen::Module=Oxygen)
    oxygen.@get "/health" function (req)
        return Dict("status" => "ok", "banner" => NOT_FOR_CLINICAL_USE_BANNER)
    end
    oxygen.@get "/sources" function (req)
        return Dict(
            "banner"  => NOT_FOR_CLINICAL_USE_BANNER,
            "sources" => _load_sources_tsv(),
        )
    end
    oxygen.@get "/disease/{mondo}" function (req, mondo::String)
        return Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
                    "mondo"=>mondo,
                    "status"=>"not_implemented")
    end
    oxygen.@post "/diagnose" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)))
        return Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
                    "echo"=>body,
                    "status"=>"not_implemented")
    end
    oxygen.@post "/simulate" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)))
        return Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
                    "echo"=>body,
                    "status"=>"not_implemented")
    end
    oxygen.@post "/treatments" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)))
        return Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
                    "echo"=>body,
                    "status"=>"not_implemented")
    end
    return nothing
end

"""
    start_server(; host="127.0.0.1", port=8000)

Register routes and start Oxygen's HTTP server.
"""
function start_server(; host::String="127.0.0.1", port::Int=8000)
    register_routes!()
    serve(host=host, port=port)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _load_sources_tsv(path::AbstractString=joinpath(
        @__DIR__, "..", "..", "..", "data", "SOURCES.tsv"))
    isfile(path) || return []
    out = Vector{Dict{String,String}}()
    open(path, "r") do io
        header = split(strip(readline(io)), '\t')
        for line in eachline(io)
            isempty(strip(line)) && continue
            row = split(line, '\t')
            length(row) < length(header) && continue
            push!(out, Dict(string(header[i]) => string(row[i])
                            for i in eachindex(header)))
        end
    end
    return out
end

end # module
