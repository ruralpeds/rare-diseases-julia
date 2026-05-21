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
using OrdinaryDiffEq
using Oxygen

using RareDiseaseCore
using RDOntology
using RDDiagnostics
using RDPathways
using RDTreatment
using RDSimulation

export
    routes, register_routes!, start_server,
    AppState, build_default_state,
    NOT_FOR_CLINICAL_USE_BANNER,
    handle_diagnose, handle_simulate, handle_treatments

const NOT_FOR_CLINICAL_USE_BANNER =
    "This service is for research and education only. Not for clinical use."

"""
    AppState

Server-side state held for the duration of a server process: the loaded
ontology graph, the bundled disease annotations used by `/diagnose`,
and the pathway network and drug library used by `/treatments`.
"""
mutable struct AppState
    ontology::OntologyGraph
    disease_annotations::Dict{String,Tuple{String,Vector{String}}}
    network::PathwayNetwork
    drugs::Vector{DrugRecord}
end

"""
    routes() -> Vector{NamedTuple}

Documented route table; single source of truth for client SDK generation.
"""
function routes()
    return [
        (method="GET",  path="/health",         desc="Liveness probe."),
        (method="GET",  path="/sources",        desc="Open data sources."),
        (method="GET",  path="/disease/{mondo}",desc="Disease record."),
        (method="GET",  path="/variant/{clinvar}", desc="Variant record."),
        (method="GET",  path="/protein/{uniprot}",  desc="Protein record."),
        (method="GET",  path="/pathway/{reactome}", desc="Pathway record."),
        (method="POST", path="/diagnose",       desc="Differential diagnosis."),
        (method="POST", path="/simulate",       desc="Run a registered model."),
        (method="POST", path="/treatments",     desc="Rank treatments for disease."),
    ]
end

# ---------------------------------------------------------------------------
# Default bundled state (mini fixtures so the server boots without
# external downloads)
# ---------------------------------------------------------------------------

"""
    build_default_state() -> AppState

Construct a small in-memory `AppState` from bundled fixtures. The
ontology is the mini HPO under `packages/RDOntology/test/fixtures/`;
diseases and the pathway network are toy data sized for demos and tests.
"""
function build_default_state()
    fixture = joinpath(@__DIR__, "..", "..", "RDOntology",
                       "test", "fixtures", "mini_hpo.obo")
    g = load_hpo(fixture)

    annotations = Dict(
        "DISEASE:SEIZ" => ("Seizure-only disease",      ["HP:0001250"]),
        "DISEASE:SEIZID" => ("Seizure + ID disease",   ["HP:0001250", "HP:0001249"]),
        "DISEASE:VIS"  => ("Vision-only disease",       ["HP:0000505"]),
    )
    information_content!(g; annotations=Dict(
        d => phs for (d, (_, phs)) in annotations
    ))

    net = PathwayNetwork()
    add_edge_undirected!(net, "SCN1A", "GABA_R")
    add_edge_undirected!(net, "GABA_R", "GAD1")
    add_edge_undirected!(net, "GAD1", "HMGCR")
    add_edge_undirected!(net, "HMGCR", "LDLR")

    drugs = [
        DrugRecord(name="lamotrigine", targets=["SCN1A"],
                   evidence_tier=APPROVED),
        DrugRecord(name="atorvastatin", targets=["HMGCR"],
                   evidence_tier=APPROVED),
        DrugRecord(name="experimental_gaba",
                   targets=["GABA_R"],
                   evidence_tier=REPURPOSING_HYPOTHESIS),
    ]
    return AppState(g, annotations, net, drugs)
end

# ---------------------------------------------------------------------------
# Pure handlers — easy to unit test without spinning up the HTTP server
# ---------------------------------------------------------------------------

"""
    handle_diagnose(state, body) -> Dict

Body shape: `{ "phenotypes_present": ["HP:..."], "phenotypes_absent": ["HP:..."] }`.
Returns ranked candidates and the not-for-clinical-use warning.
"""
function handle_diagnose(state::AppState, body::AbstractDict)
    present = String.(get(body, "phenotypes_present", String[]))
    absent  = String.(get(body, "phenotypes_absent",  String[]))
    case = PatientCase(
        phenotypes_present = [HPOId(p) for p in present],
        phenotypes_absent  = [HPOId(p) for p in absent],
    )
    dx = rank_diagnoses(case, state.ontology;
                        disease_annotations=state.disease_annotations,
                        topk=10)
    return Dict(
        "banner"  => NOT_FOR_CLINICAL_USE_BANNER,
        "warning" => dx.warning,
        "candidates" => [
            Dict(
                "disease" => c.disease,
                "name"    => c.name,
                "score"   => c.score,
                "phenotype_score" => c.phenotype_score,
                "variant_score"   => c.variant_score,
                "evidence" => c.evidence,
            )
            for c in dx.candidates
        ],
    )
end

"""
    handle_simulate(body) -> Dict

Body shape:
```
{ "model": "PAH_PKU" | "CFTR_CF" | "HBS_SCD",
  "params": { ... model-specific ... } }
```
Returns final-state values and a tiny metadata block.
"""
function handle_simulate(body::AbstractDict)
    model = String(get(body, "model", "PAH_PKU"))
    params = get(body, "params", Dict{String,Any}())
    prob, species = if model == "PAH_PKU"
        variant = Symbol(get(params, "variant", "wildtype"))
        bh4     = Float64(get(params, "bh4", 1.0))
        pah_pku_problem(; variant=variant, bh4=bh4), [:Phe, :Tyr]
    elseif model == "CFTR_CF"
        cls = Symbol(get(params, "class", "wildtype"))
        mods = (;
            ivacaftor   = Bool(get(params, "ivacaftor",   false)),
            tezacaftor  = Bool(get(params, "tezacaftor",  false)),
            elexacaftor = Bool(get(params, "elexacaftor", false)),
        )
        cftr_cf_problem(; class=cls, modulators=mods), [:ASL]
    elseif model == "HBS_SCD"
        hbf = Float64(get(params, "hbf_fraction", 0.01))
        hbs_scd_problem(; hbf_fraction=hbf), [:Mono, :Poly]
    else
        return Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
                    "error"=>"unknown model '$model'")
    end
    sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
    return Dict(
        "banner"  => NOT_FOR_CLINICAL_USE_BANNER,
        "model"   => model,
        "species" => string.(species),
        "final"   => Dict(string(s) => sol.u[end][i]
                          for (i, s) in enumerate(species)),
        "t_end"   => sol.t[end],
    )
end

"""
    handle_treatments(state, body) -> Dict

Body shape: `{ "disease_genes": ["SCN1A", ...] }`.
"""
function handle_treatments(state::AppState, body::AbstractDict)
    disease_genes = String.(get(body, "disease_genes", String[]))
    cands, warning = rank_treatments(disease_genes, state.drugs, state.network;
                                      n_bootstrap=30)
    return Dict(
        "banner"  => NOT_FOR_CLINICAL_USE_BANNER,
        "warning" => warning,
        "candidates" => [
            Dict(
                "drug"     => c.drug.name,
                "score"    => c.score,
                "proximity_d" => c.proximity_d,
                "proximity_z" => c.proximity_z,
                "tier"     => string(c.drug.evidence_tier),
                "rationale" => c.mechanism_rationale,
            )
            for c in cands
        ],
    )
end

# ---------------------------------------------------------------------------
# Oxygen wiring
# ---------------------------------------------------------------------------

"""
    register_routes!(state::AppState; oxygen=Oxygen) -> nothing
"""
function register_routes!(state::AppState; oxygen::Module=Oxygen)
    oxygen.@get "/health" function (req)
        Dict("status" => "ok", "banner" => NOT_FOR_CLINICAL_USE_BANNER)
    end
    oxygen.@get "/sources" function (req)
        Dict("banner"=>NOT_FOR_CLINICAL_USE_BANNER,
             "sources"=>_load_sources_tsv())
    end
    oxygen.@post "/diagnose" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)), Dict{String,Any})
        return handle_diagnose(state, body)
    end
    oxygen.@post "/simulate" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)), Dict{String,Any})
        return handle_simulate(body)
    end
    oxygen.@post "/treatments" function (req)
        body = JSON3.read(IOBuffer(HTTP.payload(req)), Dict{String,Any})
        return handle_treatments(state, body)
    end
    return nothing
end

"""
    start_server(; host="127.0.0.1", port=8000, state=build_default_state())

Register routes and start Oxygen's HTTP server.
"""
function start_server(; host::String="127.0.0.1", port::Int=8000,
                      state::AppState=build_default_state())
    register_routes!(state)
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
