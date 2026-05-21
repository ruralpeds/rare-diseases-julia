"""
    RDSimulation

Three-scale disease simulation, built on the Julia ecosystem:

* Molecular / reaction networks — `Catalyst.@reaction_network`
* ODE integration              — `OrdinaryDiffEq` (`Tsit5`, etc.)
* Agent-based cohorts          — `Agents.jl`

Disease-specific contributions in this package: parameter tables for
real disease variants (e.g. PAH residual-activity classes), worked
reaction networks for example diseases, and a `RunManifest`
reproducibility stamp.

Phase 9 of the build plan.

# source: BioModels
"""
module RDSimulation

using Dates
using Random
using Agents
using Catalyst
using ModelingToolkit
using OrdinaryDiffEq
using SBML
using SBMLToolkit
using RareDiseaseCore

include("models/pah_pku.jl")
include("models/sapropterin_pbpk.jl")
include("biomodels.jl")
include("cohort.jl")

export
    # Reproducibility
    RunManifest, SimulationResult,
    # PKU model + helpers
    PAH_PKU, default_pah_parameters,
    pah_residual_activity, variant_effect,
    pah_pku_problem,
    # Sapropterin PBPK
    SAPROPTERIN_PBPK, sapropterin_pbpk_problem,
    # BioModels
    load_biomodel,
    # Cohort
    Patient, build_cohort_model, run_cohort!

"""
    RunManifest

Reproducibility metadata stamped on every `SimulationResult`. Code git
sha, data hashes, RNG seed, and solver tolerances are required; anything
else goes in `extra`.
"""
Base.@kwdef struct RunManifest
    code_git_sha::String
    data_hashes::Dict{String,String}
    rng_seed::UInt64
    solver::String
    abstol::Float64
    reltol::Float64
    created_at::DateTime = now()
    extra::Dict{String,Any} = Dict{String,Any}()
end

"""
    SimulationResult

Wraps an `OrdinaryDiffEq.ODESolution` together with citations and a
`RunManifest`. Callers can still operate on `result.sol` directly to
get full SciML semantics (interpolation, sensitivity, etc.).
"""
struct SimulationResult{S}
    sol::S
    citations::Vector{String}
    manifest::RunManifest
end

end # module
