"""
    RDSimulation

Three-scale disease simulation: molecular (ODE / reaction networks),
physiological (PBPK / QSP), and cohort (agent-based).

Built on the SciML stack — `DifferentialEquations.jl`, `ModelingToolkit.jl`,
`Catalyst.jl`, `SBMLToolkit.jl`, `Agents.jl`. Heavy deps are added when the
ODE/PBPK/ABM bodies land (Phase 9); this stub fixes the public types and
the `RunManifest` reproducibility contract.

# source: BioModels
"""
module RDSimulation

using Dates
using Random
using RareDiseaseCore

export
    RunManifest, SimulationResult,
    MolecularModel, PBPKModel, CohortModel,
    simulate

"""
    RunManifest

Reproducibility metadata stamped on every `SimulationResult`. Code git-sha,
data hashes, RNG seed, and solver tolerances are required; anything else
goes in `extra`.
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
    SimulationResult{M}

Result of a simulation. `trajectories` is solver-specific; `citations` is
the BibTeX-ready list of every parameter source consumed.
"""
struct SimulationResult{M}
    model::M
    trajectories::Any
    citations::Vector{String}
    manifest::RunManifest
end

abstract type DiseaseModel end

struct MolecularModel <: DiseaseModel
    name::String
    species::Vector{Symbol}
    parameters::Dict{Symbol,Float64}
end

struct PBPKModel <: DiseaseModel
    name::String
    compartments::Vector{Symbol}
    parameters::Dict{Symbol,Float64}
end

struct CohortModel <: DiseaseModel
    name::String
    n_agents::Int
    transitions::Any
end

simulate(::DiseaseModel; kw...) =
    error("simulate not yet implemented (Phase 9)")

end # module
