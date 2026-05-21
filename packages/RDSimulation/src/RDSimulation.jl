"""
    RDSimulation

Three-scale disease simulation: molecular (ODE / reaction networks),
physiological (PBPK / QSP), and cohort (agent-based).

Phase 9 of the build plan. To keep CI lean and the dependency graph small,
this package ships a hand-rolled explicit RK4 integrator and a worked
PKU (phenylalanine hydroxylase) kinetic model. Heavier solvers from the
SciML stack can be added later as optional extensions.

# source: BioModels
"""
module RDSimulation

using Dates
using Random
using LinearAlgebra
using Statistics
using RareDiseaseCore

include("integrator.jl")
include("models/pah_pku.jl")
include("cohort.jl")

export
    RunManifest, SimulationResult,
    MolecularModel, PBPKModel, CohortModel, DiseaseModel,
    simulate, rk4,
    # PKU model exports
    PAHParameters, default_pah_parameters,
    variant_effect, pah_residual_activity,
    # Cohort
    CohortAgent, cohort_simulate

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

Result of a simulation. `t` and `u` are arrays; `citations` is the
BibTeX-ready list of every parameter source consumed.
"""
struct SimulationResult{M}
    model::M
    t::Vector{Float64}
    u::Matrix{Float64}                 # rows = species, cols = timepoints
    species::Vector{Symbol}
    citations::Vector{String}
    manifest::RunManifest
end

abstract type DiseaseModel end

struct MolecularModel <: DiseaseModel
    name::String
    species::Vector{Symbol}
    parameters::Dict{Symbol,Float64}
    rhs::Function                       # (du, u, p, t) -> nothing
end

struct PBPKModel <: DiseaseModel
    name::String
    compartments::Vector{Symbol}
    parameters::Dict{Symbol,Float64}
    rhs::Function
end

struct CohortModel <: DiseaseModel
    name::String
    n_agents::Int
    transitions::Any
end

"""
    simulate(model; u0, tspan, dt=0.01, manifest) -> SimulationResult

Integrate `model` from `tspan[1]` to `tspan[2]` with initial condition
`u0` using explicit RK4. `manifest` is the reproducibility stamp.
"""
function simulate(
    model::M;
    u0::AbstractVector{<:Real},
    tspan::Tuple{<:Real,<:Real},
    dt::Real=0.01,
    citations::Vector{String}=String[],
    manifest::RunManifest,
) where {M<:Union{MolecularModel,PBPKModel}}
    t, u = rk4(model.rhs, collect(Float64, u0), model.parameters,
               Float64(tspan[1]), Float64(tspan[2]), Float64(dt))
    species = M === MolecularModel ? model.species : model.compartments
    return SimulationResult{M}(model, t, u, species, citations, manifest)
end

end # module
