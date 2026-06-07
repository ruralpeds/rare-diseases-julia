# scripts/dev_bootstrap.jl
#
# Develop every in-tree package into the meta-environment so that
# `julia --project=.` can `using` any of them.
#
# Usage:
#   julia --project=. scripts/dev_bootstrap.jl

using Pkg

const PACKAGES = [
    "RareDiseaseCore",
    "RDDataSources",
    "RDOntology",
    "RDGenomics",
    "RDProteomics",
    "RDImmunology",
    "RDPathways",
    "RDPharmacology",
    "RDClinical",
    "RDDiagnostics",
    "RDSimulation",
    "RDTreatment",
    "RDApp",
]

paths = String[]
for p in PACKAGES
    path = joinpath(@__DIR__, "..", "packages", p)
    @info "queue dev $p" path
    push!(paths, path)
end
Pkg.develop([PackageSpec(; path=p) for p in paths])

Pkg.instantiate()
Pkg.precompile()
