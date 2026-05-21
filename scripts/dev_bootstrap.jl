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
    "RDPathways",
    "RDPharmacology",
    "RDClinical",
    "RDDiagnostics",
    "RDSimulation",
    "RDTreatment",
    "RDApp",
]

for p in PACKAGES
    path = joinpath(@__DIR__, "..", "packages", p)
    @info "dev $p" path
    Pkg.develop(PackageSpec(path=path))
end

Pkg.instantiate()
Pkg.precompile()
