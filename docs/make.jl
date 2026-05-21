using Documenter

# Dev-in all in-tree packages before loading them for docs.
using Pkg
for p in readdir(joinpath(@__DIR__, "..", "packages"); join=false)
    path = joinpath(@__DIR__, "..", "packages", p)
    Pkg.develop(PackageSpec(path=path))
end
Pkg.instantiate()

using RareDiseaseCore
using RDDataSources
using RDOntology
using RDGenomics
using RDProteomics
using RDPathways
using RDPharmacology
using RDClinical
using RDDiagnostics
using RDSimulation
using RDTreatment
using RDApp

makedocs(;
    sitename = "rare-diseases-julia",
    authors  = "rare-diseases-julia contributors",
    modules  = [
        RareDiseaseCore, RDDataSources, RDOntology, RDGenomics,
        RDProteomics, RDPathways, RDPharmacology, RDClinical,
        RDDiagnostics, RDSimulation, RDTreatment, RDApp,
    ],
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link  = "main",
        canonical  = "https://ruralpeds.github.io/rare-diseases-julia/",
        assets     = String[],
    ),
    pages = [
        "Home"            => "index.md",
        "Plan"            => "plan.md",
        "Ethics"          => "ethics.md",
        "Packages" => [
            "RareDiseaseCore" => "packages/raredisease_core.md",
            "RDDataSources"   => "packages/rd_data_sources.md",
            "RDOntology"      => "packages/rd_ontology.md",
            "RDGenomics"      => "packages/rd_genomics.md",
            "RDProteomics"    => "packages/rd_proteomics.md",
            "RDPathways"      => "packages/rd_pathways.md",
            "RDPharmacology"  => "packages/rd_pharmacology.md",
            "RDClinical"      => "packages/rd_clinical.md",
            "RDDiagnostics"   => "packages/rd_diagnostics.md",
            "RDSimulation"    => "packages/rd_simulation.md",
            "RDTreatment"     => "packages/rd_treatment.md",
            "RDApp"           => "packages/rd_app.md",
        ],
    ],
    warnonly = true,    # don't fail the build on cross-ref misses while
                        # the API surface is still solidifying
)

deploydocs(;
    repo = "github.com/ruralpeds/rare-diseases-julia.git",
    devbranch = "main",
    push_preview = true,
)
