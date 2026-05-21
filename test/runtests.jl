using Test

# Top-level integration tests live here. As phases land, this file
# exercises multi-package flows (e.g. variant → gene → protein →
# AlphaFold residue → pathway → drug). For now it asserts that every
# package can be loaded together without method-table conflicts.

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

@testset "Co-loading all packages" begin
    for p in PACKAGES
        @info "Loading $p"
        @eval using $(Symbol(p))
        @test isdefined(Main, Symbol(p))
    end
end
