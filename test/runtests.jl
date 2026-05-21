using Test

# Top-level integration tests live here. As phases land, this file
# exercises multi-package flows.

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

@testset "Cross-package integration" begin

    @testset "Co-loading all packages" begin
        for p in PACKAGES
            @info "Loading $p"
            @eval using $(Symbol(p))
            @test isdefined(Main, Symbol(p))
        end
    end

    @testset "ontology → diagnostics → treatment flow" begin
        using RareDiseaseCore
        using RDOntology
        using RDDiagnostics
        using RDPathways
        using RDTreatment
        using Random

        # 1. Load mini HPO and compute information content from a tiny
        #    annotation corpus.
        fixture = joinpath(@__DIR__, "..", "packages", "RDOntology",
                           "test", "fixtures", "mini_hpo.obo")
        g = load_hpo(fixture)
        annotations = Dict(
            "DISEASE:SEIZ" => ("Seizure-only disease",   ["HP:0001250"]),
            "DISEASE:VIS"  => ("Vision-only disease",    ["HP:0000505"]),
        )
        information_content!(g; annotations=Dict(
            d => phs for (d, (_, phs)) in annotations
        ))

        # 2. Rank diagnoses for a seizure patient.
        case = PatientCase(phenotypes_present=[HPOId("HP:0001250")])
        dx = rank_diagnoses(case, g; disease_annotations=annotations, topk=5)
        @test !isempty(dx.candidates)
        @test occursin("NOT FOR CLINICAL USE", dx.warning)
        @test dx.candidates[1].disease == "DISEASE:SEIZ"

        # 3. For the top candidate, build a toy pathway network and rank
        #    candidate treatments. The seizure-disease "module" is a
        #    single placeholder gene SCN1A; sodium-channel blockers
        #    target it directly while statins target a far node.
        net = PathwayNetwork()
        add_edge_undirected!(net, "SCN1A", "GABA_R")
        add_edge_undirected!(net, "GABA_R", "GAD1")
        # Distant node + connector so closest_distance stays finite.
        add_edge_undirected!(net, "GAD1", "HMGCR")
        add_edge_undirected!(net, "HMGCR", "LDLR")

        drugs = [
            DrugRecord(name="lamotrigine",
                       targets=["SCN1A"],
                       evidence_tier=APPROVED),
            DrugRecord(name="atorvastatin",
                       targets=["HMGCR"],
                       evidence_tier=APPROVED),
            DrugRecord(name="some_repurposing",
                       targets=["GABA_R"],
                       evidence_tier=REPURPOSING_HYPOTHESIS),
        ]
        cands, warning = rank_treatments(
            ["SCN1A"], drugs, net;
            n_bootstrap=30, rng=Xoshiro(0),
        )
        @test occursin("NOT FOR CLINICAL USE", warning)
        @test !isempty(cands)
        # The on-target approved drug must rank above the off-target
        # approved drug (same tier bonus, better proximity).
        names = [c.drug.name for c in cands]
        @test findfirst(==("lamotrigine"), names) <
              findfirst(==("atorvastatin"), names)
    end
end
