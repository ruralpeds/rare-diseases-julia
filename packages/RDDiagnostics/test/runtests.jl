using Test
using RDDiagnostics
using RDOntology
using RareDiseaseCore

const FIXTURE = joinpath(
    @__DIR__, "..", "..", "RDOntology", "test", "fixtures", "mini_hpo.obo",
)

@testset "RDDiagnostics" begin

    @testset "PatientCase construction" begin
        c = PatientCase(phenotypes_present=[HPOId("HP:0001250")])
        @test length(c.phenotypes_present) == 1
        @test isempty(c.variants)
    end

    @testset "rank_diagnoses — phenotype-only" begin
        g = load_hpo(FIXTURE)

        # Synthetic disease library aligned with the mini-HPO fixture.
        annotations = Dict(
            "DISEASE:SEIZ"  => ("Seizure-only disease",     ["HP:0001250"]),
            "DISEASE:SEIZID" => ("Seizure + ID disease",    ["HP:0001250", "HP:0001249"]),
            "DISEASE:VIS"   => ("Vision-only disease",      ["HP:0000505"]),
        )
        # IC must be computed from the same corpus we rank against.
        information_content!(g; annotations=Dict(
            d => phs for (d, (_, phs)) in annotations
        ))

        case = PatientCase(phenotypes_present=[HPOId("HP:0001250")])
        dx = rank_diagnoses(case, g; disease_annotations=annotations, topk=3)

        @test occursin("NOT FOR CLINICAL USE", dx.warning)
        @test length(dx.candidates) == 3
        # Seizure-only disease should rank above vision-only disease.
        seiz_only = findfirst(c -> c.disease == "DISEASE:SEIZ",  dx.candidates)
        vis       = findfirst(c -> c.disease == "DISEASE:VIS",   dx.candidates)
        @test seiz_only < vis    # lower index = higher rank
        @test dx.candidates[1].phenotype_score > dx.candidates[end].phenotype_score
    end

    @testset "Excluded phenotypes are filtered out of disease annotations" begin
        g = load_hpo(FIXTURE)
        annotations = Dict(
            "D1" => ("Has only seizure", ["HP:0001250"]),
        )
        information_content!(g; annotations=Dict(
            d => phs for (d, (_, phs)) in annotations
        ))
        case = PatientCase(
            phenotypes_present=[HPOId("HP:0001249")],
            phenotypes_absent=[HPOId("HP:0001250")],
        )
        dx = rank_diagnoses(case, g; disease_annotations=annotations)
        # D1's only phenotype is excluded, so it drops out.
        @test isempty(dx.candidates)
    end
end
