using Test
using RDOntology

const FIXTURE = joinpath(@__DIR__, "fixtures", "mini_hpo.obo")

@testset "RDOntology" begin

    @testset "OBO parser" begin
        g = load_hpo(FIXTURE)
        @test length(g) == 9
        @test haskey(g, "HP:0001250")
        seizure = g["HP:0001250"]
        @test seizure.name == "Seizure"
        @test seizure.namespace == :hpo
        @test "Seizures" in seizure.synonyms
        @test "Epileptic seizure" in seizure.synonyms
        @test "UMLS:C0036572" in seizure.xrefs
        @test g["HP:9999999"].is_obsolete
    end

    @testset "Ancestors / descendants / is_a" begin
        g = load_hpo(FIXTURE)
        anc = ancestors(g, "HP:0001250")
        @test "HP:0012638" in anc
        @test "HP:0000707" in anc
        @test "HP:0000118" in anc
        @test "HP:0000001" in anc
        @test !("HP:0000478" in anc)            # different subtree

        dec = descendants(g, "HP:0000118")
        @test "HP:0001250" in dec
        @test "HP:0000505" in dec
        @test !("HP:0000001" in dec)

        @test is_a(g, "HP:0001250", "HP:0000001")
        @test is_a(g, "HP:0001250", "HP:0001250")
        @test !is_a(g, "HP:0001250", "HP:0000478")
        @test !is_a(g, "HP:nonexistent", "HP:0000001")
    end

    @testset "Cross-reference resolution" begin
        g = load_hpo(FIXTURE)
        @test resolve_xref(g, "UMLS:C0036572") == "HP:0001250"
        @test resolve_xref(g, "UMLS:nonexistent") === nothing
    end

    @testset "Information content + Resnik similarity" begin
        g = load_hpo(FIXTURE)
        # Two synthetic diseases share a seizure phenotype.
        annotations = Dict(
            "DISEASE:1" => ["HP:0001250", "HP:0001249"],   # seizure + ID
            "DISEASE:2" => ["HP:0001250"],                  # seizure only
            "DISEASE:3" => ["HP:0000505"],                  # visual impairment
        )
        information_content!(g; annotations=annotations)

        # The root has p=1 → IC=0; specific terms have higher IC.
        @test term_ic(g, "HP:0000001") ≈ 0.0 atol=1e-9
        @test term_ic(g, "HP:0001250") > term_ic(g, "HP:0000118")
        @test term_ic(g, "HP:0000505") > term_ic(g, "HP:0000478")

        # Seizure ↔ seizure → MICA = seizure itself
        mica_id, mica_ic = most_informative_common_ancestor(
            g, "HP:0001250", "HP:0001250"
        )
        @test mica_id == "HP:0001250"
        @test mica_ic == term_ic(g, "HP:0001250")

        # Seizure ↔ intellectual disability share their nervous-system parent
        mica_id, _ = most_informative_common_ancestor(
            g, "HP:0001250", "HP:0001249"
        )
        @test mica_id == "HP:0012638"

        # Seizure ↔ visual impairment share only "Phenotypic abnormality"
        mica_id, _ = most_informative_common_ancestor(
            g, "HP:0001250", "HP:0000505"
        )
        @test mica_id == "HP:0000118"

        # Best-pairs Resnik: query overlapping with itself ≥ query vs unrelated
        sim_self  = phenotype_similarity(g,
            ["HP:0001250"], ["HP:0001250"]; method=:resnik)
        sim_cross = phenotype_similarity(g,
            ["HP:0001250"], ["HP:0000505"]; method=:resnik)
        @test sim_self > sim_cross

        # Lin similarity is bounded in [0, 1]
        sim_lin = phenotype_similarity(g,
            ["HP:0001250"], ["HP:0001249"]; method=:lin)
        @test 0.0 ≤ sim_lin ≤ 1.0

        # Jiang-Conrath similarity is bounded in (0, 1]
        sim_jc = phenotype_similarity(g,
            ["HP:0001250"], ["HP:0001250"]; method=:jc)
        @test 0.0 < sim_jc ≤ 1.0
    end

    @testset "Errors before IC is computed" begin
        g = load_hpo(FIXTURE)
        @test_throws ErrorException term_ic(g, "HP:0001250")
        @test_throws ErrorException phenotype_similarity(
            g, ["HP:0001250"], ["HP:0001250"]
        )
    end
end
