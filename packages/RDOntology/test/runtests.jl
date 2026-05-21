using Test
using RDOntology

@testset "RDOntology — API surface" begin
    g = OntologyGraph()
    @test isempty(g.terms)
    @test_throws ErrorException load_hpo("/nonexistent")
    @test_throws ErrorException phenotype_similarity(g, ["HP:0001250"], ["HP:0001250"])
end
