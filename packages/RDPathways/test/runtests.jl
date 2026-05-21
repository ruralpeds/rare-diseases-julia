using Test
using RDPathways

@testset "RDPathways — API surface" begin
    n = PathwayNetwork()
    @test isempty(n.pathways)
    @test_throws ErrorException load_reactome("/nonexistent")
end
