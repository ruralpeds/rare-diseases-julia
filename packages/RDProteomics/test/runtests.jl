using Test
using RDProteomics
using RareDiseaseCore

@testset "RDProteomics — API surface" begin
    r = ResidueFeature(UniProtAcc("P00439"), 408, 95.2, :helix, ["PF00351"], 3.4)
    @test r.position == 408
    @test r.secondary == :helix
    @test_throws ErrorException load_uniprot("/nonexistent")
end
