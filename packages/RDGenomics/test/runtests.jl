using Test
using RDGenomics

@testset "RDGenomics — API surface" begin
    @test GRCh38 isa GenomeAssembly
    @test GRCh37 isa GenomeAssembly
    e = AcmgEvidence()
    @test isempty(e.codes)
    @test isempty(e.rationales)
end
