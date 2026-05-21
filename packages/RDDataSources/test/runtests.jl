using Test
using RDDataSources

@testset "RDDataSources" begin
    @testset "Registry" begin
        @test "HPO" in registered_sources()
    end

    @testset "HPOSource manifest" begin
        m = manifest(HPOSource())
        @test m.name == "HPO"
        @test m.license == "CC-BY-4.0"
        @test occursin("hp.obo", first(m.urls))
        @test occursin("Köhler", m.citation)
    end

    @testset "Unimplemented dispatch errors clearly" begin
        s = HPOSource()
        @test_throws ErrorException fetch!(s; cache_dir=mktempdir())
        @test_throws ErrorException parse_source(s, String[])
    end
end
