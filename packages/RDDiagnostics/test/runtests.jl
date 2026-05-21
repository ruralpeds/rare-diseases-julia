using Test
using RDDiagnostics
using RareDiseaseCore

@testset "RDDiagnostics — API surface" begin
    c = PatientCase(phenotypes_present=[HPOId("HP:0001250")])
    @test length(c.phenotypes_present) == 1
    @test isempty(c.variants)
end
