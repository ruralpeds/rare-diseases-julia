using Test
using RDPharmacology
using RareDiseaseCore

@testset "RDPharmacology — API surface" begin
    b = Bioactivity(ChemblId("CHEMBL25"), UniProtAcc("P00439"), :IC50, 0.5, "ChEMBL")
    @test b.kind == :IC50
    @test b.value_um == 0.5
end
