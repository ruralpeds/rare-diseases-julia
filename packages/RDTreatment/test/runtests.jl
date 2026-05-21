using Test
using RDTreatment
using RareDiseaseCore

@testset "RDTreatment — API surface" begin
    c = TreatmentCandidate(
        "sapropterin",
        ChemblId("CHEMBL1201824"),
        RxCui("714583"),
        APPROVED,
        [UniProtAcc("P00439")],
        0.92,
        "BH4 cofactor restoration in residual-activity PAH variants",
        ["PMID:18762447"],
    )
    @test c.evidence_tier == APPROVED
    @test 0 ≤ c.score ≤ 1
end
