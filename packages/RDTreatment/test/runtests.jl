using Test
using Random
using RareDiseaseCore
using RDPathways
using RDTreatment

@testset "RDTreatment" begin

    @testset "TreatmentCandidate basic construction" begin
        d = DrugRecord(name="sapropterin",
                       chembl=ChemblId("CHEMBL1201824"),
                       rxcui=RxCui("714583"),
                       targets=["PAH"],
                       evidence_tier=APPROVED,
                       citations=["PMID:18762447"])
        @test d.name == "sapropterin"
        @test d.evidence_tier == APPROVED
    end

    @testset "rank_treatments uses pathway proximity + tier bonus" begin
        # Build a tiny pathway network:
        #   PAH -- TYR_PATH -- TH
        #   off-target node OFF connected far away
        n = PathwayNetwork()
        add_edge_undirected!(n, "PAH", "TYR_PATH")
        add_edge_undirected!(n, "TYR_PATH", "TH")
        add_edge_undirected!(n, "FAR1", "FAR2")
        add_edge_undirected!(n, "FAR2", "FAR3")
        add_edge_undirected!(n, "FAR3", "OFF")
        # Connect the two components weakly so closest_distance is finite
        add_edge_undirected!(n, "TH", "FAR1")

        disease = ["PAH"]
        on_target = DrugRecord(name="on",
                               targets=["TYR_PATH"],
                               evidence_tier=REPURPOSING_HYPOTHESIS)
        off_target = DrugRecord(name="off",
                                targets=["OFF"],
                                evidence_tier=REPURPOSING_HYPOTHESIS)
        approved_far = DrugRecord(name="approved-but-far",
                                  targets=["OFF"],
                                  evidence_tier=APPROVED)

        cands, warning = rank_treatments(
            disease, [on_target, off_target, approved_far], n;
            n_bootstrap=30,
            rng=Xoshiro(0),
        )

        @test occursin("NOT FOR CLINICAL USE", warning)
        # All three drugs reachable -> all three candidates
        @test length(cands) == 3

        names = [c.drug.name for c in cands]
        # Approved drug should rank first thanks to tier bonus,
        # even though its target is far from disease.
        @test names[1] == "approved-but-far"
        # On-target unapproved should outrank off-target unapproved.
        @test findfirst(==("on"),  names) <
              findfirst(==("off"), names)

        # Proximity distance must reflect on-target advantage
        on_cand  = cands[findfirst(c -> c.drug.name == "on",  cands)]
        off_cand = cands[findfirst(c -> c.drug.name == "off", cands)]
        @test on_cand.proximity_d < off_cand.proximity_d
    end

    @testset "Drugs with no reachable target are dropped" begin
        n = PathwayNetwork()
        add_edge_undirected!(n, "A", "B")
        add_node!(n, "ISOLATED")

        d = DrugRecord(name="isolated", targets=["ISOLATED"])
        cands, _ = rank_treatments(["A"], [d], n; n_bootstrap=5)
        @test isempty(cands)
    end
end
