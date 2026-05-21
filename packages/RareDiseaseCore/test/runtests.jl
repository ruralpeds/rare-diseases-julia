using Test
using Dates
using RareDiseaseCore

@testset "RareDiseaseCore" begin

    @testset "Identifier parsing — happy path" begin
        @test string(HPOId("HP:0001250")) == "HP:0001250"
        @test string(MondoId("MONDO:0009861")) == "MONDO:0009861"
        @test string(OrphaId("ORPHA:716")) == "ORPHA:716"
        @test string(OmimId("OMIM:261600")) == "OMIM:261600"
        @test string(DoId("DOID:9281")) == "DOID:9281"
        @test string(MeshId("D010661")) == "D010661"
        @test string(HgncId("HGNC:8582")) == "HGNC:8582"
        @test string(EnsemblGeneId("ENSG00000171759")) == "ENSG00000171759"
        @test string(EnsemblGeneId("ENSG00000171759.12")) == "ENSG00000171759.12"
        @test string(EnsemblTranscriptId("ENST00000553106")) == "ENST00000553106"
        @test string(RefSeqId("NM_000277.3")) == "NM_000277.3"
        @test string(UniProtAcc("P00439")) == "P00439"
        @test string(ClinVarId("12345")) == "12345"
        @test string(DbSnpId("rs5030858")) == "rs5030858"
        @test string(PubMedId("33264411")) == "33264411"
        @test string(PmcId("PMC7778952")) == "PMC7778952"
        @test string(RxCui("8629")) == "8629"
        @test string(ChemblId("CHEMBL25")) == "CHEMBL25"
        @test string(PubChemCid("2244")) == "2244"
        @test string(ReactomeId("R-HSA-71291")) == "R-HSA-71291"
        @test string(PfamId("PF00351")) == "PF00351"
        @test string(InterProId("IPR019773")) == "IPR019773"
    end

    @testset "Identifier parsing — rejection" begin
        @test_throws ArgumentError HPOId("")
        @test_throws ArgumentError HPOId("0001250")            # missing prefix
        @test_throws ArgumentError HPOId("HP:abc")             # non-digit
        @test_throws ArgumentError MondoId("MONDO:1")          # too short
        @test_throws ArgumentError UniProtAcc("not-an-acc")
        @test_throws ArgumentError EnsemblGeneId("ENSG123")    # wrong length
        @test_throws ArgumentError RefSeqId("ZZ_000001")
        @test_throws ArgumentError DbSnpId("12345")            # missing rs prefix
        @test_throws ArgumentError ClinVarId("abc")
    end

    @testset "Identifier equality + hashing" begin
        a = HPOId("HP:0001250")
        b = HPOId("HP:0001250")
        c = HPOId("HP:0001251")
        @test a == b
        @test a != c
        @test hash(a) == hash(b)
        @test hash(a) != hash(c)
        # Type-disjoint: same string body should not hash-collide across types
        @test hash(HPOId("HP:0000001")) != hash(MondoId("MONDO:0000001"))
        # Set membership
        s = Set([HPOId("HP:0001250"), HPOId("HP:0001250")])
        @test length(s) == 1
    end

    @testset "DataProvenance + cite" begin
        prov = DataProvenance(
            "HPO", "2024-12-12", "deadbeef", DateTime(2025, 1, 2, 3, 4, 5);
            url="https://hpo.jax.org",
            citation="Köhler S, et al. NAR 2021. PMID:33264411",
        )
        s = cite(prov)
        @test occursin("Köhler", s)
        @test occursin("HPO", s)
        @test occursin("2025-01-02", s)
    end

    @testset "Record construction" begin
        prov = DataProvenance("HPO", "v1", "0"^64, now())
        ph = Phenotype(HPOId("HP:0001250"), "Seizure", nothing, prov)
        @test ph.name == "Seizure"
        @test string(ph.id) == "HP:0001250"
    end

    @testset "sha256_file" begin
        mktempdir() do dir
            path = joinpath(dir, "x.txt")
            write(path, "rare-diseases-julia")
            h = sha256_file(path)
            @test length(h) == 64
            @test all(c -> c in "0123456789abcdef", h)
            # Stable across calls
            @test h == sha256_file(path)
        end
    end
end
