using Test
using RDGenomics

@testset "RDGenomics" begin

    @testset "GenomeAssembly enum" begin
        @test GRCh38 isa GenomeAssembly
        @test GRCh37 isa GenomeAssembly
        @test GRCh37 != GRCh38
    end

    @testset "HGVS coding parser" begin
        v = parse_hgvs_c("NM_000277.3:c.1241A>G")
        @test v.reference == "NM_000277.3"
        @test v.kind == :substitution
        @test v.start == 1241 && v.stop == 1241
        @test v.intronic_offset == 0
        @test v.ref == "A" && v.alt == "G"

        v = parse_hgvs_c("NM_000277.3:c.842+1G>A")
        @test v.kind == :substitution
        @test v.intronic_offset == 1

        v = parse_hgvs_c("NM_000277.3:c.842-2C>T")
        @test v.intronic_offset == -2

        v = parse_hgvs_c("NM_000277.3:c.143_144del")
        @test v.kind == :deletion
        @test (v.start, v.stop) == (143, 144)
        @test v.ref == "" && v.alt == ""

        v = parse_hgvs_c("NM_000277.3:c.143_144delAG")
        @test v.kind == :deletion && v.ref == "AG"

        v = parse_hgvs_c("NM_000277.3:c.143_144insGT")
        @test v.kind == :insertion && v.alt == "GT"

        v = parse_hgvs_c("NM_000277.3:c.143_144delinsAT")
        @test v.kind == :delins && v.alt == "AT"

        v = parse_hgvs_c("NM_000277.3:c.143dup")
        @test v.kind == :duplication

        @test_throws ArgumentError parse_hgvs_c("NM_000277.3:p.Arg408Trp")
        @test_throws ArgumentError parse_hgvs_c("no-colon-here")
        @test_throws ArgumentError parse_hgvs_c("NM_000277.3:c.junk")
        @test_throws ArgumentError parse_hgvs_c("NM_000277.3:c.5_10insAA")  # non-adjacent
    end

    @testset "HGVS protein parser" begin
        v = parse_hgvs_p("NP_000268.1:p.Arg408Trp")
        @test v.kind == :missense
        @test v.position == 408
        @test v.ref_aa == "R" && v.alt_aa == "W"

        v = parse_hgvs_p("NP_000268.1:p.R408W")
        @test v.kind == :missense && v.ref_aa == "R" && v.alt_aa == "W"

        v = parse_hgvs_p("NP_000268.1:p.Arg408*")
        @test v.kind == :nonsense && v.alt_aa == "*"

        v = parse_hgvs_p("NP_000268.1:p.Arg408Ter")
        @test v.kind == :nonsense

        v = parse_hgvs_p("NP_000268.1:p.Arg408Arg")
        @test v.kind == :synonymous

        v = parse_hgvs_p("NP_000268.1:p.Arg408del")
        @test v.kind == :deletion

        v = parse_hgvs_p("NP_000268.1:p.Arg408fs")
        @test v.kind == :frameshift

        v = parse_hgvs_p("NP_000268.1:p.(Arg408Trp)")
        @test v.kind == :missense        # predicted, parentheses tolerated

        v = parse_hgvs_p("NP_000268.1:p.Met1?")
        @test v.kind == :initiator_loss

        @test_throws ArgumentError parse_hgvs_p("NP_000268.1:c.1A>G")
        @test_throws ArgumentError parse_hgvs_p("NP_000268.1:p.Xyz408Trp")
        @test_throws ArgumentError parse_hgvs_p("NP_000268.1:p.J408W")
    end

    @testset "ACMG evidence accumulation and classification" begin
        e = AcmgEvidence()
        @test isempty(e.codes)
        push!(e, :PVS1, "canonical splice")
        push!(e, :PM2, "absent from gnomAD")
        push!(e, :PP3, "in silico support")
        @test acmg_classify(e) == :pathogenic
        @test e.rationales[:PVS1] == "canonical splice"

        e2 = AcmgEvidence()
        push!(e2, :PS1); push!(e2, :PM1); push!(e2, :PM2)
        @test acmg_classify(e2) == :likely_pathogenic

        e3 = AcmgEvidence()
        push!(e3, :BA1)
        @test acmg_classify(e3) == :benign

        e4 = AcmgEvidence()
        push!(e4, :BS1); push!(e4, :BP4)
        @test acmg_classify(e4) == :likely_benign

        e5 = AcmgEvidence()
        push!(e5, :PM2)
        @test acmg_classify(e5) == :uncertain
    end
end
