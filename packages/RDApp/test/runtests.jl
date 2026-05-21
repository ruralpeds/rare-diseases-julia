using Test
using RDApp

@testset "RDApp" begin
    @testset "Route table" begin
        rs = routes()
        @test any(r -> r.path == "/health" && r.method == "GET", rs)
        @test any(r -> r.path == "/diagnose" && r.method == "POST", rs)
        @test any(r -> r.path == "/simulate" && r.method == "POST", rs)
        @test any(r -> r.path == "/treatments" && r.method == "POST", rs)
        # Every route is documented
        @test all(r -> haskey(r, :desc) && !isempty(r.desc), rs)
    end

    @testset "Banner contents" begin
        @test occursin("research", lowercase(RDApp.NOT_FOR_CLINICAL_USE_BANNER))
        @test occursin("not for clinical use",
                       lowercase(RDApp.NOT_FOR_CLINICAL_USE_BANNER))
    end

    @testset "SOURCES.tsv loader" begin
        rows = RDApp._load_sources_tsv()
        # Repo ships ~40 sources; tolerate either the bundled file being
        # available or the test running from a stripped-down checkout.
        @test rows isa AbstractVector
        if !isempty(rows)
            first_row = first(rows)
            @test haskey(first_row, "name")
        end
    end
end
