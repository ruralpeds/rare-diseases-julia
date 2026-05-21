using Test
using RDApp
using RDApp: handle_diagnose, handle_simulate, handle_treatments,
             build_default_state

@testset "RDApp" begin
    @testset "Route table" begin
        rs = routes()
        @test any(r -> r.path == "/health" && r.method == "GET", rs)
        @test any(r -> r.path == "/diagnose" && r.method == "POST", rs)
        @test any(r -> r.path == "/simulate" && r.method == "POST", rs)
        @test any(r -> r.path == "/treatments" && r.method == "POST", rs)
        @test all(r -> haskey(r, :desc) && !isempty(r.desc), rs)
    end

    @testset "Banner contents" begin
        @test occursin("research", lowercase(RDApp.NOT_FOR_CLINICAL_USE_BANNER))
        @test occursin("not for clinical use",
                       lowercase(RDApp.NOT_FOR_CLINICAL_USE_BANNER))
    end

    @testset "Default app state loads bundled fixtures" begin
        st = build_default_state()
        @test length(st.ontology) > 0
        @test !isempty(st.disease_annotations)
        @test !isempty(st.drugs)
    end

    @testset "handle_diagnose against bundled state" begin
        st = build_default_state()
        body = Dict("phenotypes_present" => ["HP:0001250"])
        out = handle_diagnose(st, body)
        @test occursin("NOT FOR CLINICAL USE", out["warning"])
        @test !isempty(out["candidates"])
        # The seizure-only disease should top the list for a seizure case
        @test first(out["candidates"])["disease"] == "DISEASE:SEIZ"
    end

    @testset "handle_simulate dispatches to the right model" begin
        out = handle_simulate(Dict(
            "model" => "PAH_PKU",
            "params" => Dict("variant" => "classical_pku", "bh4" => 2.0),
        ))
        @test out["model"] == "PAH_PKU"
        @test haskey(out["final"], "Phe")
        @test out["final"]["Phe"] > 0

        out2 = handle_simulate(Dict(
            "model" => "CFTR_CF",
            "params" => Dict("class" => "II",
                             "tezacaftor" => true,
                             "elexacaftor" => true,
                             "ivacaftor" => true),
        ))
        @test out2["model"] == "CFTR_CF"
        @test haskey(out2["final"], "ASL")

        out3 = handle_simulate(Dict(
            "model" => "HBS_SCD",
            "params" => Dict("hbf_fraction" => 0.15),
        ))
        @test haskey(out3["final"], "Poly")

        bad = handle_simulate(Dict("model" => "nonexistent"))
        @test occursin("unknown model", bad["error"])
    end

    @testset "handle_treatments returns ranked candidates" begin
        st = build_default_state()
        out = handle_treatments(st, Dict("disease_genes" => ["SCN1A"]))
        @test occursin("NOT FOR CLINICAL USE", out["warning"])
        @test !isempty(out["candidates"])
        # First-ranked drug should be lamotrigine (approved, on-target)
        @test first(out["candidates"])["drug"] == "lamotrigine"
    end

    @testset "SOURCES.tsv loader" begin
        rows = RDApp._load_sources_tsv()
        @test rows isa AbstractVector
        if !isempty(rows)
            @test haskey(first(rows), "name")
        end
    end
end
