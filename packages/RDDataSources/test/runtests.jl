using Test
using Dates
using TOML
using RDDataSources
using RareDiseaseCore: sha256_file

# Test-only sources, defined at module scope so struct syntax is legal.

mutable struct _LocalSource <: RDDataSources.AbstractSource
    url::String
end
RDDataSources.manifest(s::_LocalSource) =
    RDDataSources.SourceManifest(
        name="LOCAL_TEST",
        urls=[s.url],
        license="CC0",
        citation="test",
    )

mutable struct _BadShaSource <: RDDataSources.AbstractSource
    url::String
    expected::String
end
RDDataSources.manifest(s::_BadShaSource) =
    RDDataSources.SourceManifest(
        name="BAD",
        urls=[s.url],
        license="CC0",
        citation="test",
        expected_sha256=[s.expected],
    )

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

    @testset "fetch! against a local file:// URL" begin
        mktempdir() do cache
            payload_dir = mktempdir()
            payload = joinpath(payload_dir, "hello.txt")
            write(payload, "rare-diseases-julia")
            url = "file://" * payload

            src = _LocalSource(url)
            files = fetch!(src; cache_dir=cache)
            @test length(files) == 1
            f = files[1]
            @test isfile(f.path)
            @test f.sha256 == sha256_file(payload)
            @test f.bytes == filesize(payload)
            @test occursin("LOCAL_TEST", f.path)

            files2 = fetch!(src; cache_dir=cache)
            @test files2[1].sha256 == f.sha256
        end
    end

    @testset "fetch! verifies sha256 when provided" begin
        mktempdir() do cache
            payload_dir = mktempdir()
            payload = joinpath(payload_dir, "v.txt")
            write(payload, "abc")
            url = "file://" * payload
            bad = _BadShaSource(url, repeat("0", 64))
            @test_throws ErrorException fetch!(bad; cache_dir=cache)
        end
    end

    @testset "update_manifest! round-trip" begin
        mktempdir() do dir
            tomlpath = joinpath(dir, "manifest.toml")
            ff = RDDataSources.FetchedFile(
                joinpath(dir, "x"), "deadbeef", "https://example/x",
                DateTime(2026, 1, 2, 3, 4, 5), 7,
            )
            update_manifest!(tomlpath, "HPO", [ff])
            data = TOML.parsefile(tomlpath)
            @test haskey(data, "sources")
            @test haskey(data["sources"], "HPO")
            entry = data["sources"]["HPO"]["files"][1]
            @test entry["sha256"] == "deadbeef"
            @test entry["url"] == "https://example/x"
            @test entry["bytes"] == 7

            ff2 = RDDataSources.FetchedFile(
                joinpath(dir, "y"), "f00", "https://example/y",
                DateTime(2026, 2, 3, 4, 5, 6), 11,
            )
            update_manifest!(tomlpath, "MONDO", [ff2])
            data2 = TOML.parsefile(tomlpath)
            @test haskey(data2["sources"], "HPO")
            @test haskey(data2["sources"], "MONDO")
        end
    end

    @testset "Unimplemented parse_source errors clearly" begin
        @test_throws ErrorException parse_source(HPOSource(), String[])
    end
end
