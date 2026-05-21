using Test
using RDApp

@testset "RDApp — route table" begin
    rs = routes()
    @test any(r -> r.path == "/health", rs)
    @test any(r -> r.path == "/diagnose" && r.method == "POST", rs)
    @test any(r -> r.path == "/simulate" && r.method == "POST", rs)
    @test occursin("research", lowercase(RDApp.NOT_FOR_CLINICAL_USE_BANNER))
end
