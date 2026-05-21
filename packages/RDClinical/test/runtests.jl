using Test
using RDClinical

@testset "RDClinical — API surface" begin
    s = AdverseEventSignal("sapropterin", "headache", 42, 1.2, 1.3)
    @test s.n_reports == 42
    @test s.prr ≈ 1.2
end
