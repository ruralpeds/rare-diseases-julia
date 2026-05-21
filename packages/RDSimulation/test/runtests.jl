using Test
using Dates
using RDSimulation

@testset "RDSimulation — API surface" begin
    m = MolecularModel("PKU-PAH", [:Phe, :Tyr], Dict(:kcat => 1.2, :Km => 0.5))
    @test m.name == "PKU-PAH"
    @test :Phe in m.species

    rm = RunManifest(
        code_git_sha="abc1234",
        data_hashes=Dict("ChEMBL" => "deadbeef"),
        rng_seed=UInt64(42),
        solver="Tsit5",
        abstol=1e-8,
        reltol=1e-6,
    )
    @test rm.rng_seed == 0x2a
    @test rm.solver == "Tsit5"
end
