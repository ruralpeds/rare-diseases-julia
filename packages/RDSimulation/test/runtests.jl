using Test
using Dates
using Random
using OrdinaryDiffEq
using Agents
using RDSimulation

@testset "RDSimulation" begin

    @testset "PAH residual activity table" begin
        @test pah_residual_activity(:wildtype) == 1.0
        @test pah_residual_activity(:null) == 0.0
        @test pah_residual_activity(:classical_pku) < pah_residual_activity(:mild_pku)
        @test_throws ArgumentError pah_residual_activity(:bogus)
    end

    @testset "PKU model: variant + treatment monotonicity via Catalyst+OrdinaryDiffEq" begin
        function steady_phe(class::Symbol; bh4::Float64=1.0)
            prob = pah_pku_problem(; variant=class, bh4=bh4, tspan=(0.0, 96.0))
            sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
            # Phe is the first species in PAH_PKU.
            return sol.u[end][1]
        end

        wt    = steady_phe(:wildtype)
        mild  = steady_phe(:mild_pku)
        class = steady_phe(:classical_pku)
        nul   = steady_phe(:null)

        @test wt < mild < class < nul

        # Sapropterin (bh4=2.0) lowers steady-state Phe in residual-activity
        # variants, but not in null (no enzyme to potentiate).
        @test steady_phe(:mild_pku; bh4=2.0)      < steady_phe(:mild_pku)
        @test steady_phe(:classical_pku; bh4=2.0) < steady_phe(:classical_pku)
        @test steady_phe(:null; bh4=2.0)          ≈ steady_phe(:null) atol=1e-9
    end

    @testset "Sapropterin PBPK — peak after absorption then decay" begin
        prob = sapropterin_pbpk_problem(; dose_mg=10.0, tspan=(0.0, 24.0))
        sol = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
        # A_central starts at 0, must rise above 0 then decline by t=24h
        central = [u[2] for u in sol.u]
        peak, peak_idx = findmax(central)
        @test peak > 0
        @test peak_idx > 1
        @test central[end] < peak
        # Mass balance: gut + central never exceeds the initial dose
        for u in sol.u
            @test u[1] + u[2] ≤ 10.0 + 1e-6
        end
    end

    @testset "RunManifest stamping" begin
        rm = RunManifest(
            code_git_sha="abc",
            data_hashes=Dict("BioModels" => "deadbeef"),
            rng_seed=UInt64(1),
            solver="Tsit5",
            abstol=1e-8, reltol=1e-6,
        )
        @test rm.solver == "Tsit5"
        @test rm.rng_seed == 0x01
    end

    @testset "Cohort via Agents.jl" begin
        function step!(agent, model)
            if agent.state == :on_diet && rand(abmrng(model)) < 0.02
                agent.state = :metabolic_crisis
            end
            return nothing
        end

        model = build_cohort_model(;
            n_agents=50,
            initial_state=:on_diet,
            agent_step! = step!,
            rng=Xoshiro(0),
        )
        run_cohort!(model, 100)
        agents = collect(allagents(model))
        @test length(agents) == 50
        ncrisis = count(a -> a.state == :metabolic_crisis, agents)
        @test 0 < ncrisis < 50
    end
end
