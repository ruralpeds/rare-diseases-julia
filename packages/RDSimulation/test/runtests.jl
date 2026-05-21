using Test
using Dates
using Random
using RDSimulation

@testset "RDSimulation" begin

    @testset "RK4 on exponential decay matches analytic" begin
        # du/dt = -k*u with u(0)=1  -> u(t) = exp(-k*t)
        k = 0.5
        f!(du, u, _, _t) = (du[1] = -k * u[1]; nothing)
        t, U = RDSimulation.rk4(f!, [1.0], nothing, 0.0, 5.0, 0.01)
        for (i, ti) in enumerate(t)
            @test isapprox(U[1, i], exp(-k * ti); atol=1e-5)
        end
        @test t[end] == 5.0
    end

    @testset "PAH residual activity table" begin
        @test pah_residual_activity(:wildtype) == 1.0
        @test pah_residual_activity(:null) == 0.0
        @test pah_residual_activity(:classical_pku) < pah_residual_activity(:mild_pku)
        @test_throws ArgumentError pah_residual_activity(:bogus)
    end

    @testset "PKU model: variant + treatment monotonicity" begin
        # Wildtype reaches lower steady-state [Phe] than classical PKU,
        # and sapropterin should reduce [Phe] in residual-activity variants
        # vs untreated.
        u0 = [0.1, 0.1]   # mmol/L
        tspan = (0.0, 48.0)  # hours

        function steady_phe(class::Symbol; bh4::Float64=1.0)
            p = variant_effect(default_pah_parameters(), class; bh4_factor=bh4)
            m = MolecularModel(p, "PKU-$class-bh4$bh4")
            res = simulate(m;
                u0=u0,
                tspan=tspan,
                dt=0.05,
                manifest=RunManifest(
                    code_git_sha="test",
                    data_hashes=Dict{String,String}(),
                    rng_seed=UInt64(0),
                    solver="RK4",
                    abstol=1e-8,
                    reltol=1e-6,
                ),
            )
            return res.u[1, end]
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

    @testset "SimulationResult shape and citations" begin
        p = variant_effect(default_pah_parameters(), :mild_pku)
        m = MolecularModel(p)
        res = simulate(m;
            u0=[0.1, 0.1],
            tspan=(0.0, 1.0),
            dt=0.1,
            citations=["Blau N. Lancet 2010. PMID:20971365"],
            manifest=RunManifest(
                code_git_sha="abc",
                data_hashes=Dict("BioModels" => "deadbeef"),
                rng_seed=UInt64(1),
                solver="RK4",
                abstol=1e-8, reltol=1e-6,
            ),
        )
        @test size(res.u, 1) == 2
        @test res.species == [:Phe, :Tyr]
        @test "Blau N. Lancet 2010. PMID:20971365" in res.citations
        @test res.manifest.solver == "RK4"
    end

    @testset "Cohort simulate" begin
        rng = Xoshiro(0)
        agents = [CohortAgent(i, :on_diet) for i in 1:50]
        # Constant per-tick transition: small probability of crisis.
        function step!(a, _t, _dt, rng)
            if a.state == :on_diet && rand(rng) < 0.02
                a.state = :metabolic_crisis
            end
            return nothing
        end
        cohort_simulate(agents, step!; tspan=(0.0, 100.0), dt=1.0, rng=rng)
        ncrisis = count(a -> a.state == :metabolic_crisis, agents)
        @test 0 < ncrisis < 50      # some, but not all, drifted to crisis
        # Every agent's history covers the full timespan
        @test all(a -> last(a.history)[1] == 100.0, agents)
    end
end
