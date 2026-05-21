# scripts/demo.jl
#
# End-to-end demo of the platform against the bundled fixtures.
#
# Run with:
#   julia --project=. scripts/demo.jl
#
# No network access required.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
for p in readdir(joinpath(@__DIR__, "..", "packages"); join=false)
    Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..", "packages", p)))
end
Pkg.instantiate()

using RareDiseaseCore
using RDOntology
using RDDiagnostics
using RDPathways
using RDTreatment
using RDSimulation
using OrdinaryDiffEq
using Random

println("=" ^ 72)
println("rare-diseases-julia — end-to-end demo")
println("NOT FOR CLINICAL USE — research and education only")
println("=" ^ 72)

# ---------------------------------------------------------------------------
# 1. Ontology — load the mini HPO and compute IC from a tiny corpus
# ---------------------------------------------------------------------------
fixture = joinpath(@__DIR__, "..", "packages", "RDOntology",
                   "test", "fixtures", "mini_hpo.obo")
g = load_hpo(fixture)
println("\n[1] Loaded HPO graph: $(length(g)) terms.")

annotations = Dict(
    "DISEASE:SEIZ"   => ("Seizure-only disease",     ["HP:0001250"]),
    "DISEASE:SEIZID" => ("Seizure + ID disease",     ["HP:0001250", "HP:0001249"]),
    "DISEASE:VIS"    => ("Vision-only disease",      ["HP:0000505"]),
)
information_content!(g; annotations=Dict(d => phs for (d, (_, phs)) in annotations))
println("    IC computed from $(length(annotations)) synthetic diseases.")

# ---------------------------------------------------------------------------
# 2. Differential diagnosis for a seizure patient
# ---------------------------------------------------------------------------
case = PatientCase(phenotypes_present=[HPOId("HP:0001250")])
dx = rank_diagnoses(case, g; disease_annotations=annotations, topk=5)
println("\n[2] Diagnosis ranking for {HP:0001250 / Seizure}:")
for (i, c) in enumerate(dx.candidates)
    println("    $(i). $(c.name) ($(c.disease))  score=$(round(c.score; digits=3))")
end

# ---------------------------------------------------------------------------
# 3. Mechanistic simulation: PAH/PKU on diet vs sapropterin
# ---------------------------------------------------------------------------
println("\n[3] PKU steady-state Phe by variant class:")
for class in (:wildtype, :mild_pku, :classical_pku, :null)
    prob = pah_pku_problem(; variant=class, tspan=(0.0, 96.0))
    sol  = solve(prob, Tsit5(); abstol=1e-9, reltol=1e-8)
    println("    $(rpad(string(class), 18)) Phe = $(round(sol.u[end][1]; digits=3))")
end

println("\n    Sapropterin effect on classical PKU:")
prob_off = pah_pku_problem(; variant=:classical_pku, bh4=1.0, tspan=(0.0, 96.0))
prob_on  = pah_pku_problem(; variant=:classical_pku, bh4=2.0, tspan=(0.0, 96.0))
sol_off = solve(prob_off, Tsit5(); abstol=1e-9, reltol=1e-8)
sol_on  = solve(prob_on,  Tsit5(); abstol=1e-9, reltol=1e-8)
println("      Phe without BH4: $(round(sol_off.u[end][1]; digits=3))")
println("      Phe with BH4 2x: $(round(sol_on.u[end][1];  digits=3))")

# ---------------------------------------------------------------------------
# 4. Treatment ranking by network proximity
# ---------------------------------------------------------------------------
net = PathwayNetwork()
add_edge_undirected!(net, "SCN1A", "GABA_R")
add_edge_undirected!(net, "GABA_R", "GAD1")
add_edge_undirected!(net, "GAD1", "HMGCR")
add_edge_undirected!(net, "HMGCR", "LDLR")

drugs = [
    DrugRecord(name="lamotrigine",       targets=["SCN1A"],
               evidence_tier=APPROVED),
    DrugRecord(name="atorvastatin",      targets=["HMGCR"],
               evidence_tier=APPROVED),
    DrugRecord(name="experimental_gaba", targets=["GABA_R"],
               evidence_tier=REPURPOSING_HYPOTHESIS),
]
cands, warning = rank_treatments(["SCN1A"], drugs, net;
                                  n_bootstrap=30, rng=Xoshiro(0))
println("\n[4] Treatment ranking for disease module {SCN1A}:")
for (i, c) in enumerate(cands)
    println("    $(i). $(c.drug.name)  tier=$(c.drug.evidence_tier)  ",
            "d=$(round(c.proximity_d; digits=2))  ",
            "score=$(round(c.score; digits=3))")
end

println("\n" * "=" ^ 72)
println("Done. " * warning)
println("=" ^ 72)
