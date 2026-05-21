using Test
using Random
using RDPathways

@testset "RDPathways" begin

    @testset "Empty network" begin
        n = PathwayNetwork()
        @test length(n) == 0
        @test !has_node(n, "X")
        @test !has_edge(n, "X", "Y")
        @test isempty(neighborhood(n, "X"))
    end

    @testset "Add nodes, edges, pathways" begin
        n = PathwayNetwork()
        add_node!(n, "PAH"; pathways=["R-HSA-71291"])
        add_node!(n, "TYR")
        add_edge_undirected!(n, "PAH", "TYR"; weight=0.7)
        @test has_node(n, "PAH") && has_node(n, "TYR")
        @test has_edge(n, "PAH", "TYR")
        @test has_edge(n, "TYR", "PAH")
        @test "R-HSA-71291" in pathways_for_gene(n, "PAH")
        @test sort(neighbors_of(n, "PAH")) == ["TYR"]
    end

    @testset "Directed and signed edges" begin
        n = PathwayNetwork()
        add_edge_directed!(n, "A", "B"; weight=0.9, sign=+1)
        add_edge_directed!(n, "B", "C"; weight=0.5, sign=-1)
        @test has_edge(n, "A", "B")
        @test !has_edge(n, "B", "A")
        @test n.signs[("A", "B")] == +1
        @test n.signs[("B", "C")] == -1
    end

    @testset "Shortest path on a small grid" begin
        n = PathwayNetwork()
        # Path: A - B - C - D, plus shortcut A - D' - D
        add_edge_undirected!(n, "A", "B")
        add_edge_undirected!(n, "B", "C")
        add_edge_undirected!(n, "C", "D")
        add_edge_undirected!(n, "A", "Dp")
        add_edge_undirected!(n, "Dp", "D")

        @test shortest_path_length(n, "A", "A") == 0
        @test shortest_path_length(n, "A", "B") == 1
        @test shortest_path_length(n, "A", "D") == 2     # via shortcut
        @test shortest_path_length(n, "A", "missing") == typemax(Int)

        p = shortest_path(n, "A", "D")
        @test length(p) == 3
        @test first(p) == "A" && last(p) == "D"
    end

    @testset "Neighborhood radius" begin
        n = PathwayNetwork()
        for (a, b) in [("A","B"),("B","C"),("C","D"),("D","E")]
            add_edge_undirected!(n, a, b)
        end
        @test neighborhood(n, "A"; radius=1) == Set(["B"])
        @test neighborhood(n, "A"; radius=2) == Set(["B","C"])
        @test neighborhood(n, "A"; radius=4) == Set(["B","C","D","E"])
    end

    @testset "Guney 2016 closest_distance" begin
        n = PathwayNetwork()
        for (a, b) in [("S1","X"),("X","T1"),("X","Y"),("Y","T2"),("S2","Y")]
            add_edge_undirected!(n, a, b)
        end
        # S1 -> T1 length 2, S1 -> T2 length 3 -> min 2
        # S2 -> T2 length 1, S2 -> T1 length 3 -> min 1
        d = closest_distance(n, ["S1","S2"], ["T1","T2"])
        @test d == (2 + 1) / 2

        # Empty source or target -> Inf
        @test isinf(closest_distance(n, String[], ["T1"]))
        @test isinf(closest_distance(n, ["S1"], String[]))

        # Unreachable source -> Inf
        add_node!(n, "ISOLATED")
        @test isinf(closest_distance(n, ["ISOLATED"], ["T1"]))
    end

    @testset "network_proximity_z" begin
        n = PathwayNetwork()
        for (a, b) in [("S","X"),("X","T")]
            add_edge_undirected!(n, a, b)
        end
        # Add filler nodes so the bootstrap has options
        for i in 1:20
            add_edge_undirected!(n, "F$i", "T")
        end
        rng = Xoshiro(0)
        r = network_proximity_z(n, ["S"], ["T"]; n_bootstrap=50, rng=rng)
        @test isfinite(r.d)
        @test isfinite(r.μ)
        @test isfinite(r.σ) || r.σ == 0.0
    end

    @testset "Loader stubs error clearly" begin
        @test_throws ErrorException load_reactome("/none")
    end
end
