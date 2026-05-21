.PHONY: help dev test test-core test-onto test-sim integration fmt clean

JULIA ?= julia
PROJECT = --project=.

help:
	@echo "Targets:"
	@echo "  dev          Pkg.develop every in-tree package into the meta env"
	@echo "  test         Run tests across every in-tree package"
	@echo "  test-core    Run RareDiseaseCore tests only"
	@echo "  test-onto    Run RDOntology tests only"
	@echo "  test-sim     Run RDSimulation tests only"
	@echo "  integration  Co-load every package (test/runtests.jl)"
	@echo "  fmt          Run JuliaFormatter on the tree"
	@echo "  clean        Remove Manifest.toml and built artifacts"

dev:
	$(JULIA) $(PROJECT) scripts/dev_bootstrap.jl

test:
	@for pkg in packages/*/; do \
		echo "==> $$pkg"; \
		$(JULIA) --project=$$pkg -e 'using Pkg; \
			for p in readdir("packages"; join=false); \
				Pkg.develop(PackageSpec(path=joinpath("packages", p))); \
			end; \
			Pkg.instantiate(); Pkg.test()' || exit 1; \
	done

test-core:
	$(JULIA) --project=packages/RareDiseaseCore -e 'using Pkg; Pkg.test()'

test-onto:
	$(JULIA) --project=packages/RDOntology -e 'using Pkg; \
		Pkg.develop(PackageSpec(path="packages/RareDiseaseCore")); \
		Pkg.test()'

test-sim:
	$(JULIA) --project=packages/RDSimulation -e 'using Pkg; \
		Pkg.develop(PackageSpec(path="packages/RareDiseaseCore")); \
		Pkg.test()'

integration: dev
	$(JULIA) $(PROJECT) -e 'include("test/runtests.jl")'

fmt:
	$(JULIA) -e 'using Pkg; Pkg.add("JuliaFormatter"); \
		using JuliaFormatter; format(".", verbose=true)'

clean:
	find . -name Manifest.toml -not -path "./.git/*" -delete
	rm -rf docs/build docs/site
	rm -rf data/raw data/interim data/processed
	mkdir -p data/raw data/interim data/processed
	touch data/raw/.gitkeep data/interim/.gitkeep data/processed/.gitkeep
