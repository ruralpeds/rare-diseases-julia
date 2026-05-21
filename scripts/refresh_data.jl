# scripts/refresh_data.jl
#
# Iterate every registered source, fetch into data/raw/<source>/<ISO-date>/,
# compute sha256, and update data/manifest.toml. Run nightly via the
# data-refresh GitHub Action (Phase 4+).
#
# Usage:
#   julia --project=. scripts/refresh_data.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using RDDataSources

cache_root = joinpath(@__DIR__, "..", "data", "raw")
mkpath(cache_root)

for name in registered_sources()
    @info "Refreshing source" name
    # Each source type is registered with its name; the registry maps name
    # to its constructor. Skip ones whose fetch! is still a stub.
    try
        T = RDDataSources._REGISTRY[name]
        src = T()
        m = RDDataSources.manifest(src)
        @info "Manifest" name url=first(m.urls) license=m.license
        # Phase 4+: actually fetch and stamp manifest.toml here.
    catch err
        @warn "skip" name err
    end
end
