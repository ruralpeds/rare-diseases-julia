# Contributing

## Ground rules

1. **No new data source without a row in `data/SOURCES.tsv`.** CI will block
   PRs that add a downloader or parser referencing an unregistered source.
2. **No closed data.** If a source's license forbids redistribution, the
   downloader must fetch on the user's machine and never write into a
   committed artifact.
3. **Every public function has a docstring** with at least one runnable
   example.
4. **Tests required.** Unit tests in the package's `test/`; cross-package
   flows in top-level `test/`.

## Package layout

Each `packages/<Name>/` is its own Julia package with `Project.toml`, `src/`,
`test/`. Add to the dev environment with:

```julia
julia --project=. -e 'using Pkg; Pkg.develop(path="packages/<Name>")'
```

## Style

- `JuliaFormatter` with the repo `.JuliaFormatter.toml` (default style).
- No `using X: *` in package code; prefer explicit imports.
- Identifier types from `RareDiseaseCore` (`HPOId`, `MondoId`, …) instead of
  raw `String` for any ID that has a registered namespace.

## Commit messages

- Imperative mood. One-line subject ≤ 72 chars, blank line, body if needed.
- Reference issue numbers in the body, not the subject.

## Release

Stable packages are released via the Julia General registry (TagBot). The
meta-repo gets a coordinated tag at the end of each milestone.
