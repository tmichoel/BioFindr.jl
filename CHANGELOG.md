# Changelog

All notable changes to this project will be documented in this file.

## [v1.2.0] - 2026-04-16

### Breaking Changes

There are no breaking changes to the public API in this release. However, users who accessed internal (unexported) functions or module-level constants by name should note:

- The internal module-level string constants `corr`, `link`, `med`, `relev`, and `pleio` (accessible as e.g. `BioFindr.corr`) have been removed
- Internal functions have been renamed to follow Julia's snake_case convention: `groupmeans` → `group_means`, `realLLR_col` → `real_llr_col`, `llrstats_col` → `llr_stats_col`
- Source files renamed for consistency: `realLLR.jl` → `real_llr.jl`, `randomLLR.jl` → `random_llr.jl` (plus corresponding test and documentation files)

### New Exports

- `LBeta`, `coerce_scitypes!`, and `generate_test_data` are now exported (they were previously accessible but unexported despite being used prominently in the documentation)

### Bug Fixes

- **`lbeta.jl`**: `logpdf` incorrectly returned `0.` for `x < 0`; corrected to `-Inf` (the log-probability of an impossible event) ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`posteriorprobs.jl`**: All five `try/catch` blocks previously caught every exception, silently masking real bugs. They now only catch `AssertionError`; all other exceptions are rethrown ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`findr.jl`**: `findr(dX, dG)` with a categorical genotype column `dG` silently produced a result matrix of `CategoricalValue` objects instead of integers. The function now validates the scitype of `dG` and uses `levelcode()` for correct integer conversion ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`utils.jl`**: `getpairs` previously caused a `MethodError` when column names in the two input data frames did not match; it now throws an informative `ErrorException` ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`utils.jl`**: `symprobs` docstring incorrectly claimed the default combination method was `"prod"`; corrected to `"none"` ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`utils.jl`**: `qvalue` now clamps q-values to the interval `[0, 1]` and replaces noisy `@info` logging with `@debug` ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))

### Performance Improvements

- **`utils.jl`** (`group_means`): `findall` is now computed once per group instead of evaluating a membership mask three times ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`supernormalization.jl`** (`supernormalize`): The quantile lookup table is now built once per column length instead of issuing one `quantile()` call per element ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- **`posteriorprobs.jl`** (`pi0est`): The λ-grid search now uses O(log n) `searchsortedlast` instead of an O(n) linear scan ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))

### Dependency Changes

- Add `CategoricalArrays` as an explicit runtime dependency (it was previously loaded as a hidden transitive dependency via DataFrames) ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- Move `Printf` from `[deps]` to `[extras]` (test-only) ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- Remove `Documenter` and `LiveServer` from runtime `[deps]` — these are documentation-only tools ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- Bump compat for `MetaGraphsNext` to include version 0.8 ([#23](https://github.com/tmichoel/BioFindr.jl/pull/23))
- Remove `Manifest.toml` and `docs/Manifest.toml` from the repository; add `Manifest.toml` to `.gitignore` ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))

### New Tests

Total test count increased from 145 to 221 ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25)):

- `test/findr_tests.jl` — comprehensive tests for all `findr()` and `findr_matrix()` overloads (the primary public API previously had zero tests)
- `test/dagfindr_tests.jl` — tests for all three `dagfindr!` methods, including cycle-freeness assertions
- `test/utils_tests.jl` — expanded with tests for `qvalue`, `globalfdr`, `globalfdr!`, `stackprobs`, `symprobs`, `combineprobs`, and `getpairs`

### Documentation

- Complete previously unfinished ("TBW") docstrings in `bayesiannets.jl` ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- Restore `push!(LOAD_PATH, "../src/")` in `docs/make.jl` to ensure the local source is used during documentation builds rather than the registered package version ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))

### CI

- Upgrade all GitHub Actions to their latest versions ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))
- Test against Julia `1` (latest stable) instead of a pinned minor version ([#25](https://github.com/tmichoel/BioFindr.jl/pull/25))

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.1.0...v1.2.0

## [v1.1.0] - 2026-04-09

### Breaking Changes

There are no breaking changes in this release.

### Changes

- Broaden function type signatures from concrete to abstract array/float types for improved interoperability (e.g., with [JuliaCall](https://juliapy.github.io/PythonCall.jl/stable/juliacall/)): `Array{T}` → `AbstractArray{T}`, `Vector{T}` → `AbstractVector{T}`, `Matrix{T}` → `AbstractMatrix{T}` across `realLLR.jl`, `posteriorprobs.jl`, `utils.jl`, `findr_matrix.jl`, and `findr_pvalues.jl` ([#22](https://github.com/tmichoel/BioFindr.jl/pull/22))
- Fix inconsistency in `findr_matrix` overloads where implementations used `T<:Real` while docstrings declared `T<:AbstractFloat`; all overloads now consistently use `T<:AbstractFloat` ([#22](https://github.com/tmichoel/BioFindr.jl/pull/22))
- Update dependency manifests
- Update Julia version in CI workflow

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.5...v1.1.0

## [v1.0.5] - 2025-05-21

- Catch error in kernel density estimation when π₀=1 by returning zero posterior probabilities
- CompatHelper: add new compat entry for ScientificTypes at version 3 (keep existing compat) ([#20](https://github.com/tmichoel/BioFindr.jl/pull/20))

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.4...v1.0.5

## [v1.0.4] - 2024-06-08

- Add option to perform causal inference for a subset of regulators without having to modify the input data

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.3...v1.0.4

## [v1.0.3] - 2024-04-27

- Added `dagfindr!` to exported functions

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.2...v1.0.3

## [v1.0.2] - 2024-04-26

- Added functions and documentation for DAG reconstruction

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.1...v1.0.2

## [v1.0.1] - 2024-04-26

- CompatHelper: add new compat entry for Graphs at version 1, keep existing compat ([#18](https://github.com/tmichoel/BioFindr.jl/pull/18))
- CompatHelper: add new compat entry for MetaGraphsNext at version 0.7, keep existing compat ([#19](https://github.com/tmichoel/BioFindr.jl/pull/19))

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v1.0.0...v1.0.1

## [v1.0.0] - 2024-02-06

- The package and repository have been renamed from **Findr.jl** to **BioFindr.jl** to comply with Julia package naming guidelines (see [JuliaRegistries/General#100261](https://github.com/JuliaRegistries/General/pull/100261))

**Full Changelog**: https://github.com/tmichoel/BioFindr.jl/compare/v0.2.0...v1.0.0

## [v0.2.0] - 2024-02-05

- Fixed bug in Test 3 (mediation test) where alternative hypothesis was selected instead of null
- Added new function to generate synthetic test data sampled from simple linear models
- Added unit tests for most functions

**Full Changelog**: https://github.com/tmichoel/Findr.jl/compare/v0.1.0...v0.2.0

## [v0.1.0] - 2023-08-18

- Initial release
- Most of the functionality of the [original findr software](https://github.com/lingfeiwang/findr) implemented
- See the [tutorials](https://tmichoel.github.io/FindrTutorials/) for usage examples and the [docs](https://tmichoel.github.io/Findr.jl/dev/) for more details

[v1.2.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.2.0
[v1.1.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.1.0
[v1.0.5]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.5
[v1.0.4]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.4
[v1.0.3]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.3
[v1.0.2]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.2
[v1.0.1]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.1
[v1.0.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.0
[v0.2.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v0.2.0
[v0.1.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v0.1.0
