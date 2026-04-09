# Changelog

All notable changes to this project will be documented in this file.

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

[v1.0.5]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.5
[v1.0.4]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.4
[v1.0.3]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.3
[v1.0.2]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.2
[v1.0.1]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.1
[v1.0.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v1.0.0
[v0.2.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v0.2.0
[v0.1.0]: https://github.com/tmichoel/BioFindr.jl/releases/tag/v0.1.0
