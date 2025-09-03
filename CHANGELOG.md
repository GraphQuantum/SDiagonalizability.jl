# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Fixed a typo in some comments in the *Aqua.jl* static analysis tests (#60).

## [0.1.2] - 2025-08-13

### Added

- Added **References** sections to docstrings for immediate readability in the REPL and in the source code without needing to open the Documenter-generated website (#47, #48).

### Changed

- Finished the docstring for `_assert_graph_has_defined_s_bandwidth` function (#55).
- Made the return type of the `_pot_kernel_1neg_eigvecs` and `_pot_nonkernel_1neg_eigvecs` functions consistent regardless of the `n` parameter passed (#51).
- Finished the docstrings for `src/eigenvector_generation.jl`, fixing some minor inaccuracies along the way (#51).

### Fixed

- Fixed an `UndefVarError` in the `_assert_graph_has_defined_s_bandwidth` function (#55).

## [0.1.1] - 2025-08-05

### Changed

- Bumped compat for *DataStructures.jl* from `0.18.15` to `0.18.15 - 0.19` (#39).
- Changed `check_spectrum_integrality` to compute the real integer eigenvalues lazily (comparison was already lazy, but taking the real part and rounding was not) (#37).

### Removed

- Removed the unnecessary and redundant `_laplacian_1neg_spectra(spec::SSpectra)` method from `src/laplacian_s_spectra.jl` (this method literally did nothing whatsoever; I left it in from a previous design approach) (#38).

## [0.1.0] - 2025-08-01

### Added

- Released the initial stable version of the package.

[unreleased]: https://github.com/GraphQuantum/SDiagonalizability.jl/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.2
[0.1.1]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
