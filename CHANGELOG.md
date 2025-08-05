# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2025-08-05

### Changed

- Bumped compat for *DataStructures.jl* from `0.18.15` to `0.18.15 - 0.19` (#39).
- Changed `check_spectrum_integrality` to compute the real integer eigenvalues lazily (comparison was already lazy, but taking the real part and rounding was not) (#37).

### Removed

- Removed the unnecessary and redundant `_laplacian_1neg_spectra(spec::SSpectra)` method from `src/laplacian_s_spectra.jl` (this method literally did nothing whatsoever; I left it in from a previous design approach) (#38).

## [0.1.0] - 2025-08-01

### Added

- Released the initial stable version of the package.

[unreleased]: https://github.com/GraphQuantum/SDiagonalizability.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/Luis-Varona/MatrixBandwidth.jl/releases/tag/v0.1.0
