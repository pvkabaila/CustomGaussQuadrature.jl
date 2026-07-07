# CustomGaussQuadrature v3.0.6

CustomGaussQuadrature 3.0.6 is a patch release that preserves the package API while tightening package metadata, making manual comparison scripts easier to run reproducibly, and improving cleanup of temporary BigFloat precision changes.

## Highlights

- Restored BigFloat precision to 256 bits on both normal and error exits from the precision-selection routines.
- Removed manual-script-only standard libraries from the package runtime dependency metadata.
- Added a separate `test/manual` environment for exploratory comparison scripts that use packages such as `QuadGK`, `PolyChaos`, and `Plots`.
- Refreshed the package metadata so the manifest is consistent with the current `Project.toml`.

## Compatibility Notes

- The package API is unchanged.
- `Pkg.test()` continues to use the package test target and does not require the manual comparison-script environment.
- Manual comparison scripts can be run after activating and instantiating the `test/manual` environment.
- The package version has been bumped from 3.0.5 to 3.0.6 in `Project.toml`.

## Validation

- The package test suite was run with `Pkg.test()`.
- All 52 tests passed.
