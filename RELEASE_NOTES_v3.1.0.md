# CustomGaussQuadrature v3.1.0

CustomGaussQuadrature 3.1.0 is a minor release that preserves the package API while tightening the default convergence tolerance used by the Stieltjes procedure and reducing the runtime of expensive Stieltjes tests.

## Highlights

- Reduced the default value of `epsilon` in `stieltjes_a_vec_b_vec_final_fn` from `1.0e-15` to `1.0e-17`.
- Kept the public package API unchanged.
- Reduced the size and cost of selected test examples.
- Removed an expensive Stieltjes gamma-pdf test whose main coverage was already provided by other tests.
- Retained coverage for 2-vector parameter specifications through the moment-determinants gamma test and the stored generalized-normal Stieltjes regression test.
- Updated the generalized-normal Stieltjes regression example to use faster nontrivial parameters.

## Compatibility Notes

- The package API is unchanged.
- Existing calls that pass `epsilon` explicitly are unaffected.
- Calls that rely on the default Stieltjes `epsilon` may use a stricter convergence criterion than in previous releases.
- The package version has been bumped from 3.0.6 to 3.1.0 in `Project.toml`.

## Validation

- The package test suite should be run with `Pkg.test()` after the release edits.
- The manifest should be refreshed or checked after the `Project.toml` version bump.
