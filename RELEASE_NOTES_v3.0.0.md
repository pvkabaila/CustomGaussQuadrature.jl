# CustomGaussQuadrature v3.0.0

CustomGaussQuadrature 3.0.0 is a major release with breaking interface changes and substantial revisions across the moment determinants and Stieltjes workflows.

## Highlights

- Revised the user-facing `which_f` specification so support endpoints and parameters can be kept in their input representation until the working arithmetic type is chosen, while weight-function names are handled separately in canonical form.
- Added typed materialization helpers such as `materialize_scalar_spec_fn`, `materialize_integer_spec_fn`, and `materialize_support_spec_fn` to preserve exact decimal input and avoid premature conversion.
- Refactored the Stieltjes include-based drivers for both built-in and user-defined weights.
- Expanded test coverage for string-based specifications, built-in weights, and user-defined Weibull examples.
- Updated README examples and package documentation to match the revised interfaces.

## Breaking Changes

- Code written against earlier `which_f` conventions should be reviewed and updated to the current three-part form: weight-function name, support interval, and parameter specification.
- For finite non-integer constants that must preserve their decimal representation, pass strings such as `"2.1"` instead of pre-rounded floating-point literals.
- User-defined moment and log-weight functions should materialize values from `which_f` inside the chosen arithmetic type `T` rather than assuming parameters were already converted.
- The user-defined Stieltjes workflow now requires `mu0` and no longer requires the caller to supply `mu1`.
- The Stieltjes include scripts now operate on `lnf_typed_fn(T, which_f, x)` and the input `which_f` separately, rather than expecting the caller to pre-build the final one-argument closure.

## Migration Notes

- If you pass decimal-valued parameters or support endpoints, prefer strings in `which_f` and convert them inside your callback with `materialize_scalar_spec_fn`.
- If you use built-in stored weights, keep `which_f[1]` in its canonical lowercase form.
- If you use the Stieltjes include scripts, define `mu0` and optionally `offset`, `j_max`, and `epsilon` before `include(...)`.
- Review user-defined Weibull or other custom examples to ensure parameters are read from `which_f[3]` inside the callback.

## Validation

- The revised package metadata has been bumped to version 3.0.0.
- Current repository tests cover string-spec support for both moment determinants and Stieltjes paths, along with the updated built-in and user-defined workflows.