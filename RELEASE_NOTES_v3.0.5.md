# CustomGaussQuadrature v3.0.5

CustomGaussQuadrature 3.0.5 is a patch release that preserves existing behavior while refining the Stieltjes include-script interface and aligning the examples and manual with the current caller-visible names.

## Highlights

- Added primary Stieltjes output names `nodes_stieltjes`, `weights_stieltjes`, `a_vec_stieltjes`, `b_vec_stieltjes`, and `nbits_stieltjes` in both stored and user-defined include scripts.
- Preserved backward compatibility by continuing to provide the older `stieltjes_*` output names as aliases.
- Updated README and User Manual examples to use the new primary Stieltjes names.
- Expanded and normalized the manual examples so that moment-determinants and Stieltjes comparisons are clearer and use consistent example-specific names.
- Kept the R JuliaConnectoR helper interface method-neutral while clarifying this design in comments.

## Compatibility Notes

- Existing code that reads `stieltjes_nodes`, `stieltjes_weights`, `stieltjes_a_vec`, `stieltjes_b_vec`, or `stieltjes_nbits` after `include(...)` continues to work.
- The scalar `r` output remains unchanged.
- The package version has been bumped from 3.0.0 to 3.0.5 in `Project.toml`.

## Validation

- The stored Stieltjes workflow was checked after the rename-and-alias changes.
- The user-defined Stieltjes workflow was checked through the R JuliaConnectoR examples.
- Repository documentation and examples were updated to match the current interface.