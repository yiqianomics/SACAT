# DASRA 0.5.0

- The relative-abundance arm now uses one zero-truncated conditional-mark
  likelihood for every taxon. Conditioning on a positive count removes the
  structural-presence probability from this likelihood, including at zero
  structural-absence probability, while preserving the existing
  present-conditional estimand and target-excluded reference correction.
- Relative-abundance standard errors now come from centered, sample-aligned
  influence contributions derived from the same conditional likelihood. The
  primary p-value remains a first-order normal Wald test.
- All-positive taxa can enter the abundance arm when their positive-mark model
  is otherwise identifiable. The public result columns and component labels are
  unchanged.
- `min_positive_samples`, `min_reference_taxa`, and
  `structural_quadrature_points` make previously fixed analysis and numerical
  controls explicit. Structural start selection now accepts the descriptive
  modes `"adaptive"` and `"full"`, while legacy `1L` and `5L` inputs remain
  supported.
- Optional ordered PSOCK workers and arm-level progress messages are available
  through `workers` and `verbose`.
- Package help, missing-data behavior, numerical controls, and the separated
  reference assumptions are now documented. The unused `nleqslv` dependency
  has been removed, and the R reference/compiled quadrature implementations
  have distinct names.

# DASRA 0.4.2

- Structural regular-interior inference now uses the same
  nuisance-orthogonalized per-sample estimating-function contributions for
  both the test numerator and its empirical sandwich variance. At an exact
  nuisance root, the numerator remains algebraically identical to the
  restricted target score.
- A covariate-adjusted structural nuisance fit is reported as unavailable when
  it is strictly dominated by the exact zero
  structural-absence-probability limit in the same detection objective.
- The model, estimands, public function arguments, result columns, and
  relative-abundance implementation are unchanged.
