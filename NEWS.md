# DASRA 0.6.0

- `plot.dasra()` provides a publication-oriented dual-component association
  profile. Structural-absence and present-conditional abundance evidence share
  an aligned component Z-statistic axis, while covariate-standardized group
  summaries retain the fitted reference and comparison labels.
- `store_plot_data = TRUE` prepares the compact descriptive summaries needed
  by the profile. Structural companion summaries reuse the configured worker
  pool, and higher-order quadrature sensitivity checks remain tied to
  `full_output`. Primary estimates, standard errors, and p-values are identical
  with either setting; the default `FALSE` retains lightweight fitted objects.
- Feature labels use only the space required by the displayed names. Component
  adjusted p-values use a shared adaptive `-log10(p)` color scale from one to a
  rounded lower endpoint, with an explicit numeric range available when plots
  need a fixed cross-study scale.
- The default plot now separates group and evidence color semantics, uses a
  continuous component-adjusted-p legend with the 0.05 landmark when it lies
  within the resolved scale, and applies manuscript-scale typography, spacing,
  and adaptive output dimensions. Optional display-only group labels support
  publication figures with long fitted group names.
- The plot follows standard R conventions: `plot(fit)` draws to the active
  device, and the `file` argument writes PDF, PNG, or SVG output directly. SVG
  output is available in R builds with Cairo support.
- BH-adjusted fits display separate structural and abundance discovery
  thresholds based on the complete fitted component families.
- Parallel analyses use ordered, cross-platform worker processes.

# DASRA 0.5.2

- Expanded the README and reference manual with guidance for interpreting the
  present-conditional, reference-centered abundance result.

# DASRA 0.5.1

- The structural fitting control is now named
  `structural_conditional_present_starts`. Named calls using the previous
  argument must be updated to the new name.
- `abundance_quadrature_points` makes the relative-abundance quadrature order
  explicit. Detailed output can include a higher-order quadrature comparison.
- The `utils` namespace is now declared explicitly in `Imports`.

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

- Updated the structural test to use nuisance-orthogonalized sample
  contributions with an empirical sandwich variance.
- Added handling for structural fits at the zero
  structural-absence-probability boundary.
