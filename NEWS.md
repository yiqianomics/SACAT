# SACAT 0.6.0

- `plot.sacat()` provides a publication-oriented dual-component association
  profile. Structural-absence and present-conditional abundance evidence share
  an aligned component Z-statistic axis, while covariate-standardized group
  summaries retain the fitted reference and comparison labels.
- `store_plot_data = TRUE` stores the compact descriptive summaries used by
  the profile, and structural companion summaries reuse the configured worker
  pool. Higher-order quadrature comparisons are controlled by `full_output`.
- Feature labels use only the space required by the displayed names. Component
  adjusted p-values use a shared adaptive `-log10(p)` color scale from one to a
  rounded lower endpoint, with an explicit numeric range available when plots
  need a fixed cross-study scale.
- The default plot separates group and evidence color semantics, uses a
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

# SACAT 0.5.2

- Expanded the README and reference manual with guidance for interpreting the
  present-conditional, reference-centered abundance result.

# SACAT 0.5.1

- The structural fitting control is named
  `structural_conditional_present_starts`.
- `abundance_quadrature_points` makes the relative-abundance quadrature order
  explicit. Detailed output can include a higher-order quadrature comparison.

# SACAT 0.5.0

- The relative-abundance arm uses one zero-truncated conditional-mark
  likelihood for every taxon. Conditioning on a positive count removes the
  structural-presence probability from this likelihood, including at zero
  structural-absence probability. The estimand is a present-conditional,
  target-excluded reference contrast.
- Relative-abundance standard errors use centered, sample-aligned
  influence contributions derived from the same conditional likelihood. The
  primary p-value remains a first-order normal Wald test.
- All-positive taxa can enter the abundance arm when their positive-mark model
  is identifiable.
- `min_positive_samples`, `min_reference_taxa`, and
  `structural_quadrature_points` control the analysis and numerical settings.
  Structural start selection accepts `"adaptive"` and `"full"`; numeric `1L`
  and `5L` inputs are also accepted.
- Optional ordered PSOCK workers and arm-level progress messages are available
  through `workers` and `verbose`.

# SACAT 0.4.2

- Updated the structural test to use nuisance-orthogonalized sample
  contributions with an empirical sandwich variance.
- Added handling for structural fits at the zero
  structural-absence-probability boundary.
