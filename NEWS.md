# DASRA 0.6.0

- `plot.dasra()` provides a publication-oriented dual-component association
  profile. Structural-absence and present-conditional abundance evidence share
  an aligned component Z-statistic axis, while covariate-standardized group
  summaries retain the fitted reference and comparison labels.
- `store_plot_data = TRUE` prepares the compact descriptive summaries needed
  by the profile. Primary estimates, standard errors, and p-values are
  identical with either setting, and the default `FALSE` retains lightweight
  fitted objects.
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
  device, and the `file` argument writes PDF, PNG, or SVG output directly.
- BH-adjusted fits now show separate, family-wide structural and abundance
  discovery gates on the component Z-statistic axis. The data-dependent
  boundaries use the complete fitted component families and are invariant to
  displayed feature selection. A finite gate appears when the component has at
  least one BH discovery.

# DASRA 0.5.2

- The README and reference manual now present the relative-abundance method as
  a tutorial workflow: define the present-conditional reference-centered
  estimand, review formation and support fields, and use detailed reference and
  numerical diagnostics when needed.
- This is a documentation-only release. Functions, arguments, defaults, output
  schemas, fitted procedures, and numerical results are unchanged.

# DASRA 0.5.1

- The structural fitting control is now named
  `structural_conditional_present_starts`, making its scope explicit. This is a
  clean public-interface rename; the fitted structural procedure, argument
  position, and defaults are unchanged. Named calls using the previous argument
  must be updated to the new name.
- `abundance_quadrature_points` makes the relative-abundance quadrature order
  explicit while preserving the validated 41-node default exactly. Detailed
  output can report a fixed-fit comparison with a strictly higher-order rule;
  this is a numerical sensitivity diagnostic, not an error bound and does not
  alter inference.
- The `utils` namespace used for parallel-worker version checks is now declared
  explicitly in `Imports`. The canonical `GPL (>= 3)` license declaration is
  unchanged and does not require a separate `LICENSE` file.

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
