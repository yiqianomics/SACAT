# SparseDOSSA2 n = 50 null calibration

This module evaluates DASRA 0.4.1, ZINQ 2.0, and MaAsLin3 1.5.3 with 50 control and 50 case samples under three no-signal SparseDOSSA2-based settings. It is a small-sample extension of the frozen SparseDOSSA2 robustness simulation and reuses its Stool template, 50-taxon panel, generator calibration, seeds, and formal method wrappers.

## Design

Each of 100 independent replications contains the same three method inputs and three null settings:

- `NULL_BALANCED`: median library depths 4,000 and 4,000;
- `NULL_BALANCED_LOW`: median library depths 1,500 and 1,500;
- `NULL_FOURFOLD`: median library depths 1,500 and 6,000.

There are no structural or abundance spikes. The balanced settings are exact global-null calibration experiments. The fourfold setting is a no-direct-effect depth-imbalance stress test; because sequencing depth changes marginal detection, it is not interpreted as a common marginal observed-prevalence null for every method.

The reported components are kept separate: DASRA structural absence and present-conditional abundance, ZINQ observed prevalence and detected-quantile abundance, and MaAsLin3 observed prevalence and detected-log abundance. Common BH and BY adjustments are recalculated within each fixed 50-taxon method-component family.

Primary raw rejection rates retain all 50 taxa in the denominator and count unavailable tests as non-rejections. Availability and the raw rejection rate conditional on availability are reported separately. Under the exact global null, the mean BH false-discovery proportion equals the probability of at least one BH rejection, so the summary labels it BH FWER.

## Reproduce

Run from the package root:

```sh
Rscript Analysis/SparseDOSSA2Simulation/n50_null/simulation.R run \
  --workers 7 --replications 100

Rscript Analysis/SparseDOSSA2Simulation/n50_null/summarize_results.R \
  Analysis/SparseDOSSA2Simulation/n50_null \
  Analysis/SparseDOSSA2Simulation/n50_null/summary/null_audit
```

The raw replication records are local generated artifacts and are ignored by Git. The frozen design and package provenance are in `design/`. Reviewer-facing results are in `summary/null_audit/`:

- `null_calibration.csv`: raw rejection, availability, available-conditional rejection, and BH/BY family error;
- `pvalue_tail.csv`: empirical p-value tails from 0.001 through 0.10;
- `positive_support.csv`: realized positive-count support without support-based filtering;
- `paired_differences.csv`: replication-paired method differences within the same broad domain;
- `audit_checks.csv`: the complete design, truth, multiplicity, and output validation record.

The analysis contract is `e567863b5a19204e64eeb0a008aca42d55b799be0f816bce59c7e4a877601b7b`. The parent runner and frozen method-reference SHA-256 values are recorded and checked by `simulation.R`.
