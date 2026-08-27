# SparseDOSSA2 robustness simulation

This module evaluates DASRA and the comparison methods from the formal simulation with a SparseDOSSA2-based, reservoir-preserving semi-synthetic DGP initialized from the correlated Stool template. The design is fixed before any method is run, and every completed replication retains the generated counts, metadata, truth, method output, and generator checks needed for review.

## Frozen design

- SparseDOSSA2 0.99.2 at GitHub SHA `26a998a6e3a5f04d6a86cce14d6d3229ca82633e`, using the pretrained `Stool` template with `new_features = FALSE`.
- 120 independent samples per group and 50 tested taxa.
- Ten settings: balanced and fourfold-depth no-direct-effect settings; structural-only, abundance-only, and joint alternatives under both depth designs; and balanced 40% and 60% same-direction abundance stress settings.
- Balanced median depths are 4,000/4,000. Fourfold median depths are 1,500/6,000 for control/case. All depths use log-scale SD 0.45 and are bounded at 300 and 30,000.
- Structural effects use a SparseDOSSA prevalence-spike coefficient of magnitude `log(2)`. A signed coefficient `sβ` implies a structural logit-absence contrast of `-sβ`. Abundance effects use a native nonzero-absolute-abundance log fold change of `log(2)`.
- The first 10 signals in each replication alternate direction. The 20- and 30-signal stress settings use the nested first 20 and 30 taxa in one direction.

Signal taxa are sampled uniformly without replacement from the frozen panel with a replication-specific seed. The 10/20/30 sets are nested. Taxa are never replaced according to observed support, convergence, p-values, or any method result.

## Generator calibration

The panel is selected with a generator-only calibration of 20,000 Stool samples at depth 1,500 using seed 2026082301. Features must have template structural-zero probability in `[0.05, 0.90]` and expected detection probability at least 0.22. Fifty taxa are selected from the 53 eligible features by decreasing unconditional expected absolute abundance, with source feature name as the stable tie-breaker. The selected minimum expected detection probability is approximately 0.225, corresponding to at least 27 expected positives among 120 baseline samples.

The abundance magnitude is also fixed by generator-only oracle mapping. Positive and negative candidate effects are assessed separately through the final reservoir construction, and the smallest candidate whose median absolute conditional log-relative effect lies in `[0.30, 0.60]` in both directions is selected. `log(1.5)` fails this gate and `log(2)` passes in both directions. The complete feature-level and candidate-level mappings are written to the design directory. `truth_oracle.csv` additionally verifies, without method output, nonzero conditional-relative and expected-detection contrasts for abundance signals and nonzero expected-detection contrasts for joint signals in both directions. These mapped effects are calibration summaries, not claims that a native absolute-abundance coefficient equals a relative-abundance coefficient.

Negative structural effects intentionally increase absence. In structural and joint settings, a small number of directly affected taxa can therefore have fewer than 20 expected positives in a group. This signal-induced support loss is retained as part of the robustness assessment; availability is reported separately. The fourfold no-direct-effect setting is interpreted as a conditional-independence robustness check under strong depth imbalance, rather than as a claim that every working logistic coefficient is exactly zero.

## Reservoir-preserving construction

Let `R0` be the normalized null SparseDOSSA2 community and `RU` the total probability of all untested features. Each active tested taxon receives a private reservoir block with `c = 1/60`. For active taxon `j`,

```
B_j  = R0_j + c RU
W_j  = spiked absolute abundance of j
WO_j = c × total untested null absolute abundance
R_j  = B_j W_j / (W_j + WO_j)
```

The unused part `B_j - R_j` returns to `Other_unmodeled`. Inactive tested probabilities remain element-for-element equal to `R0`; only an active taxon and its private reservoir block exchange mass. The final probabilities are checked for nonnegativity and unit column sums before multinomial sequencing.

## Method parity

[`frozen/formal_hpc_method_reference.R`](frozen/formal_hpc_method_reference.R) is an exact copy of the final DASRA 0.6.0 formal simulation script with SHA-256 `7597ac45e2540043f1a2c56a824768ecdf7363b30bf77c134660c61f05c62d97`. The runner verifies the hash, parses the script, confirms the final expression is `main()`, skips only that expression, and directly calls its `analyze_setting()` function.

This preserves the original DASRA, ZINQ, MaAsLin3, edgeR, DESeq2, ANCOM-BC2, LinDA, corncob, and metagenomeSeq calls, seeds, options, output extraction, and common BH/BY adjustment. Generic abundance methods are omitted only in the two structural-only settings, matching the frozen design. Common multiplicity is applied within the fixed 50-taxon method-by-component family. Package-provided adjusted values remain diagnostic fields.

## Run

Run commands from the package root so the optional project R library can be found. A different library can be supplied through `DASRA_SPARSEDOSSA_R_LIB`.

```sh
Rscript Analysis/SparseDOSSA2Simulation/simulation.R design

Rscript Analysis/SparseDOSSA2Simulation/simulation.R smoke

Rscript Analysis/SparseDOSSA2Simulation/simulation.R run \
  --workers 7 --replications 100

Rscript Analysis/SparseDOSSA2Simulation/simulation.R summarize \
  --replications 100

Rscript Analysis/SparseDOSSA2Simulation/summarize_results.R
```

The reported analysis uses the prespecified replication IDs 1 through 100. The raw directory contains exactly these 100 records, and `summary/completion.csv` records the corresponding 100-replication analysis. Each worker is limited to one numerical thread.

## Outputs

By default, the module directory is the formal analysis root. Smoke output uses an independent session-temporary directory. The completed 100-replication analysis contains the source files plus:

- `design/`: contract, settings, panel and oracle calibration tables, plus version, RemoteSha, and installed-content hashes for every required package;
- `results/raw/`: one validated RDS file for each of replications 1 through 100, including counts, metadata, direct truth, support counts, generator invariants, and all method results;
- `summary/`: replication-level metrics, Monte Carlo standard errors and intervals, availability, DASRA structural-result counts, generator invariants, positive-support summaries, and aggregate performance. `null_calibration.csv`, `prevalence_benchmark.csv`, `abundance_benchmark.csv`, and `dense_stress.csv` provide the concise reviewer-facing views; `audit_checks.csv` records their validation contract.

Incomplete work retains `checkpoints/` and `logs/` for diagnosis and exact resume. They are not part of the completed 100-replication analysis directory.

The primary rejection rates keep every tested taxon in the denominator and treat unavailable results as non-rejections. `available_conditional_raw_null_rate` is reported separately to diagnose numerical availability without replacing the primary estimand.

## Truth convention

Truth fields describe direct generator perturbations. Structural truth marks prevalence-spiked taxa, abundance truth marks abundance-spiked taxa, and observed-prevalence truth marks either kind of directly spiked taxon because both can alter sequencing detection. Inactive tested taxa are exact probability-level nulls under the reservoir construction. The 60% same-direction setting is an explicit violation of the strict null-majority reference condition and is interpreted as a stress test rather than a nominal operating point.

SparseDOSSA2 prevalence spikes directly define the observed-prevalence benchmark used for ZINQ and MaAsLin3. DASRA structural-absence results are reported separately because they target the underlying structural component rather than the same observed-prevalence estimand. The abundance components are likewise interpreted as method-specific responses to the same native abundance perturbation, not as a single common-estimand ranking.
