# Selection and screening of additional real-data analyses

The two original colorectal-cancer and *Clostridioides difficile* analyses
were already part of the project. Additional public datasets were considered
under the rules below. Twelve candidates completed the DASRA applicability
screen, and five of the passing candidates were carried through the full
multi-method workflow before the study scope was closed.

The readable candidate-level results are retained in
[`screening/dasra_candidate_screen.csv`](screening/dasra_candidate_screen.csv)
and [`screening/dasra_unavailability_reasons.csv`](screening/dasra_unavailability_reasons.csv).
The downloaded files for candidates that were not included in the final
analysis were removed from the project directory.

## Data and design eligibility

A candidate required:

- a citable public study and a public nonnegative integer count table before
  rarefaction or conversion to relative abundance;
- a recoverable original sample library size before analysis-level taxonomic
  aggregation and filtering;
- a documented binary comparison and a small, interpretable adjustment set;
- one observation per independent biological unit, or a group-blind rule for
  selecting one observation when repeated samples were available; and
- no overlap in cohort, study accession, or samples with the negative-control
  benchmark datasets.

At least 100 observations per group was preferred. Groups of 80--99 were
allowed as clearly labeled near-threshold exceptions when the other criteria
were satisfied. Zupancic obesity, the Ravel comparison, and the produce
farming candidate entered screening under this exception. No candidate passed
or failed the DASRA applicability screen because of the direction of an
estimated effect, a combined discovery count, comparison-method results, or
figure appearance.

Use of the same public repository did not constitute overlap when the
underlying studies and samples were distinct.

## DASRA applicability screen

Availability was the proportion of all retained analysis taxa for which a
requested DASRA result was formally available. A dataset was classified as
broadly non-working when availability was below 70% for structural absence,
present-conditional abundance, or the combined result.

A separate gross common-background flag was defined for a dataset in which at
least 50% of all tested taxa were BH discoveries for present-conditional
abundance. Unavailable results remained in that BH family with an operational
p-value of one. This was an operational check for an obviously unsuitable
common-background setting, not a formal proof of the model assumption.

Passing the screen means that the requested DASRA components were broadly
available on the prespecified panel. It does not imply that the dataset must
contain a statistically significant association. In particular, the Korean
hypertension and Zupancic obesity datasets passed and were retained even
though their formal combined analyses had no BH discoveries.

## Completed screening ledger

Twelve candidates completed screening. Eight passed and four failed because
at least one component had availability below 70%. None triggered the gross
common-background flag. The screening rejection fraction is therefore 4/12,
or 33.33%.

| Candidate | Screen result | Final disposition |
|---|---|---|
| Zupancic Old Order Amish obesity | Passed | Formal analysis completed |
| Ravel vaginal microbiome comparison | Passed | Formal analysis completed |
| Korean hypertension | Passed | Formal analysis completed |
| RISK pediatric Crohn's disease | Passed | Formal analysis completed |
| GEMS pediatric diarrhea | Passed | Formal analysis completed |
| Fresh-produce farming practice | Failed: component availability below 70% | Removed before formal analysis |
| TwinsUK urinary microbiome age comparison | Failed: component availability below 70% | Removed before formal analysis |
| Matched pediatric hand, foot, and mouth disease | Failed: component availability below 70% | Removed before formal analysis |
| Pre-infection cervical chlamydia comparison | Failed: component availability below 70% | Removed before formal analysis |
| Parkinson's disease | Passed | Formal workflow attempted, then excluded because three public profiles had only 1, 3, and 10 total reads and the comparison workflow did not complete reliably |
| iMSMS household-control adult sex | Passed | Not advanced after the formal study scope was closed |
| ISALA recent antibiotic exposure | Passed | Not advanced after the formal study scope was closed |

The two original analyses were not part of this 12-candidate denominator. The
final formal collection therefore contains seven datasets: the two original
analyses plus the five additional completed analyses. iMSMS and ISALA were not
advanced after scope closure and were not selected using formal comparison
outputs or figures. Parkinson's disease was handled separately as described in
the table: its screening result remains in the denominator, but its attempted
formal output was excluded and removed.

## Sources considered but not included in the completed screen

Several sources were rejected before screening because the public files or
study design did not meet the eligibility rules.

| Candidate | Reason it did not reach screening |
|---|---|
| American Gut mental-health comparison | The source collection overlapped an American Gut source already represented in the negative-control benchmark. |
| Vangay immigration cohort | The located public analysis tables were normalized or relative-abundance tables rather than nonrarefied integer counts. |
| Beijing rheumatoid-arthritis cohort | The public count object had been rarefied to exactly 10,000 reads per sample, so the original library sizes were not recoverable from that object. |
| HCHS/SOL migration contrast | The public metadata lacked the continuous years-in-the-United-States information needed to reconstruct the published age-at-relocation comparison. |
| Qiita 2229 seaweed condition | The selected public metadata and count source did not reproducibly reconstruct the intended published cohort. |
| MicroCOPD oral-wash comparison | The public microbiome object and assay metadata could not be reconciled into the prespecified independent oral-wash cohort without undocumented assumptions. |
| Qiita 2259 stickleback sex | The reconstructed fecal cohort contained 31 female and 18 male fish, below the minimum group size. |
| Qiita 1721 soil biochar | The public metadata did not establish enough independent field plots, as opposed to repeated cores from a smaller number of plots. |

A wider shortlist was prepared while the public-data search was in progress,
including 1000 Homes residence type, Goodrich obesity, HCHS/SOL adult sex,
Guangdong adult sex, MIBI adult sex, bat sex, METS geography, and wild baboon
sex. Their DASRA screens were not completed before the formal scope was
reduced, so they do not enter either the 12-candidate screening denominator or
the seven-dataset result summaries.

## How the denominators should be interpreted

- The candidate-screen rejection percentage uses all 12 candidates that
  completed DASRA screening, including rejected candidates and passing
  candidates that were not formally analyzed.
- The completed-analysis workability summary uses the seven datasets with full
  multi-method results.
- Signal-category percentages use combined DASRA discoveries from those seven
  datasets. Equal-dataset signal percentages are defined only for datasets
  with at least one such discovery.
- These are descriptive summaries of a deliberately assembled public-data
  collection. They are not estimates of prevalence across all microbiome
  studies and are not meta-analysis results.
