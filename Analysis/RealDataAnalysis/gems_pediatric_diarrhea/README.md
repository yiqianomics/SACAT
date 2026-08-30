# GEMS pediatric diarrhea

This analysis uses the public genus-count data from the Global Enteric
Multicenter Study analysis by Pop et al., *Genome Biology* 2014
([doi:10.1186/gb-2014-15-6-r76](https://doi.org/10.1186/gb-2014-15-6-r76)).
The data file is publicly available from the authors'
[MSD1000 repository](https://github.com/jnpaulson/MSD1000) as
`forserveroptim.rdata`.

## Final prepared analysis

- **Comparison:** 484 matched controls (reference) versus 508 children with
  moderate-to-severe diarrhea (MSD; comparison), 992 children in total.
- **Tested panel:** 78 published genera.
- **Adjustment model:** `~ age_z + country + group`, where age in months is
  standardized within the final cohort and country is categorical.
- **Independent unit:** one published study observation per enrolled child;
  the source analysis object does not contribute repeated observations for a
  child to this comparison.
- **Original library size:** the published `totalCounts` value, defined before
  the 5% genus-prevalence and two-group support filters.

Download `forserveroptim.rdata` from the MSD1000 repository and place it at
`raw/forserveroptim.rdata`. The cohort contains children with
moderate-to-severe diarrhea and matched controls
from Bangladesh, Kenya, Mali, and The Gambia. Controls are the reference group.
The model adjusts for age in months and country.

The script uses the published unrarefied genus counts and the published
`totalCounts` value as the original library size. It retains genera present in
at least 5% of the final cohort and observed in both comparison groups.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/gems_pediatric_diarrhea/prepare_data.R
```

The DASRA input and readable count, metadata, taxonomy, and preprocessing
tables are written to `processed/`.
