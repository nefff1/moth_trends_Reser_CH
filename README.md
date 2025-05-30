# Analyses of moth trends across 50 years in dependence of elevation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14506883.svg)](https://doi.org/10.5281/zenodo.14506883)

This repository contains codes and some data that were used in the analyses for the following manuscript:

Neff F, Chittaro Y, Korner-Nievergelt F, Litsios G, Martínez-Núñez C, Rey E, Knop E. **Contrasting 50-year trends of moth communities depending on elevation and species traits**

The analyses are based on a vast moth community dataset from Switzerland collected by Dr. Ladislaus Rezbanyai-Reser:

The following R files are included in the folder *R_Code*:

-   **R_prepare_data.R**: File used to prepare all data frames used for the analyses.

-   **R_analyses.R**: Code used to run models and produce outputs.

The following Stan code files are included in the folder *Stan_Code*:

-   **Stan_hg_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a hurdle gamma distribution.

-   **Stan_hg_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a hurdle gamma distribution.

-   **Stan_nb_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a zero-inflated negative binomial distribution.

-   **Stan_nb_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a zero-inflated negative binomial distribution.

Besides data available from the manuscript directly or from other sources indicated in the code (e.g. GBIF), the folder *Data* contains:

-   **d_samplings.txt**: Details on sampling site-year combinations used in the analyses: Spatio-temporal clusters, sampling pairs, different land-use proportions.
-   **d_nullnights.txt**: List of nights in which nothing was caught (missing from GBIF dataset).
-   **d_taxonomy.txt**: Moth species names according to the taxonomy used in the current analyses. Can be joint to GBIF data through the *taxonID* variable.
-   **d_mass.txt**: Estimated species-level biomass used to estimate community-wide biomass. Based on a set of allometric relationships.
