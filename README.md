# Analyses of moth trends across 50 years in dependence of elevation

[![](https://zenodo.org/badge/DOI/10.5281/zenodo.14506883.svg)](https://doi.org/10.5281/zenodo.14506883)

This repository contains codes and some data that were used in the analyses for the following manuscript:

Neff F, Chittaro Y, Korner-Nievergelt F, Litsios G, Martínez-Núñez C, Rey E, Knop E. **Contrasting 50-year trends of moth communities depending on elevation and species traits**

In this study, 50-year trends in moth communities (abundance, richness, biomass) in dependence of elevation and functional traits (body size, temperature niche, food specialisation, hibernation stage) were investigated. The analyses are based on a vast moth community dataset from Switzerland collected by Dr. Ladislaus Rezbanyai-Reser between 1972 and 2021. They are part of the INSECT project.

The following R files are included in the folder *R_Code*:

-   **R_prepare_data.R**: File used to prepare all data frames used for the analyses.

-   **R_analyses.R**: Code used to run models and produce outputs.

The following Stan code files are included in the folder *Stan_Code*. All code is adapted from models built through 'brms' (Bürkner et al. 2022):

-   **Stan_hg_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a hurdle gamma distribution.

-   **Stan_hg_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a hurdle gamma distribution.

-   **Stan_nb_spline_s1_r4.stan**: Stan model code of a regression model with 1 smoothing term, 4 random terms and a zero-inflated negative binomial distribution.

-   **Stan_nb_spline_s1p1_r4.stan**: Stan model code of a regression model with 1 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a zero-inflated negative binomial distribution.

Additional to species records data available from GBIF (<https://doi.org/10.15468/dl.gcagva>), the folder *Data* contains:

-   **d_samplings.txt**: Details on sampling site-year combinations used in the analyses. Rows are unique site and year combinations. Columns are:

    -   LOC: Study site (ID)

    -   A: Year of sampling (integer)

    -   spattemp_cluster: Spatio-temporal clusters of sites in years (categorical)

    -   samplingpair: sampling pairs of simultaneously opterated sites (categorical)

    -   height: elevation of study site (meters above sea level)

-   **d_nullnights.txt**: List of nights in which nothing was caught (missing from the GBIF dataset). Rows are single sampling nights (unique combination of site and date). Columns are:

    -   LOC: Study site (ID)

    -   A: Year of sampling (integer)

    -   Samplingdate: Date of sampling (Y-m-d)

    -   active_hours: Number of hours for which a trap was active (real number)

-   **d_taxonomy.txt**: Moth species names according to the taxonomy used in the current analyses. Can be joint to GBIF data through the *taxonID* variable. Each row is a species. Columns are:

    -   taxonID: taxonID as reported in the GBIF data

    -   Name_std: Species name used in the analyses

-   **d_mass.txt**: Estimated species-level biomass used to estimate community-wide biomass. Based on wingspan data (Jonko 2002–2024, Fibiger 1990, Potocký et al. 2018, Ronkay et al. 2001) and a set of allometric relationships (Kinsella et al. 2020). Each row is a species. Columns are:

    -   Name_std: Species name used in the analyses

    -   mass: estimated mass (real number; grams)

-   **d_weather.txt**: Temperature and precipitation data for all 35,847 sampling nights. Calculated from gridded daily temperature and precipitation data from MeteoSwiss ([https://www.meteoswiss.admin.ch](https://www.meteoswiss.admin.ch/); TabsD and RhiresD; licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). Rows are sampling nights (unique combination of site and date). Columns are:

    -   LOC: Study site (ID)

    -   Samplingdate: Date of sampling (Y-m-d)

    -   T_2day: Averaged mean daily temperature of the two days that include the sampling night (real number; degree Celsius)

    -   P_2day: Summed precipitation of the two days that include the sampling night (real number; mm)

-   **d_traits.txt**: Species-level trait data (body size, temperature niche, specialisation, overwintering stage). Body size, specialisation and overwintering stage data were collated and summarised from Cook et al. (2022), Fibiger (1990), Hacker & Müller (2006), Jonko et al. (2002–2024), Lepiforum e.V. (2002–2021), Mangels et al. (2017), Pearse & Altermatt (2013), Potocký et al. (2018), Ronkay et al. (2001), Steiner et al. (2014), Ziegler (2005–2024). Temperature niches were derived from GBIF distribution data (<https://doi.org/10.15468/dl.2mev52>, <https://doi.org/10.15468/dl.km9rkn>) and WorldlClim2 (Fick & Hijmans 2017). Each row is a species. Columns are:

    -   Species: Species name used in the analyses

    -   Body size: body size class (categorial; small, medium, large)

    -   Temp. niche: temperature niche class (categorical; cold, intermediate, warm)

    -   Specialisation: food specialisation class (categorical; monophagous, oligophagous, polyphagous)

    -   Overw. stage: overwintering stage (categorical; egg, egg/larva, larva, larva/pupa, pupa, adult)

**References**

Bürkner, P.-C., Gabry, J., Weber, S., Johnson, A., Modrák, M., Badr, H. S., Weber, F., Ben-Shachar, M. S., & Rabel, H. (2022). *Bayesian Regression Models using “Stan”. R package version 2.17.0.* [Computer software]. <https://cran.r-project.org/web/packages/brms/brms.pdf>

Cook, P. M., Tordoff, G. M., Davis, T. M., Parsons, M. S., Dennis, E. B., Fox, R., Botham, M. S., & Bourn, N. A. D. (2022). Traits data for the butterflies and macro-moths of Great Britain and Ireland. *Ecology*, *103*(5), e3670. <https://doi.org/10.1002/ecy.3670>

Fibiger, M. (Ed.). (1990). *Noctuidae Europaeae. Volume 1. Noctuinae I*. Entomological Press.

Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, *37*(12), 4302–4315. <https://doi.org/10.1002/joc.5086>

Hacker, H. H., & Müller, J. (2006). *Die Schmetterlinge der bayerischen Naturwaldreservate: Eine Charakterisierung der süddeutschen Waldlebensraumtypen anhand der Lepidoptera (Insecta)*. Arbeitsgemeinschaft Bayer. Entomologen.

Jonko, C. (2002–2024). *Lepidoptera Mundi*. <https://lepidoptera.eu/>

Kinsella, R. S., Thomas, C. D., Crawford, T. J., Hill, J. K., Mayhew, P. J., & Macgregor, C. J. (2020). Unlocking the potential of historical abundance datasets to study biomass change in flying insects. *Ecology and Evolution*, *10*(15), 8394–8404. <https://doi.org/10.1002/ece3.6546>

Lepiforum e.V. (2002–2021). *Lepiforum e.V. – Bestimmung von Schmetterlingen und ihren Präimaginalstadien*. <https://lepiforum.org/>

Mangels, J., Fiedler, K., Schneider, F. D., & Blüthgen, N. (2017). Diversity and trait composition of moths respond to land-use intensification in grasslands: Generalists replace specialists. *Biodiversity and Conservation*, *26*(14), 3385–3405. <https://doi.org/10.1007/s10531-017-1411-z>

Pearse, I. S., & Altermatt, F. (2013). Predicting novel trophic interactions in a non-native world. *Ecology Letters*, *16*(8), 1088–1094. <https://doi.org/10.1111/ele.12143>

Potocký, P., Bartoňová, A., Beneš, J., Zapletal, M., & Konvička, M. (2018). Life-history traits of Central European moths: Gradients of variation and their association with rarity and threats. *Insect Conservation and Diversity*, *11*(5), 493–505. <https://doi.org/10.1111/icad.12291>

Ronkay, L., Yela Garcia, J. L., & Hreblay, M. (2001). *Noctuidae Europaeae. Volume 5. Hadeninae II*. Entomological Press.

Steiner, A., Ratzel, U., Top-Jensen, M., & Fibiger, M. (Eds.). (2014). *Die Nachtfalter Deutschlands: Ein Feldführer*. BugBook Publishing.

Ziegler, H. (2005–2024). *Butterflies & Moths of Palaearctic Regions*. <https://euroleps.ch>
