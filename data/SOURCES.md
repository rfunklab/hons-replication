# Data sources and formats

## Applied Mathematics diagrams and grid

`applied_math_diagrams.npz` contains the authors' saved H1/H2 persistence intervals for the 144 Applied Mathematics threshold pairs in Figures 8, 9 and 11. The source corpus comprises 167,001 articles from the September 2021 Dimensions snapshot.

Each array is named `upper_index_lower_index` and has columns `(dimension, birth, death)`. Dimensions are 1 and 2; births and finite deaths are normalized to 0–1. Infinite deaths remain `inf` and are replaced by 1.00001 when constructing persistence images. The archive holds 9,104,947 intervals. Empty arrays have shape `(0, 3)`.

The original network calculation retained concepts within each frequency band, dated edges by their first co-occurrence and expanded cliques through tetrahedra. Vertices entered at zero. Each retained network used its own earliest retained article/concept year and the corpus's final year, 2021, to normalize time; homology used coefficients in Z/11Z. The saved intervals preserve that preparation. Reproduction starts with these intervals.

`applied_math_grid.csv` gives each array's two indices, lower/upper article-frequency bounds, retained node/edge counts and first year. The 16 × 9 grid includes a first lower-bound column needed for neighboring comparisons; the other 128 cells are selectable. The diagram archive's SHA-256 is `90a4582e197a4203ae57c8024dc323eb9efd2735fb6a0f2f446b6de0989c87e5`.

## Comparison measurements

`applied_math_survival.csv` contains finite-feature survival counts for the primary HONS selection and the ten alternatives in Figure 12. Columns identify the method, homological dimension, minimum lifetime in integer years, feature count and denominator. H1 uses retained edges; H2 uses retained triangles. The count at year zero gives the finite-feature count. Summing survival counts at years 1 and above gives total finite persistence in years because the source publication dates are integer years. Counts extend through the full observed lifetime range, although the figure displays 0–40 years.

`comparison_networks.csv` contains 1,672 network-level rows: the primary HONS selection and ten alternatives for 152 completed scientific fields. Columns preserve field and method identifiers, selected configuration, filter setting, frequency bounds, network sizes, finite H1/H2 counts, summed finite lifetimes in years and the filtration span. The `field_code` values are ANZSRC 2008 Field of Research codes with any leading zero dropped, so `102` denotes group 0102, Applied Mathematics. The table retains the fields with zero-edge HONS selections; the analysis derives the eligible set from the measurements. Comparison denominators and positive-pair rules determine the reported sample sizes.

These two tables contain the counts, lifetimes and network sizes used in the paper's comparisons. The scripts recalculate the reported densities and ratio summaries from the included counts, lifetimes and denominators. The original filter execution and persistent-homology calculations precede these inputs. The inputs for Applied Mathematics and the cross-field comparisons were prepared as described in the paper's comparison-methods appendix.

The generated `output/table_F2.csv` contains the cross-field summaries. For each measure and alternative, `n_R` counts fields with positive values for both HONS and the alternative; `n_P` counts all fields with defined values, including zeros. Median and quartile ratios are calculated on the log2 scale and transformed back. Win/tie/loss counts use all `n_P` pairs; values are tied when their absolute difference is at most `1e-12 + 1e-9 * abs(alternative_value)`.

## Workplace observations

`workplace_contacts.tsv` is the unmodified contact list `tij_InVS.dat` distributed in the [SocioPatterns workplace dataset](https://sociopatterns.org/datasets/contacts-in-a-workplace/). The bundled file has 9,827 rows and SHA-256 `ba0c81bd10830e688abb3711647e5dc6cecb3aa18e962039ed7ca9242603fd37`.

The three whitespace-separated columns are timestamp `t` in seconds and anonymous participant identifiers `i` and `j`. Each row records a contact during `[t−20, t]`. The observations were collected in 2013 in an office building in France. Department attributes are unused and are not bundled.

The dataset is distributed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/). Cite M. Génois et al., “Data on face-to-face contacts in an office building suggest a low-cost vaccination strategy based on community linkers,” *Network Science* 3, 326–347 (2015).

## Code and data

The repository's MIT license covers its code, documentation, original assets and derived numerical files, including the manuscript inputs in this directory and the analysis outputs built from them. SocioPatterns provides the workplace data under CC0. The manuscript inputs are derived numerical results from the authors' analysis; the source Dimensions article records remain subject to their separate access terms and are not included.
