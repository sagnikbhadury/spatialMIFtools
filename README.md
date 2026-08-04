# spatialMIFtools

`spatialMIFtools` is an R package for patient-level analysis of two-dimensional
multiplex tissue imaging data. It implements the reusable analysis components
from the CRC CLR/DII and pancreatic IPMN/PDAC manuscript while deliberately
excluding governed patient data and manuscript-specific fitted objects.

## Statistical principles

- The patient is the independent sampling unit.
- Images are repeated fields nested within patients; bins are repeated spatial
  observations within images.
- Composition, marginal proximity, recurrent neighborhoods, and conditional
  networks are kept as distinct estimands.
- Pathology contrasts are formed after image estimates are aggregated within
  patient or through patient-clustered GEE.
- Cross-validation learns preprocessing, neighborhood number, centers,
  imputation/scaling, and prediction models inside each training fold.

## Installation

This repository is private. Authenticated GitHub users with access can install
with:

```r
# install.packages("pak")
pak::pak("sagnikbhadury/spatialMIFtools")
```

## Minimal workflow

```r
library(spatialMIFtools)

phenotypes <- c("Epithelial", "APC", "Treg", "CTL", "HelperT")
cells <- simulate_mif_cells(n_patients = 12, images_per_patient = 2)

fit <- run_spatial_mif_pipeline(
  cells,
  phenotypes = phenotypes,
  group_levels = c("Group1", "Group2")
)

fit$adjusted_contrasts
```

The cell table must have one row per physical cell with `patient`, `image`,
`group`, `x`, `y`, and binary phenotype columns. A positive-marker event table
must first be reconstructed to physical cells; the package does not silently
treat marker events as separate cells.

## Main functions

| Task | Function |
|---|---|
| Validate physical cells | `validate_mif_cells()` |
| Area-normalized grids | `bin_mif_cells()` |
| Epithelial KDE quartiles | `add_epithelial_zones()` |
| Recurrent neighborhoods | `fit_spatial_neighborhoods()` |
| Ridge conditional networks | `estimate_ridge_networks()` |
| Cellularity adjustment | `residualize_cellularity()` |
| Patient edge contrasts | `patient_edge_contrasts()` |
| Ordered Q1--Q4 GEE | `ordered_zone_gee()` |
| Random-label proximity | `cross_pair_enrichment()` |
| Patient-label permutation | `patient_label_permutation()` |
| Fully nested prediction | `nested_lopo_prediction()` |
| Complete core workflow | `run_spatial_mif_pipeline()` |

## Python backend

The package also includes a Python backend extracted from the validated
analysis scripts. It is optional and accessed through `reticulate`:

```r
python_backend_available()
py <- spatial_mif_python()
py$partial_correlation(matrix(rnorm(500), ncol = 5), ridge = 0.15)
```

Python requirements are listed in `inst/python/requirements.txt`. The native-R
and Python implementations are tested against common synthetic inputs.

## ISPAT and GP-GHS

`spatial_methods_workbench_adapter()` calls an installed
`SpatialMethodsWorkbench` package. ISPAT and GP-GHS source code is not vendored;
their method-specific validation remains with their original packages and
publications.

## Scope

Edges are partial correlations between local phenotype-density fields. They do
not establish physical contact, directionality, signaling, cellular function,
or causality. Q1--Q4 are relative within-image epithelial-density ranks and do
not represent temporal disease stages.

## Complete implementation guide

### 1. Software requirements

The native implementation requires R 4.2 or newer. Install the optional
components needed by the analysis you intend to run:

```r
install.packages(c(
  "geepack",     # ordered patient-clustered GEE
  "reticulate",  # bundled Python backend
  "knitr",       # vignette
  "rmarkdown"
))
```

The Python backend was tested with Python 3.12 and requires NumPy, pandas,
SciPy, scikit-learn, and statsmodels:

```bash
python -m pip install -r inst/python/requirements.txt
```

Point `reticulate` to the intended environment before loading the backend:

```r
library(reticulate)
use_python("C:/path/to/python", required = TRUE)
library(spatialMIFtools)
stopifnot(python_backend_available())
```

### 2. Required physical-cell data model

The minimum cell table is:

| Column | Meaning | Required type |
|---|---|---|
| `patient` | Independent sampling unit | character/factor |
| `image` | Tissue image or ROI nested within patient | character/factor |
| `group` | Pathology or comparison group | character/factor |
| `x`, `y` | Cell centroid in the original image coordinates | finite numeric |
| phenotype columns | Supplied cell/marker annotations | binary 0/1 |

Example:

```r
phenotypes <- c("Epithelial", "APC", "Treg", "CTL", "HelperT")
validate_mif_cells(cells, phenotypes, allow_overlap = TRUE)
```

The package permits overlapping binary annotations because some source panels
assign multiple positive markers to one physical cell. If the source is a
positive-marker event export with several rows per physical cell, reconstruct
one row per physical cell before calling the package. The package intentionally
does not guess a cell identifier or collapse event rows automatically.

Do not place patient identifiers or governed raw files in an R package
repository. Read them from an approved external location at analysis time.

### 3. Phenotype mapping and provenance

Create the requested analysis compartments before gridding and retain the
mapping in the project-level analysis protocol. For the manuscript, the shared
labels were Epithelial, APC, Treg, CTL, and HelperT, but their source definitions
were cohort-specific. A statistical label such as APC should not be interpreted
as a direct functional antigen-presentation measurement unless the assay
supports that claim.

```r
# Example mapping only; adapt to the governed source schema.
cells$Epithelial <- as.integer(cells$source_label == "tumor cells")
cells$CTL <- as.integer(cells$source_label == "CD8+ T cells")
```

Avoid reannotation or proxy substitution unless it is explicitly part of the
analysis protocol. A phenotype with no positive calls is non-estimable and
should not be replaced with an unrelated marker.

### 4. Spatial gridding

```r
grids <- bin_mif_cells(
  cells,
  phenotypes = phenotypes,
  nx = 12,
  ny = 9,
  min_cells = 5,
  scale_factor = 1e4
)
```

For phenotype `f` in bin `b`, the density variable is

```text
log(1 + scale_factor * phenotype_count / bin_area).
```

`bin_area` is calculated in the original image coordinate system. Bin-center
coordinates are separately normalized to the unit square for spatial
calculations. The density is not a within-bin cell-type proportion. Changing
grid resolution changes the local scale and should be treated as a sensitivity
analysis.

### 5. Epithelial-density coordinate

```r
zones <- add_epithelial_zones(
  grids,
  epithelial = "Epithelial",
  sigma = 1.25,
  nx = 12,
  ny = 9
)
```

The function places epithelial counts on the grid, smooths the image with a
Gaussian kernel, ranks the smoothed values at occupied bins, and partitions the
ranks into Q1--Q4. Ties are resolved deterministically before quartiling.

- Q1 and Q4 mean relatively low and high epithelial density within the same
  image.
- Quartiles are not common absolute density thresholds.
- Quartiles are not stages of a longitudinal disease trajectory.
- Images with no varying epithelial surface are marked `Q0_no_epithelial`.

Repeat the analysis across prespecified grid sizes and kernel widths when the
conclusion could depend on spatial scale.

### 6. Spatial autocorrelation

```r
one_image <- subset(zones, image == unique(image)[1])
qc <- four_nn_autocorrelation(
  one_image$density__Epithelial,
  one_image[c("x", "y")],
  permutations = 199,
  seed = 20260803
)
```

The score uses four nearest bin centers and is proportional to Moran's I on a
regular four-neighbor graph. It does not include normalization by total spatial
weight and should not be labeled formal Moran's I.

### 7. Recurrent spatial neighborhoods

```r
neighborhood_fit <- fit_spatial_neighborhoods(
  zones,
  phenotypes,
  k_range = 3:7,
  nstart = 50,
  max_silhouette_bins = 6000,
  seed = 20260803
)

assigned <- assign_neighborhoods(neighborhood_fit, zones)
prevalence <- patient_neighborhood_prevalence(assigned)
neighborhood_fit$scores
```

Phenotype densities are standardized before K-means. `K` is selected by the
largest sampled silhouette coefficient. Cluster names are deliberately not
assigned by the fitting function: biological or descriptive names should be
applied after inspecting the fitted centers, so names cannot influence
assignments.

For confirmatory work, bootstrap patients within pathology and images within
patient, refit the standardizer and centers, optimally match bootstrap centers
to the reference, and report center-selection frequencies. The manuscript used
matched-center correlation at least 0.8 in at least 80% of replicates as its
prespecified reproducibility rule.

### 8. Ridge conditional networks

```r
networks <- estimate_ridge_networks(
  assigned,
  phenotypes,
  ridge = 0.15,
  minimum_bins = 12,
  adjust_cellularity = FALSE
)
```

Within each image, standardized phenotype-density covariance is regularized by
adding `ridge * I` before inversion. Off-diagonal precision elements are
converted to partial correlations. Constant phenotypes are non-estimable.

```r
pcor <- ridge_partial_correlation(
  as.matrix(assigned[paste0("density__", phenotypes)]),
  ridge = 0.15
)
```

A positive edge indicates positive conditional association between two density
fields after conditioning on the other included fields. It does not imply cell
contact, direction, signaling, or causal interaction.

### 9. Total-cellularity adjustment

```r
adjusted_networks <- estimate_ridge_networks(
  assigned,
  phenotypes,
  ridge = 0.15,
  minimum_bins = 12,
  adjust_cellularity = TRUE
)
```

For each phenotype and image, the package fits a linear nuisance regression of
log phenotype density on `log1p(n_cells)` and estimates the network from the
residual fields. This asks whether an edge remains beyond a linear association
with local total cellularity. It does not control omitted morphology, unmeasured
cell types, nonlinear cellularity effects, or acquisition batch.

### 10. Patient-level edge contrasts

Always specify contrast orientation explicitly:

```r
edge_results <- patient_edge_contrasts(
  adjusted_networks,
  group_levels = c("IPMN", "PDAC"),  # returns PDAC minus IPMN
  p_adjust = "BH",
  bootstrap_replicates = 2000,
  seed = 20260803
)
```

Image-specific edges are averaged within patient before Welch comparison.
Multiplicity is controlled over the returned edge family. The number of images
or bins is not used as the independent sample size.

### 11. Ordered Q1--Q4 inference

Estimate zone-specific networks with:

```r
zone_networks <- estimate_ridge_networks(
  assigned,
  phenotypes,
  ridge = 0.15,
  minimum_bins = 8,
  by_zone = TRUE
)
```

Then call the patient-clustered GEE separately for each prespecified edge or
outcome:

```r
target <- subset(zone_networks, edge == "Epithelial--APC")
ordered_zone_gee(
  target,
  value = "partial_correlation",
  group_levels = c("CLR", "DII")
)
```

The fitted model is `outcome ~ group + zone + group:zone`, with zone coded
1--4, an exchangeable working correlation, and patient-clustered sandwich
standard errors. The interaction coefficient is the pathology difference in
linear change per epithelial-density quartile. Apply BH correction across the
prespecified edge family outside the single-edge call.

### 12. Random-label cross-pair enrichment

```r
proximity <- cross_pair_enrichment(
  cells,
  pairs = list(c("Epithelial", "CTL"),
               c("Epithelial", "HelperT"),
               c("CTL", "HelperT")),
  radii = c(0.025, 0.05, 0.10),
  permutations = 99,
  seed = 20260803
)
```

The x and y axes are normalized separately by each image span. Observed
unordered cross-type pair counts are compared with random relabelings that
preserve phenotype margins. Enrichment is
`log2((observed + 0.5)/(expected + 0.5))`.

This statistic is a marginal random-label diagnostic. It is not formal Ripley's
cross-K because it lacks cross-K intensity normalization and boundary
correction. Average image enrichments within patient before group comparison.

### 13. Patient-label permutation

`patient_label_permutation()` accepts a statistic function that returns a named
numeric vector, making it applicable to global edges or ordered interactions:

```r
edge_statistic <- function(d) {
  result <- patient_edge_contrasts(d, c("Group1", "Group2"))
  setNames(result$statistic, result$edge)
}

perm <- patient_label_permutation(
  adjusted_networks,
  statistic = edge_statistic,
  permutations = 999,
  target = "CTL--HelperT",
  seed = 20260803
)

perm$target_p
perm$familywise_p
```

Labels are permuted at the patient level. The familywise empirical reference is
the maximum absolute statistic across the named statistic family in each
permutation.

### 14. Fully nested LOPO prediction

```r
scores <- nested_lopo_prediction(
  assigned,
  phenotypes,
  group_levels = c("IPMN", "PDAC"),
  feature_set = "combined",
  static_features = patient_network_features,
  k_range = 3:7,
  lambda = 4,
  seed = 20260803
)

attr(scores, "auc")
table(scores$selected_k)
```

Inside every training fold the function relearns:

1. phenotype-density means and standard deviations;
2. neighborhood number by silhouette score;
3. K-means centers;
4. held-out neighborhood assignments;
5. feature scaling; and
6. ridge-logistic coefficients.

`static_features` can contain one row per patient with prespecified network or
other image-aggregated summaries. They are joined only after the patient split.
The output is internal discrimination. External data are required before
diagnostic or biomarker claims.

### 15. Native-R versus Python execution

Use native R functions for a dependency-light package workflow. Use the Python
module when exact reuse of the manuscript's computational kernels or
interoperability with the original Python pipeline is desired:

```r
py <- spatial_mif_python("C:/path/to/python")

py_grids <- py$grid_image(
  cells,
  phenotypes,
  nx = 12L,
  ny = 9L,
  minimum_cells = 5L
)

py_zones <- py$add_epithelial_zones(
  py_grids,
  epithelial = "Epithelial",
  nx = 12L,
  ny = 9L,
  sigma = 1.25
)
```

The bundled module exposes `grid_image`, `add_epithelial_zones`,
`partial_correlation`, `residualize_cellularity`, `estimate_networks`, and
`cross_pair_enrichment`. The R and Python partial-correlation implementations
are tested numerically on a common matrix.

### 16. ISPAT and GP-GHS through SpatialMethodsWorkbench

Install the workbench separately and pass its method-specific arguments through
the adapter:

```r
result <- spatial_methods_workbench_adapter(
  data = workbench_input,
  method = "GPGHS",
  # method-specific parameters follow
  n_basis = 4
)
```

The adapter intentionally does not harmonize estimands across ridge, ISPAT, and
GP-GHS. Direction, magnitude-rank, and spatial-location concordance should be
reported separately because the methods condition and smooth differently.

### 17. Recommended project-level workflow

For each cohort:

1. Freeze and checksum the governed physical-cell input.
2. Record the phenotype mapping and all non-estimable annotations.
3. Validate cells and create the primary grid.
4. Construct within-image epithelial-density zones.
5. Fit and resample recurrent neighborhoods.
6. Estimate unadjusted and cellularity-adjusted image networks.
7. Aggregate images within patient or use patient-clustered GEE.
8. Control multiplicity within prespecified analysis families.
9. Run patient-label permutations for principal statistics.
10. Repeat prespecified grid and KDE scales.
11. Run fully nested patient-level prediction.
12. Use ISPAT/GP-GHS only for explicitly defined secondary comparisons.
13. Save session information, seeds, package versions, checksums, and derived
    nonidentifying tables.

Do not pool biologically different disease contrasts merely because the same
software was applied to both cohorts.

### 18. Reproducibility and package validation

Run the complete package checks locally:

```powershell
$env:SPATIALMIFTOOLS_PYTHON = 'C:\path\to\python.exe'
R CMD build .
R CMD check --no-manual spatialMIFtools_0.1.0.tar.gz
```

The checked release reports 0 errors, 0 warnings, and 0 notes. Continuous
integration repeats the R package check and installs the Python requirements.
See `CHECK_RESULTS.md` for the validated environment and equivalence tests.

### 19. Troubleshooting

**No bins returned**

Lower `min_cells`, verify the coordinate units, and confirm each image contains
finite coordinates. Do not combine unrelated images to increase bin counts.

**All zones are `Q0_no_epithelial`**

Verify that the selected epithelial column contains positive calls and that it
is a count column in the grid table.

**An edge is `NA`**

One phenotype is constant or too few bins remain in that image/zone. Treat the
edge as non-estimable; do not replace it with zero.

**K changes across LOPO folds**

Report the foldwise selection frequencies. A changing K is not software
failure, but it qualifies claims that neighborhoods are reproducible states.

**Python backend is unavailable**

Call `reticulate::py_config()`, explicitly select the environment with
`use_python()`, and install `inst/python/requirements.txt` into that same
environment.

**GEE fails to estimate the interaction**

Check that both groups contain multiple patients and that at least three ordered
zones are represented for the target summary. Do not treat bins as independent
clusters.

### 20. Reporting checklist

At minimum, report:

- patient, image, physical-cell, and retained-bin counts;
- phenotype definitions and non-estimable annotations;
- grid resolution, area transformation, minimum cells, and KDE width;
- neighborhood candidate K values, selection rule, and stability frequencies;
- ridge penalty, minimum network bins, and cellularity adjustment;
- contrast orientation, effect estimate, confidence interval, and multiplicity
  family;
- patient-level resampling or permutation design;
- all preprocessing learned inside prediction folds;
- AUC together with Brier score and calibration;
- spatial-scale sensitivity; and
- the distinction between statistical association and cellular function.
