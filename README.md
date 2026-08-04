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
