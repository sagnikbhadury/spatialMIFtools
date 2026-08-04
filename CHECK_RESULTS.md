# Package validation

Validation date: 2026-08-04

- R version: 4.5.2 on Windows 11 x64
- Source package build: passed, including vignette construction
- `R CMD check --no-manual`: **Status: OK** (0 errors, 0 warnings, 0 notes)
- Native-R runtime workflow: passed on synthetic physical-cell data
- Ordered patient-clustered GEE: passed with `geepack` 1.3.13
- Fully nested LOPO workflow: passed
- R/Python ridge partial-correlation equivalence: maximum absolute difference
  below `1e-10` on the shared test matrix
- Python gridding bridge: passed with Python 3, NumPy, pandas, and SciPy

The formal check used `_R_CHECK_FORCE_SUGGESTS_=false` because
`SpatialMethodsWorkbench` is an optional external adapter rather than a hard
package dependency. Its method code is not vendored in this repository.
