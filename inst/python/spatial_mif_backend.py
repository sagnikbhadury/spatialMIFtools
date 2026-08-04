"""Python backend extracted from the validated spatial mIF analysis pipeline."""

from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


def grid_image(frame: pd.DataFrame, phenotypes: Iterable[str], nx: int = 12,
               ny: int = 9, minimum_cells: int = 5,
               scale_factor: float = 10000.0) -> pd.DataFrame:
    """Convert physical cells with binary annotations into occupied grid bins."""
    phenotypes = list(phenotypes)
    required = {"patient", "image", "group", "x", "y", *phenotypes}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    records = []
    for (_, _), d in frame.groupby(["patient", "image"], sort=False):
        d = d.copy()
        x = pd.to_numeric(d["x"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d["y"], errors="coerce").to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        d, x, y = d.loc[ok].copy(), x[ok], y[ok]
        if d.empty:
            continue
        xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
        xr, yr = max(xmax - xmin, 1.0), max(ymax - ymin, 1.0)
        d["bin_x"] = np.minimum((nx * (x - xmin) / xr).astype(int), nx - 1)
        d["bin_y"] = np.minimum((ny * (y - ymin) / yr).astype(int), ny - 1)
        agg = {feature: "sum" for feature in phenotypes}
        agg["x"] = "count"
        grid = d.groupby(["bin_x", "bin_y"], as_index=False).agg(agg)
        grid = grid.rename(columns={"x": "n_cells"})
        grid = grid.loc[grid.n_cells >= minimum_cells].copy()
        if grid.empty:
            continue
        area = max((xr / nx) * (yr / ny), 1.0)
        grid["x"] = (grid.bin_x + 0.5) / nx
        grid["y"] = (grid.bin_y + 0.5) / ny
        for feature in phenotypes:
            grid[f"density__{feature}"] = np.log1p(grid[feature] / area * scale_factor)
            grid[f"prop__{feature}"] = grid[feature] / grid.n_cells
        grid["patient"] = str(d.patient.iloc[0])
        grid["image"] = str(d.image.iloc[0])
        grid["group"] = str(d.group.iloc[0])
        grid["image_cells"] = len(d)
        grid["bin_area"] = area
        records.append(grid)
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()


def add_epithelial_zones(grids: pd.DataFrame, epithelial: str = "Epithelial",
                         nx: int = 12, ny: int = 9,
                         sigma: float = 1.25) -> pd.DataFrame:
    """Create deterministic within-image epithelial KDE quartiles."""
    output = []
    labels = ["Q1_low", "Q2", "Q3", "Q4_high"]
    for _, d in grids.groupby("image", sort=False):
        d = d.copy()
        surface = np.zeros((ny, nx), float)
        surface[d.bin_y.to_numpy(int), d.bin_x.to_numpy(int)] = d[epithelial].to_numpy(float)
        kde = gaussian_filter(surface, sigma=sigma, mode="nearest")
        values = kde[d.bin_y.to_numpy(int), d.bin_x.to_numpy(int)]
        d["epithelial_kde"] = values
        if np.nanmax(values) <= 0 or np.nanstd(values) == 0:
            d["tumor_zone"] = "Q0_no_epithelial"
        else:
            ranks = pd.Series(values).rank(method="first")
            d["tumor_zone"] = pd.qcut(ranks, 4, labels=labels).astype(str).to_numpy()
        output.append(d)
    return pd.concat(output, ignore_index=True)


def partial_correlation(x: np.ndarray, ridge: float = 0.15) -> np.ndarray:
    """Ridge partial-correlation matrix used by the manuscript pipeline."""
    x = np.asarray(x, float)
    sd = x.std(axis=0, ddof=1)
    active = np.isfinite(sd) & (sd > 0)
    result = np.full((x.shape[1], x.shape[1]), np.nan)
    np.fill_diagonal(result, 1.0)
    if active.sum() < 2:
        return result
    z = (x[:, active] - x[:, active].mean(axis=0)) / sd[active]
    covariance = (z.T @ z) / max(len(z) - 1, 1) + ridge * np.eye(active.sum())
    precision = np.linalg.inv(covariance)
    pcor = -precision / np.sqrt(np.outer(np.diag(precision), np.diag(precision)))
    np.fill_diagonal(pcor, 1.0)
    result[np.ix_(active, active)] = pcor
    return result


def residualize_cellularity(x: np.ndarray, total: np.ndarray) -> np.ndarray:
    """Residualize every phenotype field on log1p total cellularity."""
    x = np.asarray(x, float)
    total = np.log1p(np.asarray(total, float))
    design = np.column_stack([np.ones(len(total)), total])
    coefficients = np.linalg.lstsq(design, x, rcond=None)[0]
    return x - design @ coefficients


def estimate_networks(grids: pd.DataFrame, phenotypes: Iterable[str],
                      ridge: float = 0.15, minimum_bins: int = 12,
                      adjust_cellularity: bool = False) -> pd.DataFrame:
    """Estimate long-form image-level ridge networks."""
    phenotypes = list(phenotypes)
    columns = [f"density__{p}" for p in phenotypes]
    rows = []
    for keys, d in grids.groupby(["group", "patient", "image"], sort=False):
        if len(d) < minimum_bins:
            continue
        x = d[columns].to_numpy(float)
        if adjust_cellularity:
            x = residualize_cellularity(x, d.n_cells.to_numpy(float))
        pcor = partial_correlation(x, ridge)
        for i, j in combinations(range(len(phenotypes)), 2):
            rows.append({"group": keys[0], "patient": keys[1], "image": keys[2],
                         "source": phenotypes[i], "target": phenotypes[j],
                         "edge": f"{phenotypes[i]}--{phenotypes[j]}",
                         "partial_correlation": pcor[i, j], "n_bins": len(d),
                         "cellularity_adjusted": adjust_cellularity})
    return pd.DataFrame(rows)


def cross_pair_enrichment(frame: pd.DataFrame, phenotype_a: str,
                          phenotype_b: str, radii=(0.025, 0.05, 0.10),
                          permutations: int = 99, seed: int = 20260803) -> pd.DataFrame:
    """Image-level random-label cross-pair enrichment."""
    rng = np.random.default_rng(seed)
    rows = []
    for keys, d in frame.groupby(["group", "patient", "image"], sort=False):
        xy = d[["x", "y"]].to_numpy(float)
        span = np.maximum(xy.max(axis=0) - xy.min(axis=0), 1)
        xy = (xy - xy.min(axis=0)) / span
        sparse = cKDTree(xy).sparse_distance_matrix(cKDTree(xy), max(radii),
                                                    output_type="coo_matrix")
        keep = sparse.row < sparse.col
        ii, jj, distance = sparse.row[keep], sparse.col[keep], sparse.data[keep]
        a, b = d[phenotype_a].to_numpy(bool), d[phenotype_b].to_numpy(bool)
        observed_cross = (a[ii] & b[jj]) | (b[ii] & a[jj])
        null = np.zeros((permutations, len(radii)))
        for p in range(permutations):
            pa, pb = rng.permutation(a), rng.permutation(b)
            cross = (pa[ii] & pb[jj]) | (pb[ii] & pa[jj])
            null[p] = [np.sum(cross & (distance <= r)) for r in radii]
        for j, radius in enumerate(radii):
            observed = np.sum(observed_cross & (distance <= radius))
            expected = null[:, j].mean()
            rows.append({"group": keys[0], "patient": keys[1], "image": keys[2],
                         "pair": f"{phenotype_a}--{phenotype_b}",
                         "radius_fraction": radius, "observed_pairs": observed,
                         "random_label_expected": expected,
                         "log2_enrichment": np.log2((observed + .5) / (expected + .5))})
    return pd.DataFrame(rows)
