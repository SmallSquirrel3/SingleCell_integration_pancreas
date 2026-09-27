"""
Reusable functions for single-cell QC.

Why MAD (median absolute deviation) instead of hard thresholds:
hard cutoffs (e.g. "remove everything with mt% > 10") don't account for the fact that
QC-metric distributions depend heavily on tissue type, sequencing depth, and the biology
of the specific experiment. A MAD-based approach finds outliers relative to the median of
the dataset itself, which is more robust and less likely to discard an entire rare cell
population just because of an overall shift in the distribution.

Approach reference: scverse best practices (sc-best-practices.org).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def calculate_qc_metrics(
    adata: AnnData,
    mt_pattern: str = "^MT-",
    inplace: bool = True,
) -> AnnData:
    """Compute standard RNA QC metrics: n_genes, n_counts, pct_counts_mt.

    Parameters
    ----------
    adata : AnnData
        Data with raw counts in .X.
    mt_pattern : str
        Regex for mitochondrial genes ("^MT-" for human, "^mt-" for mouse).
    inplace : bool
        Whether to modify adata in place.
    """
    if not inplace:
        adata = adata.copy()

    adata.var["mt"] = adata.var_names.str.match(mt_pattern)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )
    return adata


def mad_outlier(adata: AnnData, metric: str, n_mads: float = 5.0) -> pd.Series:
    """MAD-based outlier detection for a single QC metric.

    Returns a boolean mask: True = the cell is an outlier (a candidate for removal).
    """
    values = adata.obs[metric]
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    # 1.4826 -- the coefficient that scales MAD to the std of a normal distribution
    threshold = n_mads * mad * 1.4826
    return (values < median - threshold) | (values > median + threshold)


def apply_mad_filtering(
    adata: AnnData,
    metrics: dict[str, float] | None = None,
) -> tuple[AnnData, pd.Series]:
    """Apply MAD-based filtering across several metrics at once.

    Parameters
    ----------
    metrics : dict
        {metric_name: n_mads}. Defaults to permissive thresholds chosen so as not to lose
        rare cell populations.
    """
    if metrics is None:
        metrics = {
            "log1p_total_counts": 5,
            "log1p_n_genes_by_counts": 5,
            "pct_counts_mt": 3,
        }

    outlier_mask = pd.Series(False, index=adata.obs_names)
    for metric, n_mads in metrics.items():
        outlier_mask |= mad_outlier(adata, metric, n_mads)

    adata_filtered = adata[~outlier_mask].copy()
    return adata_filtered, outlier_mask


def protein_qc_summary(adata: AnnData, protein_key: str = "protein_expression") -> pd.DataFrame:
    """Summary of ADT (protein) counts per cell -- for CITE-seq data."""
    protein_counts = adata.obsm[protein_key]
    if hasattr(protein_counts, "toarray"):
        protein_counts = protein_counts.toarray()

    per_cell_total = protein_counts.sum(axis=1)
    return pd.DataFrame(
        {
            "min": [per_cell_total.min()],
            "median": [np.median(per_cell_total)],
            "max": [per_cell_total.max()],
            "n_proteins": [protein_counts.shape[1]],
        }
    )
