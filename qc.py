"""
Modulo QC -- Control de calidad y análisis estadístico del pipeline AutoGCMS.

Funciones principales:
    TIC / QC técnico:
        bin_tic, build_tic_matrix, pca_tic_qc, plot_pca_tic, compute_qc_summary

    Preprocesado estadístico de features:
        infer_sample_metadata, preprocess_feature_matrix

    PCA de features (análisis principal):
        pca_features, plot_pca_features, compute_pca_loadings_table

    Heatmaps:
        plot_heatmap_main, plot_heatmap_features (z-score)

    Utilidades:
        validate_minimum_data, compute_advanced_qc_summary
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.typing import NDArray

import config
from preprocessing import GCMSData

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# TIC binning
# ─────────────────────────────────────────────────────────────────────────────

def bin_tic(
    data: GCMSData,
    bin_size_min: float = config.TIC_BIN_SIZE_MIN,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Discretiza el TIC corregido en bins de RT constantes.
    Cada bin contiene la MEDIA de las intensidades del TIC en ese intervalo.
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw
    rt = data.scan_times_min

    rt_min_val, rt_max_val = rt[0], rt[-1]
    bin_edges = np.arange(rt_min_val, rt_max_val + bin_size_min, bin_size_min)
    n_bins = len(bin_edges) - 1

    binned = np.zeros(n_bins, dtype=np.float64)
    for i in range(n_bins):
        mask = (rt >= bin_edges[i]) & (rt < bin_edges[i + 1])
        if mask.any():
            binned[i] = np.mean(tic[mask])

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    return bin_centers, binned


def build_tic_matrix(
    samples: dict[str, GCMSData],
    bin_size_min: float = config.TIC_BIN_SIZE_MIN,
) -> tuple[pd.DataFrame, NDArray[np.float64]]:
    """
    Construye matriz (muestras x bins) de TIC binned, normalizada por area total.
    Cada fila suma 1.0 (proporcion del area del TIC).
    """
    if not samples:
        return pd.DataFrame(), np.array([])

    global_rt_min = min(d.scan_times_min[0] for d in samples.values())
    global_rt_max = max(d.scan_times_min[-1] for d in samples.values())
    bin_edges = np.arange(global_rt_min, global_rt_max + bin_size_min, bin_size_min)
    n_bins = len(bin_edges) - 1
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    rows = {}
    for name, data in samples.items():
        tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw
        rt = data.scan_times_min
        binned = np.zeros(n_bins, dtype=np.float64)
        for i in range(n_bins):
            mask = (rt >= bin_edges[i]) & (rt < bin_edges[i + 1])
            if mask.any():
                binned[i] = np.mean(tic[mask])
        total = binned.sum()
        if total > 0:
            binned = binned / total
        rows[name] = binned

    df = pd.DataFrame(rows).T
    df.columns = [f"bin_{c:.3f}" for c in bin_centers]
    df.index.name = "sample"
    return df, bin_centers


# ─────────────────────────────────────────────────────────────────────────────
# PCA interno (numpy SVD como fallback)
# ─────────────────────────────────────────────────────────────────────────────

def _pca_numpy(X: NDArray, n_components: int) -> dict:
    """PCA via SVD con numpy (sin dependencia sklearn)."""
    X_centered = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    total_var = (S ** 2).sum()
    explained = (S[:n_components] ** 2) / total_var
    scores = X_centered @ Vt[:n_components].T
    return {
        "scores": scores,
        "loadings": Vt[:n_components],
        "explained_variance_ratio": explained,
    }


# ─────────────────────────────────────────────────────────────────────────────
# PCA de TIC (QC tecnico)
# ─────────────────────────────────────────────────────────────────────────────

def pca_tic_qc(
    tic_matrix: pd.DataFrame,
    n_components: int = config.PCA_N_COMPONENTS,
) -> dict | None:
    """
    PCA sobre la matriz de TIC binned.

    Returns dict con scores, loadings, explained_variance_ratio, sample_names.
    None si hay menos de 2 muestras.
    """
    n_samples = len(tic_matrix)
    if n_samples < 2:
        logger.warning("  PCA TIC requiere >= 2 muestras (hay %d). Saltando.", n_samples)
        return None

    X = tic_matrix.values.astype(np.float64)
    n_comp = min(n_components, n_samples - 1, X.shape[1])

    try:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_comp)
        scores = pca.fit_transform(X)
        result = {
            "scores": scores,
            "loadings": pca.components_,
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    except ImportError:
        logger.info("  sklearn no disponible, usando PCA numpy SVD.")
        result = _pca_numpy(X, n_comp)

    result["sample_names"] = list(tic_matrix.index)
    logger.info("  PCA TIC: %d componentes, varianza explicada = %s",
                n_comp, [f"{v:.1%}" for v in result["explained_variance_ratio"]])
    return result


def plot_pca_tic(
    pca_result: dict | None,
    sample_metadata: pd.DataFrame | None = None,
    save_path: str | Path | None = None,
) -> plt.Figure | None:
    """
    Scatter PC1 vs PC2 del TIC (QC tecnico), coloreado por tipo de muestra.

    Parameters
    ----------
    pca_result      : dict de pca_tic_qc() o None
    sample_metadata : DataFrame con columnas sample_id y sample_type (opcional)
    save_path       : ruta para guardar la figura
    """
    if pca_result is None:
        return None

    scores = pca_result["scores"]
    names = pca_result["sample_names"]
    evr = pca_result["explained_variance_ratio"]

    color_map = {
        "Sample": "steelblue", "QC": "orange",
        "Blank": "red", "Control": "green", "Unknown": "gray",
    }

    colors = ["steelblue"] * len(names)
    unique_types: set = set()
    if sample_metadata is not None:
        meta_dict = dict(zip(sample_metadata["sample_id"], sample_metadata["sample_type"]))
        colors = []
        for n in names:
            st = meta_dict.get(n, "Unknown")
            unique_types.add(st)
            colors.append(color_map.get(st, "gray"))
    else:
        unique_types = {"Sample"}

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(scores[:, 0], scores[:, 1],
               c=colors, s=70, edgecolor="white", alpha=0.85, zorder=3)

    # Etiquetar puntos si son pocos o si lo pide config
    if getattr(config, "PCA_LABEL_POINTS", False) or len(names) <= 12:
        for i, name in enumerate(names):
            ax.annotate(name, (scores[i, 0], scores[i, 1]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.75)

    ax.set_xlabel(f"PC1 ({evr[0]:.1%} varianza explicada)", fontsize=11)
    ax.set_ylabel(f"PC2 ({evr[1]:.1%} varianza explicada)" if len(evr) > 1 else "PC2", fontsize=11)
    ax.set_title("PCA QC — TIC binned (control de calidad tecnico)", fontsize=12)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.grid(True, alpha=0.3)

    if unique_types:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color_map.get(t, "gray"), edgecolor="white", label=t)
            for t in sorted(unique_types)
        ]
        ax.legend(handles=legend_elements, loc="best", fontsize=9, framealpha=0.8)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# QC summary por muestra
# ─────────────────────────────────────────────────────────────────────────────

def compute_qc_summary(
    sample_results: dict,
    gcms_data: dict,
) -> pd.DataFrame:
    """Calcula metricas QC por muestra (picos, SNR, area TIC, etc.)."""
    rows = []
    for name, sr in sample_results.items():
        pt = sr.peak_table
        num_peaks = len(pt)
        median_snr = float(pt["snr"].median()) if "snr" in pt.columns and not pt["snr"].isna().all() else None
        median_width = float(pt["width_sec"].median()) if "width_sec" in pt.columns and num_peaks > 0 else None

        n_coeluted = sum(1 for sp in sr.spectra if sp.is_coeluted) if hasattr(sr, "spectra") else 0
        pct_coeluted = (n_coeluted / num_peaks * 100) if num_peaks > 0 else 0.0

        tic_area = None
        if name in gcms_data:
            data = gcms_data[name]
            tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw
            dt = np.median(np.diff(data.scan_times))
            tic_area = float(np.sum(tic) * dt)

        rt_range = float(pt["rt_min"].max() - pt["rt_min"].min()) if num_peaks > 0 else 0.0

        rows.append({
            "sample": name,
            "num_peaks": num_peaks,
            "median_snr": round(median_snr, 2) if median_snr is not None else None,
            "median_width_sec": round(median_width, 2) if median_width is not None else None,
            "pct_coeluted": round(pct_coeluted, 1),
            "tic_area": round(tic_area, 0) if tic_area is not None else None,
            "rt_range_min": round(rt_range, 2),
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Utilidades
# ─────────────────────────────────────────────────────────────────────────────

def validate_minimum_data(
    n_samples: int,
    n_features: int,
    min_samples: int = 3,
    min_features: int = 2,
    name: str = "analisis",
) -> tuple:
    """
    Comprueba si hay suficientes datos para un analisis.

    Returns
    -------
    (True, '')         si OK.
    (False, mensaje)   si faltan muestras o features.
    """
    if n_samples < min_samples:
        msg = (f"No hay suficientes muestras para {name}: "
               f"{n_samples} disponibles, se requieren >= {min_samples}.")
        logger.warning("  %s", msg)
        return False, msg
    if n_features < min_features:
        msg = (f"No hay suficientes features para {name}: "
               f"{n_features} disponibles, se requieren >= {min_features}.")
        logger.warning("  %s", msg)
        return False, msg
    return True, ""


# ─────────────────────────────────────────────────────────────────────────────
# Metadata de muestras
# ─────────────────────────────────────────────────────────────────────────────

def infer_sample_metadata(
    sample_names: list,
) -> pd.DataFrame:
    """
    Infiere metadatos (tipo de muestra, grupo) desde nombres de archivo.

    Heuristica:
    - "qc", "quality"          -> sample_type = "QC"
    - "blank", "blanco", "blk" -> sample_type = "Blank"
    - "control"                -> sample_type = "Control"
    - resto                    -> sample_type = "Sample"
    """
    rows = []
    for name in sample_names:
        name_lower = name.lower()
        if any(k in name_lower for k in ("qc", "quality")):
            sample_type = "QC"
        elif any(k in name_lower for k in ("blank", "blanco", "blk")):
            sample_type = "Blank"
        elif "control" in name_lower:
            sample_type = "Control"
        else:
            sample_type = "Sample"

        rows.append({
            "sample_id": name,
            "sample_type": sample_type,
            "group": "GROUP_A",
            "order": len(rows) + 1,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Preprocesado estadistico de features
# ─────────────────────────────────────────────────────────────────────────────

def preprocess_feature_matrix(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    impute_method: str = config.IMPUTE_METHOD,
    normalization: str = config.NORMALIZATION_METHOD,
    log_transform: str = config.LOG_TRANSFORM,
    scaling: str = config.SCALING_METHOD,
    max_missing_frac: float = config.FEATURE_MAX_MISSING_FRAC,
    enable_blank_filter: bool = config.ENABLE_BLANK_FILTER,
    blank_ratio_threshold: float = config.BLANK_RATIO_THRESHOLD,
    enable_qc_rsd_filter: bool = config.ENABLE_QC_RSD_FILTER,
    qc_rsd_threshold: float = config.QC_RSD_THRESHOLD,
) -> tuple:
    """
    Preprocesa la matriz de features con pipeline configurable.

    Pasos (en orden):
    1. Blank filtering     (solo si hay blanks en sample_metadata)
    2. Missing value filter por feature
    3. Imputacion
    4. Normalizacion
    5. Log transform
    6. Scaling
    7. QC RSD filter       (solo si hay QCs en sample_metadata)

    Parameters
    ----------
    matrix          : DataFrame (muestras x features) con posibles NaN
    sample_metadata : DataFrame con columnas sample_id y sample_type.
                      Si es None, se omiten filtros de blank y QC RSD.

    Returns
    -------
    (matrix_processed, kept_cols)
        matrix_processed : DataFrame preprocesado
        kept_cols        : lista con los nombres de columnas que sobrevivieron
    """
    X = matrix.copy()
    n_initial = X.shape[1]
    logger.info("  Preprocesado: %d muestras x %d features", X.shape[0], n_initial)

    # ── PASO 1: Blank filtering ───────────────────────────────────────────────
    if enable_blank_filter and sample_metadata is not None:
        blank_ids = sample_metadata.loc[
            sample_metadata["sample_type"] == "Blank", "sample_id"
        ].tolist()
        blank_samples = [b for b in blank_ids if b in X.index]
        if blank_samples:
            logger.info("    Blank filter: %d blanks detectados", len(blank_samples))
            X_blank = X.loc[blank_samples].fillna(0)
            X_sample = X.loc[~X.index.isin(blank_samples)].fillna(0)
            ratio = (X_sample.mean() / (X_blank.mean() + 1e-6)).fillna(0)
            keep = ratio >= blank_ratio_threshold
            X = X.loc[:, keep]
            logger.info("      Retenidas: %d features (ratio >= %.1f)",
                        int(keep.sum()), blank_ratio_threshold)

    # ── PASO 2: Missing value filter ──────────────────────────────────────────
    missing_frac = X.isna().sum() / len(X)
    keep = missing_frac <= max_missing_frac
    X = X.loc[:, keep]
    logger.info("    Missing filter (max_frac=%.2f): %d features retenidas",
                max_missing_frac, X.shape[1])

    if X.shape[1] == 0:
        logger.warning("    Todas las features eliminadas por missing filter.")
        return X, []

    # ── PASO 3: Imputacion ────────────────────────────────────────────────────
    logger.info("    Imputacion: %s", impute_method)
    if impute_method == "min_half":
        for col in X.columns:
            col_vals = X[col].dropna()
            positive = col_vals[col_vals > 0]
            min_val = float(positive.min()) if len(positive) > 0 else float(col_vals.min() if len(col_vals) > 0 else 0)
            X.loc[X[col].isna(), col] = min_val / 2.0
    elif impute_method == "median":
        X = X.fillna(X.median())
    # "none": dejar NaN

    # ── PASO 4: Normalizacion ─────────────────────────────────────────────────
    logger.info("    Normalizacion: %s", normalization)
    if normalization in ("tic", "total_sum"):
        row_sums = X.sum(axis=1).replace(0, 1e-10)
        X = X.div(row_sums, axis=0)

    # ── PASO 5: Log transform ─────────────────────────────────────────────────
    logger.info("    Log transform: %s", log_transform)
    X = X.fillna(0)
    if log_transform == "log10":
        X = np.log10(X + 1.0)
    elif log_transform == "log2":
        X = np.log2(X + 1.0)

    # ── PASO 6: Scaling ───────────────────────────────────────────────────────
    logger.info("    Scaling: %s", scaling)
    if scaling == "autoscale":
        std = X.std(axis=0).replace(0, 1e-10)
        X = (X - X.mean(axis=0)) / std
    elif scaling == "pareto":
        std = X.std(axis=0).replace(0, 1e-10)
        X = (X - X.mean(axis=0)) / np.sqrt(std)

    # ── PASO 7: QC RSD filter ─────────────────────────────────────────────────
    if enable_qc_rsd_filter and sample_metadata is not None:
        qc_ids = sample_metadata.loc[
            sample_metadata["sample_type"] == "QC", "sample_id"
        ].tolist()
        qc_samples = [q for q in qc_ids if q in X.index]
        if len(qc_samples) >= 2:
            X_qc = X.loc[qc_samples]
            rsd = (X_qc.std() / (X_qc.mean().abs() + 1e-10)) * 100
            keep = rsd <= qc_rsd_threshold
            X = X.loc[:, keep]
            logger.info("    QC RSD filter (<=%.1f%%): %d features retenidas",
                        qc_rsd_threshold, X.shape[1])

    kept_cols = list(X.columns)
    logger.info("  Preprocesado completado: %d -> %d features", n_initial, len(kept_cols))
    return X, kept_cols


# ─────────────────────────────────────────────────────────────────────────────
# PCA de features (analisis principal)
# ─────────────────────────────────────────────────────────────────────────────

def pca_features(
    matrix: pd.DataFrame,
    n_components: int = config.PCA_N_COMPONENTS,
) -> dict | None:
    """
    PCA sobre matriz de features preprocesada.

    Returns
    -------
    Dict con scores, loadings, explained_variance_ratio, sample_names, feature_names.
    None si hay menos de 2 muestras o 0 features.
    """
    n_samples, n_features = matrix.shape

    if n_samples < 2:
        logger.warning("  PCA features requiere >= 2 muestras (hay %d).", n_samples)
        return None
    if n_features < 1:
        logger.warning("  PCA features requiere >= 1 feature (hay %d).", n_features)
        return None

    X = matrix.values.astype(np.float64)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    n_comp = min(n_components, n_samples - 1, n_features)

    logger.info("  PCA features: %d x %d -> %d componentes", n_samples, n_features, n_comp)

    try:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_comp)
        scores = pca.fit_transform(X)
        result = {
            "scores": scores,
            "loadings": pca.components_,
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    except ImportError:
        logger.info("  sklearn no disponible, usando numpy SVD.")
        result = _pca_numpy(X, n_comp)

    result["sample_names"] = list(matrix.index)
    result["feature_names"] = list(matrix.columns)
    logger.info("  PCA features: varianza explicada = %s",
                [f"{v:.1%}" for v in result["explained_variance_ratio"]])
    return result


def plot_pca_features(
    pca_result: dict | None,
    sample_metadata: pd.DataFrame | None = None,
    save_path: str | Path | None = None,
) -> plt.Figure | None:
    """
    Scatter PC1 vs PC2 de features, coloreado por tipo de muestra.
    """
    if pca_result is None:
        return None

    scores = pca_result["scores"]
    names = pca_result["sample_names"]
    evr = pca_result["explained_variance_ratio"]

    color_map = {
        "Sample": "steelblue", "QC": "orange",
        "Blank": "red", "Control": "green", "Unknown": "gray",
    }

    colors = ["steelblue"] * len(names)
    unique_types: set = set()
    if sample_metadata is not None:
        meta_dict = dict(zip(sample_metadata["sample_id"], sample_metadata["sample_type"]))
        colors = []
        for n in names:
            st = meta_dict.get(n, "Unknown")
            unique_types.add(st)
            colors.append(color_map.get(st, "gray"))
    else:
        unique_types = {"Sample"}

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(scores[:, 0], scores[:, 1],
               c=colors, s=80, edgecolor="white", alpha=0.85, zorder=3)

    if getattr(config, "PCA_LABEL_POINTS", False) or len(names) <= 12:
        for i, name in enumerate(names):
            ax.annotate(name, (scores[i, 0], scores[i, 1]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.75)

    ax.set_xlabel(f"PC1 ({evr[0]:.1%} varianza explicada)", fontsize=11)
    ax.set_ylabel(f"PC2 ({evr[1]:.1%} varianza explicada)" if len(evr) > 1 else "PC2", fontsize=11)
    ax.set_title("PCA — Features alineadas (analisis quimico principal)", fontsize=12)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.grid(True, alpha=0.3)

    if unique_types:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color_map.get(t, "gray"), edgecolor="white", label=t)
            for t in sorted(unique_types)
        ]
        ax.legend(handles=legend_elements, loc="best", fontsize=9, framealpha=0.8)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def compute_pca_loadings_table(
    pca_result: dict | None,
    feat_meta: pd.DataFrame | None = None,
    top_n: int = 20,
) -> pd.DataFrame:
    """
    Genera tabla de top loadings para PC1 y PC2.

    Si feat_meta esta disponible (columnas: compound_name, mean_rt_min,
    best_match_score, formula), anade esa informacion a la tabla.

    Returns
    -------
    DataFrame con: PC, rank, feature, loading, abs_loading,
    y opcionalmente: mean_rt_min, best_match_score, formula.
    """
    if pca_result is None:
        return pd.DataFrame()

    loadings = pca_result["loadings"]
    feature_names = pca_result.get("feature_names", [])

    rows = []
    for pc_idx in range(min(2, len(loadings))):
        pc_name = f"PC{pc_idx + 1}"
        load_vals = loadings[pc_idx]
        top_idx = np.argsort(np.abs(load_vals))[::-1][:top_n]

        for rank, idx in enumerate(top_idx):
            fname = feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
            row: dict = {
                "PC": pc_name,
                "rank": rank + 1,
                "feature": fname,
                "loading": round(float(load_vals[idx]), 6),
                "abs_loading": round(float(abs(load_vals[idx])), 6),
            }
            if feat_meta is not None and "compound_name" in feat_meta.columns:
                match = feat_meta[feat_meta["compound_name"] == fname]
                if not match.empty:
                    row["mean_rt_min"] = match.iloc[0].get("mean_rt_min", "")
                    row["best_match_score"] = match.iloc[0].get("best_match_score", "")
                    row["formula"] = match.iloc[0].get("formula", "")
            rows.append(row)

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Heatmaps
# ─────────────────────────────────────────────────────────────────────────────

def _select_top_features(matrix: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Selecciona las top_n features mas variables (por varianza)."""
    variances = matrix.var(axis=0)
    top_features = variances.nlargest(min(top_n, matrix.shape[1])).index
    return matrix[top_features]


def _order_samples_by_type(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
) -> pd.DataFrame:
    """Reordena filas (muestras) por tipo de muestra si hay metadata."""
    if sample_metadata is None:
        return matrix
    meta_dict = dict(zip(sample_metadata["sample_id"], sample_metadata["sample_type"]))
    type_order = {"QC": 0, "Blank": 1, "Control": 2, "Sample": 3, "Unknown": 4}
    ordered = sorted(
        matrix.index,
        key=lambda x: (type_order.get(meta_dict.get(x, "Unknown"), 9), x),
    )
    return matrix.loc[ordered]


def _add_group_separators(
    ax: plt.Axes,
    col_names: list,
    sample_metadata: pd.DataFrame | None,
) -> None:
    """Anade lineas verticales entre grupos de muestras distintos."""
    if sample_metadata is None:
        return
    meta_dict = dict(zip(sample_metadata["sample_id"], sample_metadata["sample_type"]))
    prev_type = None
    for i, name in enumerate(col_names):
        curr_type = meta_dict.get(name, "Unknown")
        if prev_type is not None and curr_type != prev_type:
            ax.axvline(i - 0.5, color="black", linewidth=1.5, zorder=4)
        prev_type = curr_type


def plot_heatmap_main(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    top_n: int = 50,
    save_path: str | Path | None = None,
) -> plt.Figure | None:
    """
    Heatmap principal tipo 'paper' (sin z-score).

    - Filas    = metabolitos/features (top N por varianza)
    - Columnas = muestras, ordenadas por tipo de muestra
    - Valor    = abundancia log-normalizada (tal como llega en matrix)
    - Separadores visuales entre grupos de muestras distintos

    Parameters
    ----------
    matrix          : DataFrame preprocesado (muestras x features)
    sample_metadata : DataFrame con sample_id y sample_type
    top_n           : cuantos features mostrar
    save_path       : ruta para guardar
    """
    if matrix.shape[0] < 2 or matrix.shape[1] < 2:
        logger.warning("  Heatmap principal: datos insuficientes.")
        return None

    X_top = _select_top_features(matrix, top_n)
    X_ordered = _order_samples_by_type(X_top, sample_metadata)

    # Clustering de features (columnas de X_ordered = filas del heatmap)
    use_clustering = getattr(config, "HEATMAP_CLUSTERING", True)
    if use_clustering:
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import pdist
            col_link = linkage(
                pdist(X_ordered.T.fillna(0).values, metric="euclidean"), method="ward"
            )
            col_order = dendrogram(col_link, no_plot=True)["leaves"]
            X_ordered = X_ordered.iloc[:, col_order]
        except Exception:
            use_clustering = False

    # Transponer: filas=features, columnas=muestras
    X_plot = X_ordered.T

    n_feat = X_plot.shape[0]
    n_samp = X_plot.shape[1]
    fig_w = max(8, n_samp * 1.1)
    fig_h = max(5, n_feat * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vals = X_plot.values.astype(float)
    # Usar percentiles para evitar que outliers dominen la escala de color
    vmin = np.nanpercentile(vals, 2)
    vmax = np.nanpercentile(vals, 98)
    im = ax.imshow(vals, aspect="auto", cmap="YlOrRd",
                   vmin=vmin, vmax=vmax, interpolation="nearest")

    ax.set_xticks(range(n_samp))
    ax.set_yticks(range(n_feat))
    ax.set_xticklabels(list(X_plot.columns), rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(list(X_plot.index), fontsize=8)
    ax.set_xlabel("Muestras", fontsize=10)
    ax.set_ylabel(f"Features (top {n_feat} por varianza)", fontsize=10)
    clust_label = " + clustering" if use_clustering else ""
    ax.set_title(
        f"Heatmap principal — Features alineadas (log-normalizado{clust_label})", fontsize=11
    )

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Abundancia (log)", fontsize=9)

    _add_group_separators(ax, list(X_plot.columns), sample_metadata)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_heatmap_features(
    matrix: pd.DataFrame,
    top_n: int = 50,
    sample_metadata: pd.DataFrame | None = None,
    save_path: str | Path | None = None,
) -> plt.Figure | None:
    """
    Heatmap exploratorio con z-score por feature.

    - Filas    = metabolitos/features (top N por varianza)
    - Columnas = muestras, ordenadas por tipo de muestra
    - Valor    = z-score por feature (para visualizar patrones relativos)
    - Clustering jerarquico si scipy disponible

    Parameters
    ----------
    matrix          : DataFrame preprocesado (muestras x features)
    top_n           : cuantos features mostrar
    sample_metadata : DataFrame con sample_id y sample_type
    save_path       : ruta para guardar
    """
    if matrix.shape[0] < 2 or matrix.shape[1] < 2:
        logger.warning("  Heatmap z-score: datos insuficientes.")
        return None

    X_top = _select_top_features(matrix, top_n)
    X_ordered = _order_samples_by_type(X_top, sample_metadata)

    # Z-score por columna (feature)
    std = X_ordered.std(axis=0).replace(0, 1e-10)
    X_zscore = (X_ordered - X_ordered.mean(axis=0)) / std

    use_clustering = getattr(config, "HEATMAP_CLUSTERING", True)
    if use_clustering:
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import pdist
            col_link = linkage(
                pdist(X_zscore.T.fillna(0).values, metric="euclidean"), method="ward"
            )
            col_order = dendrogram(col_link, no_plot=True)["leaves"]
            X_zscore = X_zscore.iloc[:, col_order]
        except Exception:
            use_clustering = False

    # Transponer: filas=features, columnas=muestras
    X_plot = X_zscore.T

    n_feat = X_plot.shape[0]
    n_samp = X_plot.shape[1]
    fig_w = max(8, n_samp * 1.1)
    fig_h = max(5, n_feat * 0.28)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(X_plot.values.astype(float), aspect="auto",
                   cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")

    ax.set_xticks(range(n_samp))
    ax.set_yticks(range(n_feat))
    ax.set_xticklabels(list(X_plot.columns), rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(list(X_plot.index), fontsize=8)
    ax.set_xlabel("Muestras", fontsize=10)
    ax.set_ylabel(f"Features (top {n_feat} por varianza)", fontsize=10)
    clust_label = " + clustering" if use_clustering else ""
    ax.set_title(f"Heatmap exploratorio — z-score{clust_label}", fontsize=11)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("z-score", fontsize=9)

    _add_group_separators(ax, list(X_plot.columns), sample_metadata)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# QC summary avanzado (antes / despues preproceso)
# ─────────────────────────────────────────────────────────────────────────────

def compute_advanced_qc_summary(
    matrix_raw: pd.DataFrame,
    matrix_processed: pd.DataFrame,
    sample_metadata: pd.DataFrame | None = None,
    preprocessing_stats: dict | None = None,
) -> pd.DataFrame:
    """
    Resumen QC por muestra: features detectadas, missing, intensidad total.

    Parameters
    ----------
    matrix_raw        : DataFrame crudo (muestras x features)
    matrix_processed  : DataFrame preprocesado
    sample_metadata   : DataFrame con sample_id y sample_type (opcional)
    preprocessing_stats : dict con estadisticas del preprocesado (opcional, ignorado)
    """
    rows = []
    meta_dict: dict = {}
    if sample_metadata is not None:
        meta_dict = dict(zip(sample_metadata["sample_id"], sample_metadata["sample_type"]))

    for sample in matrix_raw.index:
        raw_row = matrix_raw.loc[sample]
        proc_row = (
            matrix_processed.loc[sample]
            if sample in matrix_processed.index
            else pd.Series(dtype=float)
        )

        rows.append({
            "sample_id": sample,
            "sample_type": meta_dict.get(sample, "Unknown"),
            "n_features_raw": int(raw_row.notna().sum()),
            "n_features_processed": int(proc_row.notna().sum()) if len(proc_row) > 0 else 0,
            "pct_missing_raw": round(float(raw_row.isna().mean() * 100), 2),
            "sum_intensity_raw": round(float(raw_row.fillna(0).sum()), 2),
            "sum_intensity_processed": round(float(proc_row.fillna(0).sum()), 4) if len(proc_row) > 0 else 0.0,
        })

    df = pd.DataFrame(rows)

    # Fila de totales globales
    df = pd.concat([df, pd.DataFrame([{
        "sample_id": "GLOBAL_SUMMARY",
        "sample_type": "—",
        "n_features_raw": matrix_raw.shape[1],
        "n_features_processed": matrix_processed.shape[1],
        "pct_missing_raw": round(float(matrix_raw.isna().mean().mean() * 100), 2),
        "sum_intensity_raw": round(float(matrix_raw.fillna(0).values.sum()), 2),
        "sum_intensity_processed": round(float(matrix_processed.fillna(0).values.sum()), 4),
    }])], ignore_index=True)

    return df
