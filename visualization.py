"""
Modulo 6 -- Visualizacion y reportes.

Funciones de plotting para QC del pipeline.
Incluye visualizaciones single-sample y multi-sample.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numpy.typing import NDArray

import config
from preprocessing import GCMSData
from peak_detection import PeakResult


def plot_tic_preprocessing(
    data: GCMSData,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """
    Visualiza el TIC antes y después del preprocesado (3 paneles).

    Panel 1: TIC raw
    Panel 2: TIC suavizado
    Panel 3: TIC corregido (baseline-subtracted)
    """
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(f"Preprocesado TIC — {data.name}", fontsize=14, fontweight="bold")
    t = data.scan_times_min

    # Panel 1 - Raw
    axes[0].plot(t, data.tic_raw, color="gray", linewidth=0.5)
    axes[0].set_ylabel("Intensidad")
    axes[0].set_title("TIC Original")

    # Panel 2 - Smoothed
    if data.tic_smooth is not None:
        axes[1].plot(t, data.tic_raw, color="lightgray", linewidth=0.3, label="Raw")
        axes[1].plot(t, data.tic_smooth, color="blue", linewidth=0.7, label="Suavizado")
        axes[1].legend(loc="upper right", fontsize=8)
    axes[1].set_ylabel("Intensidad")
    axes[1].set_title(f"Savitzky-Golay (window={config.SAVGOL_WINDOW}, order={config.SAVGOL_POLYORDER})")

    # Panel 3 - Baseline corrected
    if data.tic_corrected is not None:
        axes[2].plot(t, data.tic_corrected, color="darkgreen", linewidth=0.7)
    axes[2].set_ylabel("Intensidad")
    axes[2].set_xlabel("Tiempo de retención (min)")
    axes[2].set_title("Corrección de baseline (AsLS)")

    for ax in axes:
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_peaks(
    data: GCMSData,
    peak_result: PeakResult,
    save_path: str | Path | None = None,
    annotate_top_n: int = 20,
) -> plt.Figure:
    """
    Visualiza el TIC corregido con los picos detectados marcados.
    Anota los top N picos más intensos con su RT.
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw
    t = data.scan_times_min
    peaks = peak_result.indices

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.plot(t, tic, color="darkgreen", linewidth=0.5, label="TIC corregido")
    ax.plot(
        t[peaks], tic[peaks],
        "rv", markersize=5, alpha=0.7,
        label=f"Picos detectados (n={len(peaks)})",
    )

    # Anotar los picos más intensos
    if len(peaks) > 0:
        sorted_by_intensity = np.argsort(tic[peaks])[::-1]
        for rank, idx in enumerate(sorted_by_intensity[:annotate_top_n]):
            peak_idx = peaks[idx]
            ax.annotate(
                f"{t[peak_idx]:.2f}",
                xy=(t[peak_idx], tic[peak_idx]),
                xytext=(0, 10),
                textcoords="offset points",
                fontsize=6,
                ha="center",
                color="red",
                rotation=45,
            )

    ax.set_xlabel("Tiempo de retención (min)")
    ax.set_ylabel("Intensidad")
    ax.set_title(f"Detección de picos — {data.name} ({peak_result.n_peaks} picos)")
    ax.legend(loc="upper right")
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_peak_quality(
    peak_result: PeakResult,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """
    Gráficos de calidad de los picos detectados:
    - Histograma de intensidades
    - Histograma de S/N
    - Histograma de anchos
    - Scatter RT vs Intensidad
    """
    df = peak_result.table
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Control de calidad — Picos detectados", fontsize=13, fontweight="bold")

    # Intensidades
    axes[0, 0].hist(df["intensity"], bins=50, color="steelblue", edgecolor="white")
    axes[0, 0].set_xlabel("Intensidad")
    axes[0, 0].set_ylabel("Frecuencia")
    axes[0, 0].set_title("Distribución de intensidades")

    # S/N
    if "snr" in df.columns and df["snr"].notna().any():
        axes[0, 1].hist(df["snr"].dropna(), bins=50, color="coral", edgecolor="white")
        axes[0, 1].axvline(config.PEAK_SNR_THRESHOLD, color="red", linestyle="--", label=f"Umbral={config.PEAK_SNR_THRESHOLD}")
        axes[0, 1].legend()
    axes[0, 1].set_xlabel("S/N Ratio")
    axes[0, 1].set_ylabel("Frecuencia")
    axes[0, 1].set_title("Distribución de S/N")

    # Anchos
    axes[1, 0].hist(df["width_sec"], bins=50, color="mediumpurple", edgecolor="white")
    axes[1, 0].set_xlabel("Ancho del pico (s)")
    axes[1, 0].set_ylabel("Frecuencia")
    axes[1, 0].set_title("Distribución de anchos")

    # RT vs Intensidad
    axes[1, 1].scatter(df["rt_min"], df["intensity"], s=8, alpha=0.6, color="darkgreen")
    axes[1, 1].set_xlabel("Tiempo de retención (min)")
    axes[1, 1].set_ylabel("Intensidad")
    axes[1, 1].set_title("RT vs Intensidad")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# =====================================================================
# Visualizaciones multi-muestra (Modulo 6 ampliado)
# =====================================================================


def plot_tic_overlay(
    samples: dict[str, GCMSData],
    save_path: str | Path | None = None,
    normalize: bool = False,
) -> plt.Figure:
    """
    Overlay de TICs corregidos de todas las muestras.

    Parameters
    ----------
    normalize : si True, normaliza cada TIC por su area total (cada curva suma 1).
    """
    fig, ax = plt.subplots(figsize=(16, 6))

    cmap = plt.colormaps["tab20"].resampled(max(len(samples), 1))
    for i, (name, data) in enumerate(samples.items()):
        tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw
        if normalize:
            total = np.sum(tic)
            tic = tic / total if total > 0 else tic
        ax.plot(data.scan_times_min, tic, linewidth=0.5, color=cmap(i), label=name, alpha=0.7)

    ax.set_xlabel("Tiempo de retencion (min)")
    ylabel = "Proporcion (normalizado)" if normalize else "Intensidad"
    ax.set_ylabel(ylabel)
    title_suffix = " (normalizado por area)" if normalize else ""
    ax.set_title(f"TIC Overlay -- {len(samples)} muestras{title_suffix}")
    if not normalize:
        ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    if len(samples) <= 20:
        ax.legend(loc="upper right", fontsize=6, ncol=2)
    else:
        ax.legend(loc="upper right", fontsize=5, ncol=3)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_match_score_histogram(
    identification_table: pd.DataFrame,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """
    Histograma de match scores para evaluar calidad de identificaciones.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Calidad de las identificaciones", fontsize=13, fontweight="bold")

    scores = identification_table["match_score_nist"]
    identified = scores[scores > 0]

    axes[0].hist(identified, bins=40, color="steelblue", edgecolor="white", alpha=0.8)
    threshold_nist = config.NIST_MATCH_THRESHOLD
    axes[0].axvline(threshold_nist, color="red", linestyle="--", linewidth=1.5,
                    label=f"Umbral={threshold_nist}")
    axes[0].set_xlabel("Match Score (0-999)")
    axes[0].set_ylabel("Frecuencia")
    axes[0].set_title(f"Distribucion de Match Scores (n={len(identified)})")
    axes[0].legend()

    n_total = len(identification_table)
    n_id = (identification_table["compound_name"] != "No match").sum()
    n_no_id = n_total - n_id
    axes[1].pie(
        [n_id, n_no_id],
        labels=[f"Identificados\n(n={n_id})", f"No identificados\n(n={n_no_id})"],
        colors=["#2ecc71", "#e74c3c"],
        autopct="%1.1f%%",
        startangle=90,
    )
    axes[1].set_title("Tasa de identificacion")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_data_matrix_heatmap(
    data_matrix: pd.DataFrame,
    save_path: str | Path | None = None,
    max_features: int = 80,
) -> plt.Figure:
    """
    Heatmap de la matriz de datos (muestras x metabolitos).
    """
    if data_matrix.shape[1] > max_features:
        n_present = data_matrix.notna().sum().sort_values(ascending=False)
        top_cols = n_present.head(max_features).index
        matrix_plot = data_matrix[top_cols]
        title_suffix = f" (top {max_features} features por frecuencia)"
    else:
        matrix_plot = data_matrix
        title_suffix = ""

    matrix_log = np.log10(matrix_plot.fillna(0) + 1)

    fig_height = max(4, len(matrix_plot) * 0.4)
    fig_width = max(10, len(matrix_plot.columns) * 0.15)
    fig, ax = plt.subplots(figsize=(min(fig_width, 30), min(fig_height, 15)))

    sns.heatmap(
        matrix_log,
        ax=ax,
        cmap="YlOrRd",
        xticklabels=True,
        yticklabels=True,
        linewidths=0.2,
        linecolor="white",
        cbar_kws={"label": "log10(Intensidad + 1)"},
    )

    ax.set_title(f"Matriz de datos{title_suffix}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Metabolitos")
    ax.set_ylabel("Muestras")
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=8)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_sample_summary(
    sample_stats: pd.DataFrame,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """
    Resumen QC por muestra: numero de picos, identificaciones, intensidades.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle("Resumen de QC por muestra", fontsize=13, fontweight="bold")

    n = len(sample_stats)
    x = range(n)
    labels = sample_stats["sample"].tolist()

    axes[0, 0].bar(x, sample_stats["n_peaks"], color="steelblue", edgecolor="white")
    axes[0, 0].set_ylabel("N picos")
    axes[0, 0].set_title("Picos detectados por muestra")
    axes[0, 0].set_xticks(list(x))
    axes[0, 0].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)

    if "n_identified" in sample_stats.columns:
        axes[0, 1].bar(x, sample_stats["n_identified"], color="#2ecc71", edgecolor="white",
                       label="Identificados")
        n_not_id = sample_stats["n_peaks"] - sample_stats["n_identified"]
        axes[0, 1].bar(x, n_not_id, bottom=sample_stats["n_identified"],
                       color="#e74c3c", edgecolor="white", label="No identificados")
        axes[0, 1].legend(fontsize=8)
    axes[0, 1].set_ylabel("N picos")
    axes[0, 1].set_title("Identificaciones por muestra")
    axes[0, 1].set_xticks(list(x))
    axes[0, 1].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)

    axes[1, 0].bar(x, sample_stats["mean_intensity"], color="coral", edgecolor="white")
    axes[1, 0].set_ylabel("Intensidad media")
    axes[1, 0].set_title("Intensidad media de picos")
    axes[1, 0].set_xticks(list(x))
    axes[1, 0].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    axes[1, 0].ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    if "pct_missing" in sample_stats.columns:
        axes[1, 1].bar(x, sample_stats["pct_missing"], color="mediumpurple", edgecolor="white")
        axes[1, 1].set_ylabel("% Missing")
        axes[1, 1].set_title("Features no detectados (%)")
        axes[1, 1].set_xticks(list(x))
        axes[1, 1].set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    else:
        axes[1, 1].text(0.5, 0.5, "Sin datos de\nmissing values",
                        ha="center", va="center", transform=axes[1, 1].transAxes)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_spectrum(
    mz: NDArray,
    intensities: NDArray,
    title: str = "Espectro de masas",
    top_n_annotate: int = 8,
    save_path: str | Path | None = None,
) -> plt.Figure:
    """
    Visualiza un espectro de masas individual (barras verticales).
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.vlines(mz, 0, intensities, color="navy", linewidth=0.8)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensidad")
    ax.set_title(title)

    if len(intensities) > 0:
        top_idx = intensities.argsort()[-top_n_annotate:]
        for idx in top_idx:
            ax.annotate(
                f"{mz[idx]:.0f}",
                xy=(mz[idx], intensities[idx]),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center", fontsize=7, color="red",
            )

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig
