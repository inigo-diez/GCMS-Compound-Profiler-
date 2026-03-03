"""
Módulo 2 — Detección automática de picos sobre el TIC preprocesado.

Funciones:
    - detect_peaks: detección adaptativa con scipy.signal.find_peaks.
    - filter_peaks_snr: filtra picos por ratio señal/ruido dinámico local.
    - build_peak_table: genera un DataFrame con información de cada pico.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.signal import find_peaks, peak_widths

import config
from preprocessing import GCMSData

logger = logging.getLogger(__name__)


@dataclass
class PeakResult:
    """Resultado de la detección de picos."""

    indices: NDArray[np.intp]
    properties: dict
    table: pd.DataFrame
    n_peaks: int = 0


def _local_percentile_baseline(
    tic: NDArray[np.float64],
    window: int = 41,
    q: float = 20.0,
) -> NDArray[np.float64]:
    """Baseline local robusta usando percentil en ventana deslizante."""
    n = len(tic)
    half = max(1, window // 2)
    baseline = np.zeros(n, dtype=np.float64)

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        baseline[i] = np.percentile(tic[lo:hi], q)

    return baseline


def _estimate_noise(
    tic: NDArray[np.float64],
    window: int = 41,
) -> NDArray[np.float64]:
    """Estimación de ruido local (MAD robusto) por ventana deslizante."""
    n = len(tic)
    half = max(1, window // 2)
    noise = np.zeros(n, dtype=np.float64)

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        seg = tic[lo:hi]
        med = np.median(seg)
        mad = np.median(np.abs(seg - med))
        sigma = 1.4826 * mad
        noise[i] = sigma

    min_floor = max(1e-8, float(np.std(tic)) * 1e-3)
    noise[noise < min_floor] = min_floor
    return noise


def _log_tic_diagnostics(tic: NDArray[np.float64]) -> None:
    """Logs diagnósticos del TIC para depurar problemas de no detección."""
    if len(tic) == 0:
        logger.warning("  TIC vacío: no se puede detectar picos.")
        return

    logger.info("  TIC stats: min=%.3f, max=%.3f, mean=%.3f, std=%.3f",
                float(np.min(tic)), float(np.max(tic)), float(np.mean(tic)), float(np.std(tic)))

    q = np.percentile(tic, [5, 25, 50, 75, 95, 99])
    logger.info(
        "  TIC percentiles: p05=%.3f, p25=%.3f, p50=%.3f, p75=%.3f, p95=%.3f, p99=%.3f",
        float(q[0]), float(q[1]), float(q[2]), float(q[3]), float(q[4]), float(q[5]),
    )


def _analyze_prominence_distribution(tic: NDArray[np.float64]) -> None:
    """Analiza distribución de prominencias sin umbral para diagnóstico."""
    quick_peaks, quick_props = find_peaks(tic, prominence=0)
    if len(quick_peaks) == 0:
        logger.warning("  Diagnóstico prominencias: no hay máximos locales en la señal.")
        return

    prominences = np.asarray(quick_props.get("prominences", []), dtype=np.float64)
    if len(prominences) == 0:
        logger.warning("  Diagnóstico prominencias: scipy no devolvió prominencias.")
        return

    qp = np.percentile(prominences, [25, 50, 75, 90, 95, 99])
    logger.info(
        "  Prominencias (%d máximos): p25=%.3f, p50=%.3f, p75=%.3f, p90=%.3f, p95=%.3f, p99=%.3f",
        len(prominences), float(qp[0]), float(qp[1]), float(qp[2]), float(qp[3]), float(qp[4]), float(qp[5]),
    )


def detect_peaks(
    tic: NDArray[np.float64],
    min_prominence: float = config.PEAK_MIN_PROMINENCE,
    min_distance: int = config.PEAK_MIN_DISTANCE,
    min_width: float = config.PEAK_MIN_WIDTH,
) -> tuple[NDArray[np.intp], dict]:
    """
    Detección adaptativa de picos sobre TIC preprocesado.

    Estrategia:
    - calcula umbrales robustos iniciales (prominence/height) a partir de la señal,
    - intenta detección en varios niveles (de estricto a laxo),
    - si no detecta, reduce progresivamente thresholds.
    """
    if len(tic) == 0:
        return np.array([], dtype=np.intp), {}

    _log_tic_diagnostics(tic)
    _analyze_prominence_distribution(tic)

    noise = _estimate_noise(tic, window=41)
    baseline = _local_percentile_baseline(tic, window=41, q=20)

    noise_global = float(np.median(noise))
    p95 = float(np.percentile(tic, 95))
    p75 = float(np.percentile(tic, 75))
    dyn_range = max(float(np.max(tic) - np.min(tic)), 1e-8)

    adaptive_prom = max(noise_global * 3.0, dyn_range * 0.01)
    base_height = max(p75 + noise_global * 2.0, p95 * 0.3)

    distance_steps = sorted({max(1, int(min_distance * s)) for s in (1.0, 0.75, 0.5)})
    width_steps = [max(1.0, float(min_width) * s) for s in (1.0, 0.75, 0.5)]

    # Prioriza valor de config, luego thresholds adaptativos cada vez más laxos
    prom_candidates = [float(min_prominence), adaptive_prom]
    prom_candidates.extend(adaptive_prom * s for s in (0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
    prom_candidates = sorted({max(noise_global * 0.3, p) for p in prom_candidates}, reverse=True)

    last_peaks = np.array([], dtype=np.intp)
    last_props: dict = {}

    for i, prom in enumerate(prom_candidates, start=1):
        # altura acompaña a la prominencia, pero nunca supera p95 para evitar bloqueo
        height = min(max(base_height * (prom / max(adaptive_prom, 1e-8)), p75), p95)

        for dist in distance_steps:
            for width in width_steps:
                peaks, props = find_peaks(
                    tic,
                    prominence=prom,
                    height=height,
                    distance=dist,
                    width=width,
                    rel_height=0.5,
                )

                logger.info(
                    "  Intento %d | prom>=%.3f, height>=%.3f, distance>=%d, width>=%.2f -> %d picos",
                    i, prom, height, dist, width, len(peaks),
                )

                last_peaks, last_props = peaks, props
                if len(peaks) > 0:
                    return peaks, props

    logger.warning("  No se detectaron picos tras búsqueda adaptativa.")
    return last_peaks, last_props


def filter_peaks_snr(
    tic: NDArray[np.float64],
    peaks: NDArray[np.intp],
    properties: dict,
    snr_threshold: float = config.PEAK_SNR_THRESHOLD,
) -> tuple[NDArray[np.intp], dict]:
    """
    Filtra picos por S/N dinámico basado en baseline y ruido local.
    """
    if len(peaks) == 0:
        return peaks, {**properties, "snr": np.array([], dtype=np.float64)}

    noise = _estimate_noise(tic, window=41)
    baseline = _local_percentile_baseline(tic, window=41, q=20)

    signal = tic[peaks] - baseline[peaks]
    signal[signal < 0] = 0.0
    snr = signal / noise[peaks]

    mask = snr >= snr_threshold

    # Relax adaptativo si elimina todo
    if not np.any(mask) and len(peaks) > 0:
        dynamic_threshold = max(1.5, float(np.percentile(snr, 75) * 0.7))
        logger.info(
            "  Relax S/N adaptativo: umbral inicial %.2f -> %.2f",
            snr_threshold, dynamic_threshold,
        )
        mask = snr >= dynamic_threshold

    filtered_peaks = peaks[mask]

    filtered_props = {}
    for key, val in properties.items():
        if isinstance(val, np.ndarray) and len(val) == len(peaks):
            filtered_props[key] = val[mask]
        else:
            filtered_props[key] = val

    filtered_props["snr"] = snr[mask]

    n_removed = len(peaks) - len(filtered_peaks)
    logger.info(
        "  Filtro S/N: eliminados %d de %d picos, quedan %d",
        n_removed, len(peaks), len(filtered_peaks),
    )
    return filtered_peaks, filtered_props


def build_peak_table(
    data: GCMSData,
    peaks: NDArray[np.intp],
    properties: dict,
) -> pd.DataFrame:
    """
    Construye un DataFrame con la información de cada pico detectado.

    Columns: peak_id, scan_idx, rt_sec, rt_min, intensity, prominence,
             width_scans, width_sec, snr
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw

    if len(peaks) == 0:
        return pd.DataFrame(columns=[
            "peak_id", "scan_idx", "rt_sec", "rt_min", "intensity", "prominence",
            "width_scans", "width_sec", "snr",
        ])

    widths_result = peak_widths(tic, peaks, rel_height=0.5)
    widths_scans = widths_result[0]

    dt = np.median(np.diff(data.scan_times)) if len(data.scan_times) > 1 else 0.0
    widths_sec = widths_scans * dt

    rows = []
    for i, peak_idx in enumerate(peaks):
        row = {
            "peak_id": i + 1,
            "scan_idx": int(peak_idx),
            "rt_sec": float(data.scan_times[peak_idx]),
            "rt_min": float(data.scan_times_min[peak_idx]),
            "intensity": float(tic[peak_idx]),
            "prominence": float(properties.get("prominences", np.full(len(peaks), np.nan))[i]),
            "width_scans": float(widths_scans[i]),
            "width_sec": float(widths_sec[i]),
            "snr": float(properties.get("snr", np.full(len(peaks), np.nan))[i]),
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    logger.info("  Tabla de picos construida: %d picos", len(df))
    return df


def find_and_filter_peaks(data: GCMSData) -> PeakResult:
    """
    Pipeline completo de detección de picos:
    1. Detección adaptativa con find_peaks
    2. Filtrado por S/N dinámico
    3. Construcción de tabla resumen
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw

    peaks, props = detect_peaks(tic)
    peaks, props = filter_peaks_snr(tic, peaks, props)
    table = build_peak_table(data, peaks, props)

    result = PeakResult(
        indices=peaks,
        properties=props,
        table=table,
        n_peaks=len(peaks),
    )
    logger.info("  Resultado final: %d picos detectados y filtrados.", result.n_peaks)
    return result
