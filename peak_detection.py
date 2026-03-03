"""
Módulo 2 — Detección automática de picos sobre el TIC preprocesado.

Funciones:
    - detect_peaks: detecta picos con scipy.signal.find_peaks.
    - filter_peaks_snr: filtra picos por ratio señal/ruido.
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

    indices: NDArray[np.intp]             # Índices de los picos en el array del TIC
    properties: dict                       # Propiedades de scipy.signal.find_peaks
    table: pd.DataFrame                    # Tabla resumen de picos
    n_peaks: int = 0


def _estimate_noise(tic: NDArray[np.float64], window: int = 50) -> NDArray[np.float64]:
    """
    Estima el nivel de ruido local usando la desviación estándar en ventanas
    deslizantes alrededor de cada punto.
    """
    n = len(tic)
    noise = np.zeros(n)
    half = window // 2
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half)
        noise[i] = np.std(tic[lo:hi])
    # Evitar divisiones por cero
    noise[noise < 1.0] = 1.0
    return noise


def detect_peaks(
    tic: NDArray[np.float64],
    min_prominence: float = config.PEAK_MIN_PROMINENCE,
    min_distance: int = config.PEAK_MIN_DISTANCE,
    min_width: float = config.PEAK_MIN_WIDTH,
) -> tuple[NDArray[np.intp], dict]:
    """
    Detecta picos en el TIC preprocesado.

    Parameters
    ----------
    tic : TIC corregido (baseline-subtracted)
    min_prominence : prominencia mínima
    min_distance : distancia mínima entre picos (scans)
    min_width : ancho mínimo del pico (scans)

    Returns
    -------
    (indices, properties) de scipy.signal.find_peaks.
    """
    peaks, props = find_peaks(
        tic,
        prominence=min_prominence,
        distance=min_distance,
        width=min_width,
        rel_height=0.5,
    )
    logger.info("  find_peaks detectó %d picos (prominence>%.0f, distance>%d, width>%.1f)",
                len(peaks), min_prominence, min_distance, min_width)
    return peaks, props


def filter_peaks_snr(
    tic: NDArray[np.float64],
    peaks: NDArray[np.intp],
    properties: dict,
    snr_threshold: float = config.PEAK_SNR_THRESHOLD,
) -> tuple[NDArray[np.intp], dict]:
    """
    Filtra picos por ratio señal/ruido.

    Returns
    -------
    (filtered_peaks, filtered_properties)
    """
    noise = _estimate_noise(tic)
    snr = tic[peaks] / noise[peaks]

    mask = snr >= snr_threshold
    filtered_peaks = peaks[mask]

    # Filtrar propiedades también
    filtered_props = {}
    for key, val in properties.items():
        if isinstance(val, np.ndarray) and len(val) == len(peaks):
            filtered_props[key] = val[mask]
        else:
            filtered_props[key] = val

    # Guardar SNR en las propiedades
    filtered_props["snr"] = snr[mask]

    n_removed = len(peaks) - len(filtered_peaks)
    if n_removed > 0:
        logger.info("  Filtro S/N (>%.1f): eliminados %d picos, quedan %d",
                     snr_threshold, n_removed, len(filtered_peaks))
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

    # Calcular anchos con más precisión
    widths_result = peak_widths(tic, peaks, rel_height=0.5)
    widths_scans = widths_result[0]

    # Calcular ancho en segundos
    dt = np.median(np.diff(data.scan_times))  # tiempo entre scans
    widths_sec = widths_scans * dt

    rows = []
    for i, peak_idx in enumerate(peaks):
        row = {
            "peak_id": 0,  # se asigna después de ordenar
            "scan_idx": int(peak_idx),
            "rt_sec": data.scan_times[peak_idx],
            "rt_min": data.scan_times_min[peak_idx],
            "intensity": tic[peak_idx],
            "prominence": properties.get("prominences", [None])[i] if "prominences" in properties else None,
            "width_scans": widths_scans[i],
            "width_sec": widths_sec[i],
            "snr": properties.get("snr", [None])[i] if "snr" in properties else None,
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # Ordenar por intensidad (desc) y asignar peak_id como ranking
    if not df.empty:
        df = df.sort_values("intensity", ascending=False).reset_index(drop=True)
        df["peak_id"] = range(1, len(df) + 1)

    logger.info("  Tabla de picos construida: %d picos (ordenados por intensidad)", len(df))
    return df


def find_and_filter_peaks(
    data: GCMSData,
    max_peaks: int | None = config.MAX_PEAKS,
) -> PeakResult:
    """
    Pipeline completo de detección de picos:
    1. Detección con find_peaks
    2. Filtrado por S/N
    3. Construcción de tabla resumen (ordenada por intensidad desc)
    4. Recorte opcional a max_peaks

    Parameters
    ----------
    data : GCMSData preprocesado (con tic_corrected)
    max_peaks : máximo de picos a retener (None = sin límite)

    Returns
    -------
    PeakResult con índices, propiedades y tabla.
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw

    peaks, props = detect_peaks(tic)
    peaks, props = filter_peaks_snr(tic, peaks, props)
    table = build_peak_table(data, peaks, props)

    # Recortar a max_peaks si se especifica (la tabla ya está ordenada por intensidad)
    if max_peaks is not None and len(table) > max_peaks:
        logger.info("  Recortando de %d a %d picos (MAX_PEAKS=%d)", len(table), max_peaks, max_peaks)
        table = table.head(max_peaks).reset_index(drop=True)
        table["peak_id"] = range(1, len(table) + 1)

    # Reconstruir indices alineados con la tabla final
    final_indices = table["scan_idx"].values.astype(np.intp)

    result = PeakResult(
        indices=final_indices,
        properties=props,
        table=table,
        n_peaks=len(table),
    )
    logger.info("  Resultado final: %d picos detectados y filtrados.", result.n_peaks)
    return result
