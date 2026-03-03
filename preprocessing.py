"""
Módulo 1 — Lectura y preprocesado de archivos .cdf (ANDI-MS / netCDF).

Funciones:
    - load_cdf: lee un archivo .cdf y devuelve tiempos, TIC y datos espectrales.
    - smooth_tic: aplica suavizado Savitzky-Golay al TIC.
    - correct_baseline: corrección de baseline con Asymmetric Least Squares (AsLS).
    - preprocess: pipeline completo de preprocesado para un archivo .cdf.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
import netCDF4 as nc
import numpy as np
from numpy.typing import NDArray
from scipy.signal import savgol_filter
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

import config

logger = logging.getLogger(__name__)


@dataclass
class GCMSData:
    """Contenedor de datos GC-MS para una muestra."""

    filepath: Path
    name: str
    scan_times: NDArray[np.float64]        # Tiempo de cada scan (segundos)
    scan_times_min: NDArray[np.float64]    # Tiempo en minutos
    tic_raw: NDArray[np.float64]           # TIC original
    tic_smooth: NDArray[np.float64] | None = None
    tic_corrected: NDArray[np.float64] | None = None
    scan_index: NDArray[np.int32] | None = None
    point_count: NDArray[np.int32] | None = None
    mass_values: NDArray[np.float32] | None = None
    intensity_values: NDArray[np.float32] | None = None
    n_scans: int = 0
    mass_range: tuple[float, float] = (0.0, 0.0)


def load_cdf(filepath: str | Path) -> GCMSData:
    """
    Lee un archivo .cdf (ANDI-MS) y extrae todos los datos relevantes.

    Parameters
    ----------
    filepath : ruta al archivo .cdf

    Returns
    -------
    GCMSData con tiempos, TIC e información espectral.
    """
    filepath = Path(filepath)
    logger.info("Leyendo archivo: %s", filepath.name)

    ds = nc.Dataset(str(filepath), "r")
    try:
        scan_times = np.asarray(ds.variables["scan_acquisition_time"][:], dtype=np.float64)
        tic = np.asarray(ds.variables["total_intensity"][:], dtype=np.float64)
        scan_index = np.asarray(ds.variables["scan_index"][:], dtype=np.int32)
        point_count = np.asarray(ds.variables["point_count"][:], dtype=np.int32)
        mass_values = np.asarray(ds.variables["mass_values"][:], dtype=np.float32)
        intensity_values = np.asarray(ds.variables["intensity_values"][:], dtype=np.float32)

        mass_min = float(mass_values.min())
        mass_max = float(mass_values.max())
    finally:
        ds.close()

    data = GCMSData(
        filepath=filepath,
        name=filepath.stem,
        scan_times=scan_times,
        scan_times_min=scan_times / 60.0,
        tic_raw=tic,
        scan_index=scan_index,
        point_count=point_count,
        mass_values=mass_values,
        intensity_values=intensity_values,
        n_scans=len(scan_times),
        mass_range=(mass_min, mass_max),
    )
    logger.info(
        "  %d scans, tiempo %.1f–%.1f s (%.2f–%.2f min), rango m/z %.1f–%.1f",
        data.n_scans,
        scan_times[0], scan_times[-1],
        data.scan_times_min[0], data.scan_times_min[-1],
        mass_min, mass_max,
    )
    return data


def smooth_tic(
    tic: NDArray[np.float64],
    window: int = config.SAVGOL_WINDOW,
    polyorder: int = config.SAVGOL_POLYORDER,
) -> NDArray[np.float64]:
    """Aplica suavizado Savitzky-Golay al TIC."""
    if window % 2 == 0:
        window += 1
    smoothed = savgol_filter(tic, window_length=window, polyorder=polyorder)
    # No permitir valores negativos tras el suavizado
    smoothed = np.maximum(smoothed, 0.0)
    return smoothed


def _asls_baseline(
    y: NDArray[np.float64],
    lam: float = config.BASELINE_LAMBDA,
    p: float = config.BASELINE_P,
    n_iter: int = 10,
) -> NDArray[np.float64]:
    """
    Asymmetric Least Squares baseline estimation (Eilers & Boelens, 2005).

    Parameters
    ----------
    y : señal de entrada
    lam : parámetro de suavidad (mayor = baseline más suave)
    p : asimetría (0 < p < 1, valores pequeños siguen los valles)
    n_iter : iteraciones
    """
    n = len(y)
    # Matriz de diferencias de segundo orden (n x n)
    D = diags([1, -2, 1], [0, 1, 2], shape=(n - 2, n), dtype=np.float64, format="csc")
    DTD = D.T.dot(D)
    w = np.ones(n)
    for _ in range(n_iter):
        W = diags(w, 0, shape=(n, n), format="csc")
        Z = W + lam * DTD
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y <= z)
    return z


def correct_baseline(
    tic: NDArray[np.float64],
    lam: float = config.BASELINE_LAMBDA,
    p: float = config.BASELINE_P,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Corrige la baseline del TIC usando AsLS.

    Returns
    -------
    (tic_corrected, baseline)
    """
    baseline = _asls_baseline(tic, lam=lam, p=p)
    corrected = tic - baseline
    corrected = np.maximum(corrected, 0.0)
    return corrected, baseline


def preprocess(filepath: str | Path) -> GCMSData:
    """
    Pipeline completo de preprocesado para un archivo .cdf:
    1. Lectura del archivo
    2. Suavizado Savitzky-Golay
    3. Corrección de baseline (AsLS)

    Returns
    -------
    GCMSData con tic_smooth y tic_corrected rellenados.
    """
    data = load_cdf(filepath)

    logger.info("  Aplicando suavizado Savitzky-Golay (window=%d, order=%d)",
                config.SAVGOL_WINDOW, config.SAVGOL_POLYORDER)
    data.tic_smooth = smooth_tic(data.tic_raw)

    logger.info("  Corrigiendo baseline (AsLS, lambda=%.0e, p=%.3f)",
                config.BASELINE_LAMBDA, config.BASELINE_P)
    data.tic_corrected, _ = correct_baseline(data.tic_smooth)

    logger.info("  Preprocesado completado.")
    return data


def get_spectrum(data: GCMSData, scan_idx: int) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """
    Extrae el espectro de masas de un scan específico.

    Parameters
    ----------
    data : GCMSData con los datos espectrales
    scan_idx : índice del scan (0-based)

    Returns
    -------
    (mz_array, intensity_array) para ese scan.
    """
    start = data.scan_index[scan_idx]
    n_points = data.point_count[scan_idx]
    end = start + n_points
    mz = data.mass_values[start:end]
    intensities = data.intensity_values[start:end]
    return mz, intensities
