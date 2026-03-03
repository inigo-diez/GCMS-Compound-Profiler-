"""
Modulo 3 -- Extraccion de espectros de masas por pico.

Funciones:
    - extract_apex_spectrum: extrae el espectro en el apex de un pico.
    - extract_averaged_spectrum: promedia scans alrededor del apex.
    - detect_coelution: detecta indicios de co-elucion (picos asimetricos/hombros).
    - extract_all_spectra: extrae espectros para todos los picos detectados.
    - export_msp: exporta espectros a formato .msp para busqueda en NIST.
"""
from __future__ import annotations

import builtins
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import skew

import config
from preprocessing import GCMSData, get_spectrum
from peak_detection import PeakResult

logger = logging.getLogger(__name__)


# Exponer `config` en builtins para notebooks que usan `config.*`
# tras `from spectra import ...` sin `import config` explícito.
if not hasattr(builtins, "config"):
    builtins.config = config


@dataclass
class SpectrumData:
    """Espectro de masas extraido de un pico."""

    peak_id: int
    scan_idx: int
    rt_sec: float
    rt_min: float
    mz: NDArray[np.float32]
    intensities: NDArray[np.float32]
    is_coeluted: bool = False
    coelution_note: str = ""
    subpeak: int = 0  # 0=simple, 1=subpico izquierdo, 2=subpico derecho

    @property
    def base_peak_mz(self) -> float:
        """m/z del fragmento mas intenso."""
        if len(self.intensities) == 0:
            return 0.0
        return float(self.mz[np.argmax(self.intensities)])

    @property
    def base_peak_intensity(self) -> float:
        """Intensidad del fragmento mas intenso."""
        if len(self.intensities) == 0:
            return 0.0
        return float(np.max(self.intensities))

    @property
    def n_fragments(self) -> int:
        """Numero de fragmentos en el espectro."""
        return len(self.mz)


def extract_apex_spectrum(
    data: GCMSData,
    scan_idx: int,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """
    Extrae el espectro de masas exacto en el scan del apex.

    Parameters
    ----------
    data : GCMSData con datos espectrales
    scan_idx : indice del scan (apex del pico)

    Returns
    -------
    (mz, intensities)
    """
    return get_spectrum(data, scan_idx)


def extract_averaged_spectrum(
    data: GCMSData,
    scan_idx: int,
    n_scans: int = config.APEX_SCANS_AVERAGE,
) -> tuple[NDArray[np.float32], NDArray[np.float32]]:
    """
    Extrae y promedia espectros de varios scans alrededor del apex.

    Esto mejora la calidad del espectro reduciendo ruido espectral.

    Parameters
    ----------
    data : GCMSData con datos espectrales
    scan_idx : indice del scan central (apex)
    n_scans : numero de scans a promediar (e.g. 3 = apex +/- 1)

    Returns
    -------
    (mz_common, intensities_averaged)
    """
    half = n_scans // 2
    start = max(0, scan_idx - half)
    end = min(data.n_scans, scan_idx + half + 1)

    # Recoger todos los espectros del rango
    all_mz: list[NDArray[np.float32]] = []
    all_int: list[NDArray[np.float32]] = []
    for idx in range(start, end):
        mz, ints = get_spectrum(data, idx)
        all_mz.append(mz)
        all_int.append(ints)

    if len(all_mz) == 1:
        return all_mz[0], all_int[0]

    # Unir todos los m/z unicos y construir una matriz comun
    mz_union = np.unique(np.concatenate(all_mz))

    # Interpolar intensidades de cada scan al grid comun
    intensity_matrix = np.zeros((len(all_mz), len(mz_union)), dtype=np.float32)
    for i, (mz, ints) in enumerate(zip(all_mz, all_int)):
        # Mapear cada m/z al indice mas cercano en mz_union
        indices = np.searchsorted(mz_union, mz)
        indices = np.clip(indices, 0, len(mz_union) - 1)
        np.add.at(intensity_matrix[i], indices, ints)

    # Promediar
    avg_intensities = np.mean(intensity_matrix, axis=0).astype(np.float32)

    # Eliminar fragmentos con intensidad ~0
    mask = avg_intensities > 0
    return mz_union[mask], avg_intensities[mask]


def detect_coelution(
    data: GCMSData,
    scan_idx: int,
    peak_width_scans: float,
) -> tuple[bool, str]:
    """
    Detecta indicios de co-elucion analizando la forma del pico.

    Criterios:
    1. Asimetria (skewness) del perfil del pico
    2. Presencia de hombros (sub-picos dentro del rango del pico)

    Parameters
    ----------
    data : GCMSData preprocesado
    scan_idx : indice del apex del pico
    peak_width_scans : ancho estimado del pico en scans

    Returns
    -------
    (is_coeluted, note)
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw

    # Definir la region del pico (2x ancho a cada lado del apex)
    half_window = max(int(peak_width_scans * 2), 5)
    start = max(0, scan_idx - half_window)
    end = min(len(tic), scan_idx + half_window + 1)
    segment = tic[start:end]

    if len(segment) < 5:
        return False, ""

    notes = []
    is_coeluted = False

    # 1. Asimetria
    seg_skewness = float(skew(segment))
    if abs(seg_skewness) > 1.5:
        is_coeluted = True
        direction = "cola derecha" if seg_skewness > 0 else "cola izquierda"
        notes.append(f"asimetria={seg_skewness:.2f} ({direction})")

    # 2. Hombros: buscar sub-picos en el segmento
    sub_peaks, _ = find_peaks(segment, prominence=np.max(segment) * 0.1)
    if len(sub_peaks) > 1:
        is_coeluted = True
        notes.append(f"{len(sub_peaks)} sub-picos detectados")

    note = "; ".join(notes) if notes else ""
    return is_coeluted, note


def _gaussian(x, a, mu, sigma):
    """Modelo de 1 Gaussiana."""
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _double_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    """Modelo de suma de 2 Gaussianas."""
    return (_gaussian(x, a1, mu1, sigma1) + _gaussian(x, a2, mu2, sigma2))


def _compute_bic(n, rss, k):
    """Bayesian Information Criterion: BIC = n*ln(RSS/n) + k*ln(n)."""
    if rss <= 0 or n <= 0:
        return np.inf
    return n * np.log(rss / n) + k * np.log(n)


def detect_coelution_gaussian(
    data: GCMSData,
    scan_idx: int,
    peak_width_scans: float,
    bic_delta: float = config.COELUTION_BIC_DELTA,
) -> tuple[bool, int, int, str]:
    """
    Detecta co-elucion ajustando 1 vs 2 Gaussianas y comparando BIC.

    Parameters
    ----------
    data : GCMSData preprocesado
    scan_idx : indice del apex del pico
    peak_width_scans : ancho estimado del pico en scans
    bic_delta : delta minimo de BIC para preferir 2 Gaussianas

    Returns
    -------
    (is_coeluted, mu1_scan_idx, mu2_scan_idx, note)
    Si no co-eluido, mu1_scan_idx = mu2_scan_idx = scan_idx.
    """
    tic = data.tic_corrected if data.tic_corrected is not None else data.tic_raw

    # Ventana alrededor del pico
    half_window = max(int(peak_width_scans * 2.5), 8)
    start = max(0, scan_idx - half_window)
    end = min(len(tic), scan_idx + half_window + 1)
    segment = tic[start:end].astype(np.float64)
    n = len(segment)

    if n < 8:
        return False, scan_idx, scan_idx, ""

    x = np.arange(n, dtype=np.float64)
    apex_local = scan_idx - start  # posicion del apex en coordenadas locales
    y_max = float(segment.max())
    sigma_est = max(peak_width_scans / 2.35, 1.0)

    # --- Ajuste 1 Gaussiana ---
    try:
        p0_1g = [y_max, float(apex_local), sigma_est]
        bounds_1g = ([0, 0, 0.5], [y_max * 3, n - 1, n / 2])
        popt_1g, _ = curve_fit(_gaussian, x, segment, p0=p0_1g, bounds=bounds_1g, maxfev=2000)
        pred_1g = _gaussian(x, *popt_1g)
        rss_1g = float(np.sum((segment - pred_1g) ** 2))
        bic_1g = _compute_bic(n, rss_1g, k=3)
    except (RuntimeError, ValueError):
        return False, scan_idx, scan_idx, "ajuste 1G fallo"

    # --- Ajuste 2 Gaussianas ---
    try:
        offset = max(peak_width_scans / 3, 2.0)
        p0_2g = [
            y_max * 0.6, max(apex_local - offset, 1.0), sigma_est * 0.8,
            y_max * 0.4, min(apex_local + offset, n - 2.0), sigma_est * 0.8,
        ]
        bounds_lo = [0, 0, 0.5, 0, 0, 0.5]
        bounds_hi = [y_max * 3, n - 1, n / 2, y_max * 3, n - 1, n / 2]
        popt_2g, _ = curve_fit(
            _double_gaussian, x, segment, p0=p0_2g,
            bounds=(bounds_lo, bounds_hi), maxfev=5000,
        )
        pred_2g = _double_gaussian(x, *popt_2g)
        rss_2g = float(np.sum((segment - pred_2g) ** 2))
        bic_2g = _compute_bic(n, rss_2g, k=6)
    except (RuntimeError, ValueError):
        return False, scan_idx, scan_idx, "ajuste 2G fallo"

    # --- Comparar BIC ---
    is_coeluted = bic_2g < bic_1g - bic_delta

    if is_coeluted:
        # Obtener posiciones de los subpicos (en coordenadas absolutas de scan)
        mu1_local = popt_2g[1]
        mu2_local = popt_2g[4]
        # Asegurar que mu1 < mu2
        if mu1_local > mu2_local:
            mu1_local, mu2_local = mu2_local, mu1_local
        mu1_abs = int(np.clip(round(start + mu1_local), 0, len(tic) - 1))
        mu2_abs = int(np.clip(round(start + mu2_local), 0, len(tic) - 1))

        # Verificar separacion minima (al menos 2 scans)
        if abs(mu1_abs - mu2_abs) < 2:
            return False, scan_idx, scan_idx, f"BIC favorece 2G pero separacion insuficiente (delta_BIC={bic_1g - bic_2g:.1f})"

        note = (f"co-elucion Gaussiana: delta_BIC={bic_1g - bic_2g:.1f}, "
                f"mu1=scan {mu1_abs} ({data.scan_times_min[mu1_abs]:.2f} min), "
                f"mu2=scan {mu2_abs} ({data.scan_times_min[mu2_abs]:.2f} min)")
        return True, mu1_abs, mu2_abs, note
    else:
        return False, scan_idx, scan_idx, ""


def extract_all_spectra(
    data: GCMSData,
    peak_result: PeakResult,
    average_scans: int = config.APEX_SCANS_AVERAGE,
    check_coelution: bool = True,
    coelution_avg_scans: int = config.SPECTRA_AVG_SCANS,
) -> list[SpectrumData]:
    """
    Extrae espectros de masas para todos los picos detectados.

    Si co-elucion Gaussiana detecta 2 subpicos, genera 2 SpectrumData
    con el mismo peak_id pero subpeak=1 y subpeak=2.

    Parameters
    ----------
    data : GCMSData preprocesado
    peak_result : resultado de la deteccion de picos
    average_scans : scans a promediar para picos simples
    check_coelution : si True, evalua co-elucion Gaussiana
    coelution_avg_scans : scans a promediar en cada subpico co-eluido

    Returns
    -------
    Lista de SpectrumData. Puede ser > n_peaks si hay subpicos.
    """
    spectra: list[SpectrumData] = []
    df = peak_result.table

    logger.info("  Extrayendo espectros de %d picos (promediando %d scans)...",
                len(df), average_scans)

    n_coeluted = 0
    n_subpeaks = 0
    for _, row in df.iterrows():
        peak_id = int(row["peak_id"])
        scan_idx = int(row["scan_idx"])
        rt_sec = float(row["rt_sec"])
        rt_min = float(row["rt_min"])
        width_scans = float(row["width_scans"])

        # Detectar co-elucion Gaussiana
        is_coeluted = False
        coelution_note = ""
        mu1_idx = scan_idx
        mu2_idx = scan_idx

        if check_coelution:
            is_coeluted, mu1_idx, mu2_idx, coelution_note = detect_coelution_gaussian(
                data, scan_idx, width_scans,
            )
            if is_coeluted:
                n_coeluted += 1

        if is_coeluted:
            # Generar 2 subpicos
            for sub, sub_scan_idx in [(1, mu1_idx), (2, mu2_idx)]:
                if coelution_avg_scans > 1:
                    mz, intensities = extract_averaged_spectrum(data, sub_scan_idx, coelution_avg_scans)
                else:
                    mz, intensities = extract_apex_spectrum(data, sub_scan_idx)

                spectrum = SpectrumData(
                    peak_id=peak_id,
                    scan_idx=sub_scan_idx,
                    rt_sec=data.scan_times[sub_scan_idx],
                    rt_min=data.scan_times_min[sub_scan_idx],
                    mz=mz,
                    intensities=intensities,
                    is_coeluted=True,
                    coelution_note=coelution_note,
                    subpeak=sub,
                )
                spectra.append(spectrum)
            n_subpeaks += 2
        else:
            # Pico simple
            if average_scans > 1:
                mz, intensities = extract_averaged_spectrum(data, scan_idx, average_scans)
            else:
                mz, intensities = extract_apex_spectrum(data, scan_idx)

            spectrum = SpectrumData(
                peak_id=peak_id,
                scan_idx=scan_idx,
                rt_sec=rt_sec,
                rt_min=rt_min,
                mz=mz,
                intensities=intensities,
                is_coeluted=False,
                coelution_note="",
                subpeak=0,
            )
            spectra.append(spectrum)

    logger.info("  %d espectros extraidos (%d picos co-eluidos -> %d subpicos).",
                len(spectra), n_coeluted, n_subpeaks)
    return spectra


def _normalize_spectrum(
    mz: NDArray[np.float32],
    intensities: NDArray[np.float32],
    max_intensity: float = 999.0,
) -> tuple[NDArray[np.float32], NDArray[np.int32]]:
    """Normaliza intensidades a escala 0-999 (formato NIST)."""
    if len(intensities) == 0 or np.max(intensities) == 0:
        return mz, np.zeros(len(mz), dtype=np.int32)
    normalized = (intensities / np.max(intensities)) * max_intensity
    return mz, np.round(normalized).astype(np.int32)


def export_msp(
    spectra: list[SpectrumData],
    output_path: str | Path,
    sample_name: str = "Unknown",
) -> Path:
    """
    Exporta espectros a formato .msp (NIST/AMDIS compatible).

    Cada espectro se escribe como un bloque con:
    - NAME, RT, COMMENT, NUM PEAKS
    - Listado de pares (m/z, intensidad normalizada)

    Parameters
    ----------
    spectra : lista de SpectrumData
    output_path : ruta del archivo .msp de salida
    sample_name : nombre de la muestra

    Returns
    -------
    Path al archivo exportado.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        for sp in spectra:
            mz_norm, int_norm = _normalize_spectrum(sp.mz, sp.intensities)

            # Filtrar fragmentos con intensidad 0 tras normalizar
            mask = int_norm > 0
            mz_out = mz_norm[mask]
            int_out = int_norm[mask]

            sub_suffix = f"_{chr(96 + sp.subpeak)}" if sp.subpeak > 0 else ""
            f.write(f"NAME: {sample_name}_peak{sp.peak_id:04d}{sub_suffix}\n")
            f.write(f"RETENTIONTIME: {sp.rt_sec:.2f}\n")
            f.write(f"RETENTIONINDEX: {sp.rt_min:.4f}\n")
            comment_parts = [f"RT={sp.rt_min:.2f}min"]
            if sp.is_coeluted:
                comment_parts.append(f"COELUTION: {sp.coelution_note}")
            f.write(f"COMMENT: {'; '.join(comment_parts)}\n")
            f.write(f"NUM PEAKS: {len(mz_out)}\n")
            for m, i in zip(mz_out, int_out):
                f.write(f"{m:.1f} {i}\n")
            f.write("\n")

    logger.info("  Exportados %d espectros a %s", len(spectra), output_path)
    return output_path
