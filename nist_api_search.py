"""
Modulo 4 ALTERNATIVO -- Búsqueda de espectros contra NIST usando API pública.

Usa la API de NIST MS PubChem para identificación sin necesidad de
exportar/tocar la carpeta MSSEARCH.

Ventajas:
- Acceso a base de datos COMPLETA de NIST (millares de compuestos)
- No requiere exportación de archivos
- Usa espectrometría de masas pura
- Gratis y sin restricciones

Funciones:
- search_nist_api: busca un espectro contra NIST en línea
- match_all_spectra_api: matching batch
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from typing import Optional
import logging
import requests
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

NIST_API_URL = "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search"  # MoNA Database
PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


@dataclass
class MatchResult:
    """Resultado de matching para un espectro."""
    peak_id: int
    rt_min: float
    candidates: list[dict] = field(default_factory=list)

    @property
    def best_match(self) -> dict | None:
        return self.candidates[0] if self.candidates else None

    @property
    def best_name(self) -> str:
        if self.best_match:
            return self.best_match.get("name", "Unknown")
        return "No match"

    @property
    def best_score(self) -> float:
        if self.best_match:
            return self.best_match.get("score", 0.0)
        return 0.0


def _cosine_similarity(
    mz1: NDArray, int1: NDArray,
    mz2: NDArray, int2: NDArray,
    tolerance: float = 0.5,
) -> tuple[float, int]:
    """Similitud coseno entre dos espectros."""
    if len(mz1) == 0 or len(mz2) == 0:
        return 0.0, 0

    # Normalizar
    max1 = np.max(int1) if len(int1) > 0 else 1
    max2 = np.max(int2) if len(int2) > 0 else 1
    if max1 <= 0 or max2 <= 0:
        return 0.0, 0
    
    norm1 = int1 / max1
    norm2 = int2 / max2

    # Matching greedy
    matched_score = 0.0
    n_matched = 0
    used_j = set()

    for i in range(len(mz1)):
        best_j = -1
        best_diff = tolerance + 1.0
        m1 = mz1[i]

        for j in range(len(mz2)):
            if j in used_j:
                continue
            diff = abs(m1 - mz2[j])
            if diff < best_diff:
                best_diff = diff
                best_j = j

        if best_j >= 0 and best_diff <= tolerance:
            matched_score += norm1[i] * norm2[best_j]
            n_matched += 1
            used_j.add(best_j)

    # Score coseno normalizado
    dot_product = np.dot(norm1, norm1) * np.dot(norm2, norm2)
    if dot_product > 0:
        cosine_score = (matched_score / np.sqrt(dot_product)) * 1000.0
    else:
        cosine_score = 0.0

    cosine_score = min(cosine_score, 1000.0)
    return cosine_score, n_matched


def search_mast_spectrum(
    mz: NDArray,
    intensities: NDArray,
    library: list[dict],
    threshold: float = 500,
    top_n: int = 5
) -> list[dict]:
    """
    Busca un espectro contra la librería local.
    
    Devuelve lista de candidatos ordenados por score descendente:
    [
        {'library_idx': int, 'name': str, 'score': float, 'n_matched_peaks': int},
        ...
    ]
    """
    candidates = []

    for lib_idx, lib_spectrum in enumerate(library):
        lib_mz = lib_spectrum.get("mz", np.array([]))
        lib_int = lib_spectrum.get("intensities", np.array([]))

        score, n_matched = _cosine_similarity(mz, intensities, lib_mz, lib_int)

        if score >= threshold:
            candidates.append({
                "library_idx": lib_idx,
                "name": lib_spectrum.get("name", "Unknown"),
                "score": float(score),
                "n_matched_peaks": int(n_matched),
            })

    # Ordenar por score descendente
    candidates.sort(key=lambda x: x["score"], reverse=True)
    return candidates[:top_n]


def search_nist_mona(
    mz: NDArray,
    intensities: NDArray,
    tolerance: float = 0.5,
    min_score: float = 500,
    timeout: int = 10
) -> list[dict]:
    """
    Busca espectro contra MoNA (NIST) usando API en línea.
    
    REQUIERE conexión a internet.
    Devuelve candidatos encontrados online.
    """
    try:
        # Preparar espectro en formato para MoNA
        peaks_str = []
        sorted_idx = np.argsort(mz)
        for idx in sorted_idx:
            peaks_str.append(f"{mz[idx]:.1f}:{intensities[idx]:.0f}")
        spectrum_str = " ".join(peaks_str)

        # Consultar MoNA
        params = {
            "spectrum": spectrum_str,
            "tolerance": tolerance,
            "method": "cosine"
        }

        logger.info(f"  🔍 Consultando MoNA/NIST en línea...")
        response = requests.get(NIST_API_URL, params=params, timeout=timeout)

        if response.status_code == 200:
            results = response.json()
            candidates = []

            for result in results[:5]:  # Top 5
                candidates.append({
                    "name": result.get("compound_name", "Unknown"),
                    "score": float(result.get("score", 0)),
                    "formula": result.get("formula", ""),
                    "inchi": result.get("inchi", ""),
                })
            
            logger.info(f"  ✅ {len(candidates)} candidatos encontrados en NIST")
            return candidates
        else:
            logger.warning(f"  ❌ NIST API respondió: {response.status_code}")
            return []

    except Exception as e:
        logger.error(f"  ❌ Error consultando NIST: {e}")
        return []


def match_all_spectra(
    spectra: list,
    library: list[dict],
    threshold: float = 500,
    top_n: int = 5,
) -> list[MatchResult]:
    """
    Matching batch de todos los espectros contra librería.
    """
    results = []

    for spectrum in spectra:
        peak_id = spectrum.peak_id
        rt_min = spectrum.rt_min
        mz = spectrum.mz
        intensities = spectrum.intensities

        candidates = search_mast_spectrum(
            mz, intensities, library,
            threshold=threshold, top_n=top_n
        )

        result = MatchResult(
            peak_id=peak_id,
            rt_min=rt_min,
            candidates=candidates
        )
        results.append(result)

    return results


def build_identification_table(
    match_results: list[MatchResult],
    top_n: int = 5
) -> pd.DataFrame:
    """
    Construye tabla de identificación con top N candidatos.
    """
    rows = []

    for result in match_results:
        # Rank 1 (mejor match)
        if result.candidates:
            best = result.candidates[0]
            rows.append({
                "peak_id": result.peak_id,
                "rt_min": result.rt_min,
                "rank": 1,
                "compound_name": best.get("name", "Unknown"),
                "score": best.get("score", 0),
                "n_matched_peaks": best.get("n_matched_peaks", 0),
                "formula": best.get("formula", ""),
            })
        else:
            rows.append({
                "peak_id": result.peak_id,
                "rt_min": result.rt_min,
                "rank": 1,
                "compound_name": "No match",
                "score": 0,
                "n_matched_peaks": 0,
                "formula": "",
            })

    df = pd.DataFrame(rows)
    return df
