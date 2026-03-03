"""
Modulo 4 -- Matching de espectros contra libreria NIST exportada (.msp).

Estrategia: similitud coseno (Python puro + NumPy) contra un .msp exportado
desde NIST MS Search. No requiere dependencias externas como matchms o numba.

Optimizado para GC-MS con ionizacion electronica (EI):
    - m/z enteros (unit mass resolution)
    - Intensidades normalizadas antes de comparar
    - Tolerancia de 0.5 Da para matching de fragmentos
    - Top N candidatos por pico con score coseno >= umbral

Funciones principales:
    - load_nist_library: carga .msp exportado de NIST.
    - match_all_spectra: matching batch de todos los espectros.
    - build_identification_table: genera tabla con top N candidatos/pico.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

import config
from spectra import SpectrumData

logger = logging.getLogger(__name__)


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


# =====================================================================
# Carga de libreria .msp
# =====================================================================

def _find_msp_files(search_path: Path) -> list[Path]:
    """Busca archivos .msp/.MSP en un directorio o devuelve el archivo directo."""
    if search_path.is_file() and search_path.suffix.lower() == ".msp":
        return [search_path]
    if search_path.is_dir():
        files = list(search_path.glob("*.msp")) + list(search_path.glob("*.MSP"))
        # Deduplicar por nombre en minúsculas
        seen = set()
        unique = []
        for f in files:
            key = str(f).lower()
            if key not in seen:
                seen.add(key)
                unique.append(f)
        return sorted(unique)
    return []


def load_nist_library(msp_path: str | Path | None = None) -> list[dict]:
    """
    Carga una libreria NIST en formato .msp.

    Acepta:
    - Un archivo .msp concreto
    - Un directorio con archivos .MSP (se cargan TODOS)
    - None: busca en NIST_MSP_PATH y luego en NIST_DIR

    Cada espectro se devuelve como dict con claves:
        name, cas, formula, mw, mz (ndarray), intensities (ndarray)
    """
    if msp_path is None:
        msp_path = getattr(config, "NIST_MSP_PATH", None)

    if msp_path is not None:
        msp_path = Path(msp_path)
        msp_files = _find_msp_files(msp_path)
    else:
        msp_files = []

    # Si no encontramos nada, buscar en NIST_DIR
    if not msp_files:
        nist_dir = config.NIST_DIR
        msp_files = _find_msp_files(nist_dir)
        if not msp_files:
            mainlib = nist_dir / "MAINLIB"
            if mainlib.is_dir():
                msp_files = _find_msp_files(mainlib)

    if not msp_files:
        logger.error(
            "No se encontro libreria .msp. Exporta MAINLIB desde NIST MS Search "
            "y configura NIST_MSP_PATH en config.py"
        )
        return []

    # Cargar todos los archivos .msp encontrados
    all_spectra: list[dict] = []
    for msp_file in msp_files:
        loaded = _load_single_msp(msp_file)
        all_spectra.extend(loaded)

    logger.info("  Libreria NIST: %d espectros de referencia cargados de %d archivo(s).",
                len(all_spectra), len(msp_files))
    return all_spectra


def _load_single_msp(msp_path: Path) -> list[dict]:
    """Carga un archivo .msp individual."""
    if not msp_path.exists():
        logger.warning("  Archivo no encontrado: %s", msp_path)
        return []

    file_size_mb = msp_path.stat().st_size / (1024 * 1024)
    logger.info("  Cargando: %s (%.1f MB)", msp_path.name, file_size_mb)

    spectra: list[dict] = []
    current: dict | None = None

    with open(msp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                if current and current.get("_mz"):
                    _finalize_spectrum(current)
                    spectra.append(current)
                current = None
                continue

            if current is None:
                current = {"_mz": [], "_int": []}

            if ":" in line and not line[0].isdigit():
                key, _, val = line.partition(":")
                key = key.strip().upper()
                val = val.strip()
                if key == "NAME":
                    current["name"] = val
                elif key in ("CASNO", "CAS#", "CAS", "CAS NO"):
                    current["cas"] = val
                elif key == "MW":
                    current["mw"] = val
                elif key == "FORMULA":
                    current["formula"] = val
                elif key == "NIST#":
                    current["nist_id"] = val
                elif key == "NUM PEAKS":
                    current["num_peaks"] = int(val)
            elif line[0].isdigit():
                # Parsear pares m/z intensidad (formato: "mz int; mz int;" o "mz int")
                line = line.replace(";", " ")
                parts = line.replace("\t", " ").split()
                i = 0
                while i + 1 < len(parts):
                    try:
                        mz_val = float(parts[i])
                        int_val = float(parts[i + 1])
                        current["_mz"].append(mz_val)
                        current["_int"].append(int_val)
                        i += 2
                    except ValueError:
                        i += 1

        # Ultimo espectro si no termina con linea vacia
        if current and current.get("_mz"):
            _finalize_spectrum(current)
            spectra.append(current)

    logger.info("  Cargados %d espectros de referencia.", len(spectra))
    return spectra


def _finalize_spectrum(sp: dict) -> None:
    """Convierte listas a arrays y normaliza intensidades (unit mass, 0-999)."""
    mz = np.array(sp.pop("_mz"), dtype=np.float32)
    intensities = np.array(sp.pop("_int"), dtype=np.float32)

    # Redondear m/z a enteros (GC-MS EI = unit mass resolution)
    mz = np.round(mz).astype(np.float32)

    # Combinar fragmentos con mismo m/z entero (sumar intensidades)
    unique_mz, inverse = np.unique(mz, return_inverse=True)
    combined_int = np.zeros(len(unique_mz), dtype=np.float32)
    np.add.at(combined_int, inverse, intensities)

    # Normalizar intensidades a 0-999
    max_int = combined_int.max()
    if max_int > 0:
        combined_int = (combined_int / max_int) * 999.0

    sp["mz"] = unique_mz
    sp["intensities"] = combined_int


# =====================================================================
# Motor de similitud coseno (Python puro + NumPy)
# =====================================================================

def _cosine_similarity(
    mz1: NDArray, int1: NDArray,
    mz2: NDArray, int2: NDArray,
    tolerance: float = 0.5,
) -> tuple[float, int]:
    """
    Similitud coseno entre dos espectros de masas.

    Optimizado para GC-MS EI con tolerancia amplia (0.5 Da).
    Devuelve (score, n_matched_peaks).
    """
    if len(mz1) == 0 or len(mz2) == 0:
        return 0.0, 0

    # Normalizar intensidades
    max1 = np.max(int1)
    max2 = np.max(int2)
    if max1 <= 0 or max2 <= 0:
        return 0.0, 0
    norm1 = int1 / max1
    norm2 = int2 / max2

    # Matching greedy de fragmentos dentro de tolerancia
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
            if diff <= tolerance and diff < best_diff:
                best_diff = diff
                best_j = j
        if best_j >= 0:
            used_j.add(best_j)
            matched_score += norm1[i] * norm2[best_j]
            n_matched += 1

    # Normalizar por las normas L2
    norm_a = np.sqrt(np.sum(norm1 ** 2))
    norm_b = np.sqrt(np.sum(norm2 ** 2))
    if norm_a > 0 and norm_b > 0:
        score = float(matched_score / (norm_a * norm_b))
    else:
        score = 0.0

    return score, n_matched


# =====================================================================
# Matching de un espectro contra la libreria
# =====================================================================

def match_spectrum(
    query: SpectrumData,
    library: list[dict],
    top_n: int = config.NIST_TOP_N,
    threshold: float = 0.7,
    min_matched_peaks: int = 3,
) -> MatchResult:
    """
    Busca un espectro query contra toda la libreria.

    Parameters
    ----------
    query : SpectrumData
    library : lista de dicts (de load_nist_library)
    top_n : numero de candidatos a devolver
    threshold : score coseno minimo (0.0 - 1.0)
    min_matched_peaks : minimo de fragmentos coincidentes

    Returns
    -------
    MatchResult con hasta top_n candidatos.
    """
    if not library:
        return MatchResult(peak_id=query.peak_id, rt_min=query.rt_min)

    # Preparar query: redondear m/z y normalizar
    q_mz = np.round(query.mz).astype(np.float32)
    q_int = query.intensities.astype(np.float32)

    # Filtrar fragmentos con intensidad > 0
    mask = q_int > 0
    q_mz = q_mz[mask]
    q_int = q_int[mask]

    if len(q_mz) < 3:
        return MatchResult(peak_id=query.peak_id, rt_min=query.rt_min)

    scores: list[tuple[float, int, int]] = []  # (score, n_matched, lib_idx)

    for i, ref in enumerate(library):
        ref_mz = ref.get("mz", np.array([]))
        ref_int = ref.get("intensities", np.array([]))
        if len(ref_mz) < 3:
            continue

        score, n_matched = _cosine_similarity(q_mz, q_int, ref_mz, ref_int)

        if score >= threshold and n_matched >= min_matched_peaks:
            scores.append((score, n_matched, i))

    scores.sort(key=lambda x: x[0], reverse=True)

    candidates = []
    for score, n_matched, idx in scores[:top_n]:
        ref = library[idx]
        candidates.append({
            "name": ref.get("name", "Unknown"),
            "cas": ref.get("cas", ""),
            "formula": ref.get("formula", ""),
            "mw": ref.get("mw", ""),
            "score": score,
            "score_nist": int(score * 999),
            "rmatch_nist": 0,
            "probability": round(score * 100, 1),
            "n_matched_peaks": n_matched,
            "library_idx": idx,
        })

    return MatchResult(peak_id=query.peak_id, rt_min=query.rt_min, candidates=candidates)


# =====================================================================
# Matching batch
# =====================================================================

def match_all_spectra(
    spectra: list[SpectrumData],
    library: list[dict],
    top_n: int = config.NIST_TOP_N,
    threshold: float = 0.7,
    sample_name: str = "query",
) -> list[MatchResult]:
    """
    Matching batch de todos los espectros contra la libreria .msp.

    Parameters
    ----------
    spectra : lista de SpectrumData extraidos de los picos
    library : lista de dicts (de load_nist_library)
    top_n : candidatos por espectro
    threshold : score coseno minimo (0.0 - 1.0)
    sample_name : nombre de la muestra (para logs)

    Returns
    -------
    Lista de MatchResult, uno por espectro.
    """
    if not library:
        logger.warning("  Sin libreria .msp. No se pueden identificar picos.")
        logger.warning("  Exporta MAINLIB desde NIST MS Search y configura NIST_MSP_PATH en config.py")
        return [MatchResult(peak_id=sp.peak_id, rt_min=sp.rt_min) for sp in spectra]

    logger.info("  Metodo: similitud coseno (tolerance=0.5 Da)")
    logger.info("  Libreria: %d espectros de referencia", len(library))
    logger.info("  Query: %d espectros de %s", len(spectra), sample_name)
    logger.info("  Umbral: score >= %.2f, min 3 fragmentos coincidentes", threshold)

    results: list[MatchResult] = []
    n_identified = 0

    for i, sp in enumerate(spectra):
        result = match_spectrum(sp, library, top_n=top_n, threshold=threshold)
        results.append(result)
        if result.best_match:
            n_identified += 1
        if (i + 1) % 10 == 0 or (i + 1) == len(spectra):
            logger.info("    Procesados %d/%d espectros (%d identificados)...",
                        i + 1, len(spectra), n_identified)

    logger.info("  Matching completado: %d/%d picos identificados (score >= %.2f).",
                n_identified, len(spectra), threshold)
    return results


# =====================================================================
# Tabla de identificaciones
# =====================================================================

def build_identification_table(
    spectra: list[SpectrumData],
    matches: list[MatchResult],
    top_n: int = config.NIST_TOP_N,
) -> pd.DataFrame:
    """
    Construye tabla de identificaciones con los top N candidatos por pico.

    Cada pico genera hasta top_n filas (una por candidato), con un campo
    'rank' que indica la posicion (1=mejor match, 2=segundo, etc.).
    """
    rows = []
    for sp, match in zip(spectra, matches):
        base_row = {
            "peak_id": sp.peak_id,
            "subpeak": sp.subpeak,
            "rt_min": round(sp.rt_min, 4),
            "rt_sec": round(sp.rt_sec, 2),
            "base_peak_mz": sp.base_peak_mz,
            "n_fragments": sp.n_fragments,
        }
        if not match.candidates:
            rows.append({
                **base_row,
                "rank": 1,
                "compound_name": "No match",
                "cas": "",
                "formula": "",
                "match_score": 0.0,
                "match_score_nist": 0,
                "n_matched_peaks": 0,
                "n_candidates": 0,
            })
        else:
            for rank, cand in enumerate(match.candidates[:top_n], start=1):
                rows.append({
                    **base_row,
                    "rank": rank,
                    "compound_name": cand.get("name", "Unknown"),
                    "cas": cand.get("cas", ""),
                    "formula": cand.get("formula", ""),
                    "match_score": round(cand.get("score", 0.0), 4),
                    "match_score_nist": cand.get("score_nist", 0),
                    "n_matched_peaks": cand.get("n_matched_peaks", 0),
                    "n_candidates": len(match.candidates),
                })

    df = pd.DataFrame(rows)
    n_identified = df.loc[df["rank"] == 1, "compound_name"].ne("No match").sum()
    logger.info("  Tabla: %d picos, %d identificados (top %d por pico).",
                len(spectra), n_identified, top_n)
    return df
