"""
Modulo 5 -- Alineamiento de picos entre muestras y construccion de la matriz de datos.

Funciones:
    - align_peaks: alinea picos entre multiples muestras por RT.
    - build_data_matrix: construye la matriz muestras x metabolitos.
    - export_matrix: exporta a CSV y Excel.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import NDArray

import config
from peak_detection import PeakResult
from spectra import SpectrumData
from nist_search import MatchResult

logger = logging.getLogger(__name__)


@dataclass
class AlignedFeature:
    """Un metabolito/feature alineado entre muestras."""

    feature_id: int
    mean_rt_min: float
    compound_name: str
    cas: str
    formula: str
    best_score: float
    # Mapeo: sample_name -> {rt_min, intensity, area, peak_id}
    sample_data: dict[str, dict] = field(default_factory=dict)

    @property
    def n_detected(self) -> int:
        """Numero de muestras donde se detecto este feature."""
        return len(self.sample_data)


@dataclass
class SampleResult:
    """Resultado completo del procesamiento de una muestra."""

    name: str
    peak_table: pd.DataFrame
    spectra: list[SpectrumData]
    matches: list[MatchResult]


def align_peaks(
    sample_results: dict[str, SampleResult],
    rt_tolerance: float = config.RT_TOLERANCE_MIN,
) -> list[AlignedFeature]:
    """
    Alinea picos entre multiples muestras por tiempo de retencion.

    Algoritmo:
    1. Recoger todos los picos de todas las muestras.
    2. Ordenar por RT.
    3. Agrupar picos dentro de la tolerancia de RT.
    4. Para cada grupo, asignar la mejor identificacion.

    Parameters
    ----------
    sample_results : dict de nombre_muestra -> SampleResult
    rt_tolerance : tolerancia de RT en minutos

    Returns
    -------
    Lista de AlignedFeature.
    """
    logger.info("  Alineando picos de %d muestras (tolerancia=%.3f min)...",
                len(sample_results), rt_tolerance)

    # 1. Recoger todos los picos
    all_peaks: list[dict] = []
    for sample_name, result in sample_results.items():
        # Construir lookup de matches por peak_id (soporta subpicos)
        match_lookup: dict[int, MatchResult] = {}
        for m in result.matches:
            if m.peak_id not in match_lookup or m.best_score > match_lookup[m.peak_id].best_score:
                match_lookup[m.peak_id] = m

        for i, row in result.peak_table.iterrows():
            peak_id = int(row["peak_id"])

            # Buscar la identificacion de este pico
            compound = "Unknown"
            cas = ""
            formula = ""
            score = 0.0
            match = match_lookup.get(peak_id)
            if match and match.best_match:
                compound = match.best_name
                cas = match.best_match.get("cas", "")
                formula = match.best_match.get("formula", "")
                score = match.best_score

            all_peaks.append({
                "sample": sample_name,
                "peak_id": peak_id,
                "rt_min": float(row["rt_min"]),
                "intensity": float(row["intensity"]),
                "compound": compound,
                "cas": cas,
                "formula": formula,
                "score": score,
            })

    if not all_peaks:
        logger.warning("  No hay picos para alinear.")
        return []

    # 2. Ordenar por RT
    all_peaks.sort(key=lambda x: x["rt_min"])

    # 3. Agrupar por RT (greedy clustering)
    features: list[AlignedFeature] = []
    used = set()

    for i, peak in enumerate(all_peaks):
        if i in used:
            continue

        # Iniciar nuevo grupo
        group = [peak]
        group_indices = {i}
        group_rts = [peak["rt_min"]]

        # Buscar otros picos dentro de la tolerancia
        for j in range(i + 1, len(all_peaks)):
            if j in used:
                continue
            other = all_peaks[j]
            # Si el RT se aleja demasiado del centro del grupo, parar
            mean_rt = np.mean(group_rts)
            if other["rt_min"] - mean_rt > rt_tolerance * 2:
                break
            if abs(other["rt_min"] - mean_rt) <= rt_tolerance:
                # No repetir la misma muestra en un grupo
                if other["sample"] not in {p["sample"] for p in group}:
                    group.append(other)
                    group_indices.add(j)
                    group_rts.append(other["rt_min"])

        used.update(group_indices)

        # Determinar la mejor identificacion del grupo
        best_peak = max(group, key=lambda x: x["score"])
        compound_name = best_peak["compound"]
        cas = best_peak["cas"]
        formula_str = best_peak["formula"]
        best_score = best_peak["score"]

        # Si no hay identificacion, usar RT como nombre
        if compound_name == "Unknown" or compound_name == "No match":
            mean_rt = np.mean(group_rts)
            compound_name = f"Unknown_RT{mean_rt:.2f}"

        feature = AlignedFeature(
            feature_id=len(features) + 1,
            mean_rt_min=float(np.mean(group_rts)),
            compound_name=compound_name,
            cas=cas,
            formula=formula_str,
            best_score=best_score,
        )

        # Rellenar datos por muestra
        for p in group:
            feature.sample_data[p["sample"]] = {
                "rt_min": p["rt_min"],
                "intensity": p["intensity"],
                "peak_id": p["peak_id"],
            }

        features.append(feature)

    # Ordenar features por RT
    features.sort(key=lambda f: f.mean_rt_min)
    # Re-numerar
    for i, f in enumerate(features):
        f.feature_id = i + 1

    logger.info("  Alineamiento completado: %d features de %d picos totales.",
                len(features), len(all_peaks))
    return features


def build_data_matrix(
    features: list[AlignedFeature],
    sample_names: list[str],
    value_type: str = "intensity",
) -> pd.DataFrame:
    """
    Construye la matriz de datos (muestras x metabolitos).

    Parameters
    ----------
    features : lista de AlignedFeature
    sample_names : lista ordenada de nombres de muestras
    value_type : "intensity" (default)

    Returns
    -------
    DataFrame con muestras en filas y metabolitos en columnas.
    """
    logger.info("  Construyendo matriz de datos (%d muestras x %d features)...",
                len(sample_names), len(features))

    # Crear nombres de columnas unicos
    col_names = []
    seen_names: dict[str, int] = {}
    for f in features:
        name = f.compound_name
        if name in seen_names:
            seen_names[name] += 1
            name = f"{name}_{seen_names[name]}"
        else:
            seen_names[name] = 1
        col_names.append(name)

    # Construir la matriz
    matrix = np.full((len(sample_names), len(features)), np.nan)

    for j, feature in enumerate(features):
        for i, sample in enumerate(sample_names):
            if sample in feature.sample_data:
                matrix[i, j] = feature.sample_data[sample][value_type]

    df = pd.DataFrame(matrix, index=sample_names, columns=col_names)
    df.index.name = "sample"

    n_missing = int(np.isnan(matrix).sum())
    n_total = matrix.size
    pct_missing = (n_missing / n_total) * 100 if n_total > 0 else 0
    logger.info("  Matriz construida: %d x %d, missing values: %d (%.1f%%)",
                len(sample_names), len(features), n_missing, pct_missing)
    return df


def build_feature_metadata(features: list[AlignedFeature]) -> pd.DataFrame:
    """
    Construye tabla de metadatos de los features (metabolitos).

    Util para documentar que hay en cada columna de la matriz.
    """
    rows = []
    for f in features:
        rows.append({
            "feature_id": f.feature_id,
            "compound_name": f.compound_name,
            "mean_rt_min": round(f.mean_rt_min, 4),
            "cas": f.cas,
            "formula": f.formula,
            "best_match_score": round(f.best_score, 4),
            "n_samples_detected": f.n_detected,
        })
    return pd.DataFrame(rows)


def export_matrix(
    data_matrix: pd.DataFrame,
    feature_metadata: pd.DataFrame,
    output_dir: str | Path = config.OUTPUT_DIR,
    prefix: str = "autogcms",
) -> dict[str, Path]:
    """
    Exporta la matriz de datos y metadatos a CSV y Excel.

    Returns
    -------
    Dict con las rutas de los archivos exportados.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}

    # CSV compatible con Excel ES (sep=';', encoding utf-8-sig)
    csv_matrix = output_dir / f"{prefix}_data_matrix.csv"
    data_matrix.to_csv(csv_matrix, sep=";", encoding="utf-8-sig")
    paths["matrix_csv"] = csv_matrix
    logger.info("  Matriz exportada: %s", csv_matrix)

    csv_meta = output_dir / f"{prefix}_feature_metadata.csv"
    feature_metadata.to_csv(csv_meta, sep=";", encoding="utf-8-sig", index=False)
    paths["metadata_csv"] = csv_meta
    logger.info("  Metadatos exportados: %s", csv_meta)

    # Excel con formato (freeze_panes + auto-width)
    xlsx_path = output_dir / f"{prefix}_results.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        data_matrix.to_excel(writer, sheet_name="Data Matrix")
        feature_metadata.to_excel(writer, sheet_name="Feature Metadata", index=False)
        _format_xlsx_sheets(writer)
    paths["excel"] = xlsx_path
    logger.info("  Excel exportado: %s", xlsx_path)

    return paths


def _format_xlsx_sheets(writer: pd.ExcelWriter) -> None:
    """Aplica freeze_panes y auto-width a todas las hojas del writer."""
    for sheet_name in writer.sheets:
        ws = writer.sheets[sheet_name]
        # Freeze: B2 si hay index (Data Matrix), A2 si no (metadata)
        if sheet_name == "Data Matrix":
            ws.freeze_panes = "B2"
        else:
            ws.freeze_panes = "A2"
        # Auto-ajuste de ancho de columnas (capado a 40)
        for col_cells in ws.columns:
            max_len = 0
            for cell in col_cells:
                val = str(cell.value) if cell.value is not None else ""
                max_len = max(max_len, len(val))
            col_letter = col_cells[0].column_letter
            ws.column_dimensions[col_letter].width = min(max_len + 2, 40)


def export_xlsx_report(
    report_dir: str | Path,
    tables: dict[str, pd.DataFrame],
    filename: str = "autogcms_report.xlsx",
) -> Path:
    """
    Exporta un XLSX multi-hoja con formato para el reporte final.

    Parameters
    ----------
    report_dir : directorio de salida
    tables : dict nombre_hoja -> DataFrame
    filename : nombre del archivo

    Returns
    -------
    Path al archivo exportado.
    """
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    xlsx_path = report_dir / filename

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for sheet_name, df in tables.items():
            has_index = sheet_name == "Data Matrix"
            df.to_excel(writer, sheet_name=sheet_name[:31], index=has_index)
        _format_xlsx_sheets(writer)

    logger.info("  Reporte XLSX exportado: %s (%d hojas)", xlsx_path, len(tables))
    return xlsx_path
