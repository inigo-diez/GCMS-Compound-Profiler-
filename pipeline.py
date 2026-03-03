"""
Pipeline orquestador -- Ejecuta los modulos del pipeline AutoGCMS.

Funciones:
    - run_single_sample: modulos 1-4 para un archivo .cdf
    - run_full_pipeline: pipeline completo end-to-end (modulos 1-6)
    - list_cdf_files: lista archivos .cdf en un directorio
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config
from preprocessing import preprocess, GCMSData
from peak_detection import find_and_filter_peaks, PeakResult
from spectra import extract_all_spectra, export_msp, SpectrumData
from nist_search import (
    load_nist_library, match_all_spectra, build_identification_table, MatchResult,
)
from matrix_builder import (
    SampleResult, align_peaks, build_data_matrix, build_feature_metadata,
    export_matrix, export_xlsx_report,
)
from visualization import (
    plot_tic_preprocessing, plot_peaks, plot_peak_quality,
    plot_tic_overlay, plot_match_score_histogram,
    plot_data_matrix_heatmap, plot_sample_summary, plot_spectrum,
)
from qc import (
    build_tic_matrix, pca_tic_qc, compute_qc_summary, plot_pca_tic,
)


def setup_logging(level: str = config.LOGGING_LEVEL) -> None:
    """Configura logging para el pipeline."""
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s | %(name)-20s | %(levelname)-7s | %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )


def run_single_sample(
    filepath: str | Path,
    library: list | None = None,
    save_plots: bool = True,
    save_msp: bool = True,
) -> tuple[GCMSData, PeakResult, list[SpectrumData], list[MatchResult], pd.DataFrame]:
    """
    Ejecuta los modulos 1-4 para un solo archivo .cdf.

    Parameters
    ----------
    filepath : ruta al archivo .cdf
    library : libreria de referencia (si None, no hace matching)
    save_plots : si True, guarda los plots en OUTPUT_DIR
    save_msp : si True, exporta espectros a .msp

    Returns
    -------
    (GCMSData, PeakResult, spectra, matches, identification_table)
    """
    filepath = Path(filepath)
    logger = logging.getLogger("pipeline")
    logger.info("=" * 60)
    logger.info("Procesando: %s", filepath.name)
    logger.info("=" * 60)

    output_dir = config.OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Modulo 1 -- Preprocesado
    logger.info("--- Modulo 1: Preprocesado ---")
    data = preprocess(filepath)

    # Modulo 2 -- Deteccion de picos
    logger.info("--- Modulo 2: Deteccion de picos ---")
    peaks = find_and_filter_peaks(data)

    # Modulo 3 -- Extraccion de espectros (con co-elucion Gaussiana)
    logger.info("--- Modulo 3: Extraccion de espectros ---")
    spectra = extract_all_spectra(data, peaks)

    if save_msp:
        msp_path = output_dir / f"{data.name}_spectra.msp"
        export_msp(spectra, msp_path, sample_name=data.name)

    # Modulo 4 -- Matching contra NIST
    logger.info("--- Modulo 4: Matching ---")
    if spectra:
        matches = match_all_spectra(spectra, library, sample_name=data.name)
        id_table = build_identification_table(spectra, matches)
    else:
        logger.warning(
            "No se extrajeron espectros para %s; se omite el matching.", data.name,
        )
        matches = []
        id_table = pd.DataFrame(columns=[
            "peak_id",
            "subpeak",
            "rt_min",
            "rank",
            "compound_name",
            "formula",
            "cas",
            "match_score",
            "n_matched_peaks",
        ])

    # Guardar outputs
    if save_plots:
        plot_tic_preprocessing(data, save_path=output_dir / f"{data.name}_preprocessing.png")
        plot_peaks(data, peaks, save_path=output_dir / f"{data.name}_peaks.png")
        plot_peak_quality(peaks, save_path=output_dir / f"{data.name}_peak_quality.png")
        logger.info("Plots guardados en: %s", output_dir)

    # Guardar tablas (CSV compatible con Excel ES)
    peaks.table.to_csv(output_dir / f"{data.name}_peaks.csv",
                       sep=";", encoding="utf-8-sig", index=False)
    id_table.to_csv(output_dir / f"{data.name}_identifications.csv",
                    sep=";", encoding="utf-8-sig", index=False)
    logger.info("Tablas guardadas en: %s", output_dir)

    n_id = id_table.loc[id_table["rank"] == 1, "compound_name"].ne("No match").sum()
    logger.info("Resultado: %d picos, %d espectros, %d identificados.",
                peaks.n_peaks, len(spectra), n_id)
    return data, peaks, spectra, matches, id_table


def run_full_pipeline(
    data_dir: str | Path | None = None,
    nist_msp_path: str | Path | None = None,
    save_plots: bool = True,
    library: list | None = None,
) -> dict:
    """
    Pipeline completo end-to-end (modulos 1-6).

    Parameters
    ----------
    data_dir : directorio con archivos .cdf (default: config.DATA_DIR)
    nist_msp_path : ruta al archivo .msp de la libreria NIST (o None)
    save_plots : guardar graficos
    library : libreria NIST precargada (si None, se carga automaticamente)

    Returns
    -------
    Dict con todos los resultados:
        - sample_results: dict[name -> SampleResult]
        - features: lista de AlignedFeature
        - data_matrix: DataFrame
        - feature_metadata: DataFrame
        - sample_stats: DataFrame
        - qc_summary: DataFrame
        - export_paths: dict de rutas exportadas
    """
    logger = logging.getLogger("pipeline")
    logger.info("=" * 60)
    logger.info("AUTOGCMS -- Pipeline completo")
    logger.info("=" * 60)

    output_dir = config.OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Listar archivos
    cdf_files = list_cdf_files(data_dir)
    logger.info("Encontrados %d archivos .cdf en %s", len(cdf_files),
                data_dir or config.DATA_DIR)

    if not cdf_files:
        logger.error("No se encontraron archivos .cdf!")
        return {}

    # Cargar libreria NIST (si no se paso una ya cargada)
    if library is None:
        logger.info("--- Cargando libreria de referencia ---")
        library = load_nist_library(nist_msp_path)
    else:
        logger.info("--- Libreria NIST precargada: %d espectros ---", len(library))
    if not library:
        logger.warning("No se cargo libreria .msp. Los picos no seran identificados.")
        logger.info("Exporta MAINLIB desde NIST MS Search y configura NIST_MSP_PATH en config.py")

    # Procesar cada muestra (modulos 1-4)
    sample_results: dict[str, SampleResult] = {}
    all_gcms_data: dict[str, GCMSData] = {}
    all_id_tables: list[pd.DataFrame] = []

    for cdf_path in cdf_files:
        data, peaks, spectra, matches, id_table = run_single_sample(
            cdf_path, library=library, save_plots=save_plots,
        )
        sample_results[data.name] = SampleResult(
            name=data.name,
            peak_table=peaks.table,
            spectra=spectra,
            matches=matches,
        )
        all_gcms_data[data.name] = data
        id_table_with_sample = id_table.copy()
        id_table_with_sample.insert(0, "sample", data.name)
        all_id_tables.append(id_table_with_sample)

    # Modulo 5 -- Alineamiento y matriz
    logger.info("=" * 60)
    logger.info("--- Modulo 5: Alineamiento y matriz ---")
    features = align_peaks(sample_results)

    sample_names = [sr.name for sr in sample_results.values()]
    data_matrix = build_data_matrix(features, sample_names)
    feature_metadata = build_feature_metadata(features)

    export_paths = export_matrix(data_matrix, feature_metadata, output_dir)

    # Tabla combinada de todas las identificaciones
    all_ids = pd.DataFrame()
    if all_id_tables:
        all_ids = pd.concat(all_id_tables, ignore_index=True)
        all_ids_path = output_dir / "autogcms_all_identifications.csv"
        all_ids.to_csv(all_ids_path, sep=";", encoding="utf-8-sig", index=False)
        export_paths["all_identifications"] = all_ids_path
        logger.info("  Tabla combinada de identificaciones: %s", all_ids_path)

    # ─────────────────────────────────────────────────────────────
    # Modulo 6 -- Reportes, QC y visualizacion -> final_report/
    # ─────────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("--- Modulo 6: Reportes, QC y visualizacion ---")

    report_dir = output_dir / "final_report"
    report_dir.mkdir(parents=True, exist_ok=True)

    # --- A) QC summary ---
    logger.info("  Generando QC summary...")
    qc_summary = compute_qc_summary(sample_results, all_gcms_data)

    # Añadir n_identified y pct_missing al qc_summary
    n_features_total = len(features)
    for i, (name, sr) in enumerate(sample_results.items()):
        n_identified = sum(1 for m in sr.matches if m.best_match is not None)
        qc_summary.loc[i, "n_identified"] = n_identified
        if n_features_total > 0 and name in data_matrix.index:
            n_present = data_matrix.loc[name].notna().sum()
            qc_summary.loc[i, "pct_missing"] = round(
                ((n_features_total - n_present) / n_features_total) * 100, 1,
            )

    qc_summary.to_csv(report_dir / "qc_summary.csv",
                       sep=";", encoding="utf-8-sig", index=False)

    # --- B) TIC overlay ---
    if save_plots:
        logger.info("  Generando TIC overlay...")
        plot_tic_overlay(all_gcms_data, save_path=report_dir / "tic_overlay.png")
        # Version normalizada
        plot_tic_overlay(all_gcms_data, save_path=report_dir / "tic_overlay_normalized.png",
                         normalize=True)

    # --- C) PCA QC TIC binned ---
    logger.info("  Calculando PCA QC (TIC binned)...")
    tic_matrix, bin_centers = build_tic_matrix(all_gcms_data)
    pca_result = pca_tic_qc(tic_matrix)
    if save_plots and pca_result is not None:
        plot_pca_tic(pca_result, save_path=report_dir / "pca_tic_qc.png")

    # --- D) Match score histogram ---
    if save_plots and not all_ids.empty:
        logger.info("  Generando histograma de match scores...")
        plot_match_score_histogram(all_ids, save_path=report_dir / "match_scores.png")

    # --- E) Heatmap features ---
    if save_plots and not data_matrix.empty:
        logger.info("  Generando heatmap de features...")
        plot_data_matrix_heatmap(data_matrix, save_path=report_dir / "heatmap_features.png")

    # --- F) PCA features (si hay >=2 muestras y >=2 features) ---
    pca_features_result = None
    if len(sample_names) >= 2 and not data_matrix.empty and data_matrix.shape[1] >= 2:
        logger.info("  Calculando PCA de features...")
        # Preparar matriz: log1p + fillna(0)
        feature_matrix_log = np.log1p(data_matrix.fillna(0))
        feature_df = pd.DataFrame(
            feature_matrix_log.values,
            index=data_matrix.index,
            columns=data_matrix.columns,
        )
        pca_features_result = pca_tic_qc(feature_df)
        if save_plots and pca_features_result is not None:
            plot_pca_tic(pca_features_result, save_path=report_dir / "pca_features.png")

    # --- G) Sample summary plot ---
    # Construir sample_stats compatible con plot_sample_summary
    sample_stats_rows = []
    for name, sr in sample_results.items():
        n_peaks = len(sr.peak_table)
        n_identified = sum(1 for m in sr.matches if m.best_match is not None)
        mean_int = float(sr.peak_table["intensity"].mean()) if n_peaks > 0 else 0
        max_int = float(sr.peak_table["intensity"].max()) if n_peaks > 0 else 0

        pct_missing = 0.0
        if n_features_total > 0 and name in data_matrix.index:
            n_present = data_matrix.loc[name].notna().sum()
            pct_missing = ((n_features_total - n_present) / n_features_total) * 100

        sample_stats_rows.append({
            "sample": name,
            "n_peaks": n_peaks,
            "n_identified": n_identified,
            "mean_intensity": mean_int,
            "max_intensity": max_int,
            "pct_missing": round(pct_missing, 1),
        })

    sample_stats = pd.DataFrame(sample_stats_rows)
    sample_stats.to_csv(report_dir / "sample_stats.csv",
                        sep=";", encoding="utf-8-sig", index=False)
    export_paths["sample_stats"] = report_dir / "sample_stats.csv"

    if save_plots and len(sample_stats) > 1:
        plot_sample_summary(sample_stats, save_path=report_dir / "sample_summary.png")

    # --- H) XLSX reporte multi-hoja ---
    logger.info("  Generando reporte XLSX multi-hoja...")
    report_tables = {"QC Summary": qc_summary}
    if not all_ids.empty:
        report_tables["Identifications"] = all_ids
    if not data_matrix.empty:
        report_tables["Data Matrix"] = data_matrix
        report_tables["Feature Metadata"] = feature_metadata
    report_tables["Sample Stats"] = sample_stats

    xlsx_report_path = export_xlsx_report(report_dir, report_tables)
    export_paths["xlsx_report"] = xlsx_report_path

    # --- Resumen final ---
    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETADO")
    logger.info("  Muestras procesadas: %d", len(sample_results))
    logger.info("  Features alineados: %d", len(features))
    logger.info("  Dimensiones de la matriz: %s", data_matrix.shape)
    logger.info("  Reporte final en: %s", report_dir)
    logger.info("  Archivos generados:")
    for f in sorted(report_dir.glob("*")):
        logger.info("    %s", f.name)
    logger.info("=" * 60)

    return {
        "sample_results": sample_results,
        "gcms_data": all_gcms_data,
        "features": features,
        "data_matrix": data_matrix,
        "feature_metadata": feature_metadata,
        "sample_stats": sample_stats,
        "qc_summary": qc_summary,
        "pca_tic": pca_result,
        "pca_features": pca_features_result,
        "export_paths": export_paths,
    }


def list_cdf_files(data_dir: str | Path | None = None) -> list[Path]:
    """Lista todos los archivos .cdf en un directorio."""
    data_dir = Path(data_dir) if data_dir else config.DATA_DIR
    files = sorted(data_dir.glob("*.CDF")) + sorted(data_dir.glob("*.cdf"))
    seen = set()
    unique = []
    for f in files:
        key = str(f).lower()
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique


if __name__ == "__main__":
    setup_logging()
    results = run_full_pipeline()
