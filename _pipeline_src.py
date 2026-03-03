# ═══════════════════════════════════════════════════════════════════════════════
# PIPELINE COMPLETO — módulos 1-6 para TODAS las muestras seleccionadas
# Se ejecuta automáticamente al pulsar "▶ Confirmar" en el selector de arriba.
# También se puede ejecutar manualmente (selected_cdf_paths debe estar definido).
# ═══════════════════════════════════════════════════════════════════════════════
import traceback
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import IPython, config
from pathlib import Path

from pipeline import run_single_sample
from matrix_builder import (
    SampleResult, align_peaks,
    build_data_matrix, build_feature_metadata,
    export_matrix, export_xlsx_report,
)
import qc


def _run_pipeline():
    ns = IPython.get_ipython().user_ns

    # ─── Comprobar que hay archivos seleccionados ─────────────────────────────
    selected_cdf_paths = ns.get('selected_cdf_paths', [])
    if not selected_cdf_paths:
        print("ℹ Sin archivos .cdf seleccionados. Ejecuta primero la celda del selector.")
        return

    print(f"▶ Procesando {len(selected_cdf_paths)} archivo(s)...")

    # ─── Directorios de salida ────────────────────────────────────────────────
    output_dir  = Path(config.OUTPUT_DIR)
    report_dir  = output_dir / "report"
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # ─── Cargar librería NIST (una vez) ──────────────────────────────────────
    nist_library = ns.get('_nist_library_cache', None)
    if nist_library is None:
        print("  Cargando librería NIST…")
        try:
            from nist_search import load_nist_library
            nist_library = load_nist_library(config.NIST_MSP_PATH)
            ns['_nist_library_cache'] = nist_library
            print(f"  Librería NIST cargada: {len(nist_library)} espectros")
        except Exception as exc:
            print(f"  ⚠ No se pudo cargar la librería NIST: {exc}")
            nist_library = []

    # ═════════════════════════════════════════════════════════════════════════
    # PASO 1 — Procesar cada muestra
    # ═════════════════════════════════════════════════════════════════════════
    print("\n─── PASO 1: Procesamiento por muestra ───")
    sample_results = {}
    gcms_data_dict = {}

    for cdf_path in selected_cdf_paths:
        name = Path(cdf_path).stem
        print(f"  [{name}]")
        try:
            gcms_data, peaks, spectra, matches, _id_tbl = run_single_sample(
                cdf_path, nist_library
            )
            result = SampleResult(
                name=name,
                peak_table=peaks.table,
                spectra=spectra,
                matches=matches,
            )
            sample_results[name] = result
            gcms_data_dict[name] = gcms_data
            print(f"    {peaks.n_peaks} picos detectados")
        except Exception:
            print(f"    ✗ Error procesando {name}:")
            traceback.print_exc()

    if not sample_results:
        print("✗ Ninguna muestra procesó correctamente. Abortando.")
        return

    # ─── Inferir metadatos de muestras ────────────────────────────────────────
    sample_names = list(sample_results.keys())
    sample_meta  = qc.infer_sample_metadata(sample_names)
    print(f"\n  Tipos de muestra detectados:")
    print(sample_meta[["sample_id", "sample_type"]].to_string(index=False))

    # ═════════════════════════════════════════════════════════════════════════
    # PASO 1.5 — PCA de TIC (control de calidad técnico)
    # ═════════════════════════════════════════════════════════════════════════
    pca_tic_result = None
    if getattr(config, 'ENABLE_PCA_TIC', True) and len(gcms_data_dict) >= 2:
        print("\n─── PASO 1.5: PCA de TIC (QC técnico) ───")
        try:
            tic_matrix, _ = qc.build_tic_matrix(gcms_data_dict)
            ok, msg = qc.validate_minimum_data(
                len(tic_matrix), tic_matrix.shape[1],
                min_samples=getattr(config, 'MIN_SAMPLES_FOR_PCA', 3),
                min_features=2,
                name="PCA TIC",
            )
            if ok:
                pca_tic_result = qc.pca_tic_qc(tic_matrix)
                fig_tic = qc.plot_pca_tic(
                    pca_tic_result,
                    sample_metadata=sample_meta,
                    save_path=figures_dir / "pca_tic_qc.png",
                )
                if fig_tic:
                    plt.close(fig_tic)
                    print(f"  PCA TIC guardada → {figures_dir / 'pca_tic_qc.png'}")
            else:
                print(f"  ⚠ {msg}")
        except Exception:
            print("  ⚠ Error en PCA TIC:")
            traceback.print_exc()

    # ═════════════════════════════════════════════════════════════════════════
    # PASO 2 — Alineamiento de picos y construcción de la matriz
    # ═════════════════════════════════════════════════════════════════════════
    print("\n─── PASO 2: Alineamiento y construcción de la matriz ───")
    features     = align_peaks(sample_results)
    data_matrix  = build_data_matrix(features, sample_names)
    feat_meta_df = build_feature_metadata(features)

    print(f"  Matriz: {data_matrix.shape[0]} muestras × {data_matrix.shape[1]} features")
    print(f"  Features con ID NIST: {(feat_meta_df['best_match_score'] > 0).sum()}")

    export_matrix(data_matrix, feat_meta_df, output_dir=output_dir, prefix="autogcms")

    # ═════════════════════════════════════════════════════════════════════════
    # PASO 3 — Preprocesado estadístico
    # ═════════════════════════════════════════════════════════════════════════
    print("\n─── PASO 3: Preprocesado estadístico ───")
    X_raw = data_matrix.copy()

    ok_min, msg_min = qc.validate_minimum_data(
        X_raw.shape[0], X_raw.shape[1],
        min_samples=getattr(config, 'MIN_SAMPLES_FOR_PCA', 3),
        min_features=getattr(config, 'MIN_FEATURES_FOR_ANALYSIS', 2),
        name="análisis de features",
    )
    if not ok_min:
        print(f"  ⚠ {msg_min}")
        X_proc    = X_raw.copy()
        kept_cols = list(X_raw.columns)
    else:
        X_proc, kept_cols = qc.preprocess_feature_matrix(X_raw, sample_metadata=sample_meta)
        print(f"  Features tras preprocesado: {len(kept_cols)} / {X_raw.shape[1]}")

    # ─── PASO 3a — QC summary ────────────────────────────────────────────────
    print("\n─── PASO 3a: QC summary ───")
    qc_basic = qc.compute_qc_summary(sample_results, gcms_data_dict)
    qc_adv   = qc.compute_advanced_qc_summary(X_raw, X_proc, sample_metadata=sample_meta)
    print(qc_basic.to_string(index=False))

    # ─── PASO 3b — PCA de features ───────────────────────────────────────────
    pca_feat_result = None
    if getattr(config, 'ENABLE_PCA_FEATURES', True) and ok_min:
        print("\n─── PASO 3b: PCA de features ───")
        try:
            ok_pca, msg_pca = qc.validate_minimum_data(
                X_proc.shape[0], X_proc.shape[1],
                min_samples=getattr(config, 'MIN_SAMPLES_FOR_PCA', 3),
                min_features=2,
                name="PCA features",
            )
            if ok_pca:
                pca_feat_result = qc.pca_features(X_proc)
                fig_feat = qc.plot_pca_features(
                    pca_feat_result,
                    sample_metadata=sample_meta,
                    save_path=figures_dir / "pca_features.png",
                )
                if fig_feat:
                    plt.close(fig_feat)
                    print(f"  PCA features guardada → {figures_dir / 'pca_features.png'}")
            else:
                print(f"  ⚠ {msg_pca}")
        except Exception:
            print("  ⚠ Error en PCA features:")
            traceback.print_exc()

    # ─── PASO 3c — Tabla de loadings ─────────────────────────────────────────
    loadings_df = pd.DataFrame()
    if pca_feat_result is not None:
        top_n_load  = getattr(config, 'PCA_TOP_LOADINGS_N', 20)
        loadings_df = qc.compute_pca_loadings_table(
            pca_feat_result, feat_meta=feat_meta_df, top_n=top_n_load
        )
        if not loadings_df.empty:
            print(f"\n  Top {top_n_load} loadings (PC1+PC2):")
            print(loadings_df.head(10).to_string(index=False))

    # ─── PASO 3d — Heatmap principal (log-normalizado) ───────────────────────
    if getattr(config, 'ENABLE_HEATMAP_MAIN', True) and ok_min:
        print("\n─── PASO 3d: Heatmap principal ───")
        try:
            ok_hm, msg_hm = qc.validate_minimum_data(
                X_proc.shape[0], X_proc.shape[1],
                min_samples=getattr(config, 'MIN_SAMPLES_FOR_HEATMAP', 2),
                min_features=2,
                name="heatmap principal",
            )
            if ok_hm:
                fig_hm = qc.plot_heatmap_main(
                    X_proc,
                    sample_metadata=sample_meta,
                    top_n=getattr(config, 'HEATMAP_TOP_N', 50),
                    save_path=figures_dir / "heatmap_main.png",
                )
                if fig_hm:
                    plt.close(fig_hm)
                    print(f"  Heatmap principal guardado → {figures_dir / 'heatmap_main.png'}")
            else:
                print(f"  ⚠ {msg_hm}")
        except Exception:
            print("  ⚠ Error en heatmap principal:")
            traceback.print_exc()

    # ─── PASO 3e — Heatmap z-score exploratorio ──────────────────────────────
    if getattr(config, 'ENABLE_HEATMAP_ZSCORE', True) and ok_min:
        print("\n─── PASO 3e: Heatmap z-score ───")
        try:
            ok_hz, msg_hz = qc.validate_minimum_data(
                X_proc.shape[0], X_proc.shape[1],
                min_samples=getattr(config, 'MIN_SAMPLES_FOR_HEATMAP', 2),
                min_features=2,
                name="heatmap z-score",
            )
            if ok_hz:
                fig_hz = qc.plot_heatmap_features(
                    X_proc,
                    top_n=getattr(config, 'HEATMAP_TOP_N', 50),
                    sample_metadata=sample_meta,
                    save_path=figures_dir / "heatmap_zscore.png",
                )
                if fig_hz:
                    plt.close(fig_hz)
                    print(f"  Heatmap z-score guardado → {figures_dir / 'heatmap_zscore.png'}")
            else:
                print(f"  ⚠ {msg_hz}")
        except Exception:
            print("  ⚠ Error en heatmap z-score:")
            traceback.print_exc()

    # ═════════════════════════════════════════════════════════════════════════
    # PASO 4 — Exportar reporte XLSX consolidado
    # ═════════════════════════════════════════════════════════════════════════
    print("\n─── PASO 4: Exportar reporte XLSX ───")
    try:
        tables = {
            "Data Matrix":      X_raw,
            "Data Matrix Proc": X_proc,
            "Feature Metadata": feat_meta_df,
            "Sample Metadata":  sample_meta,
            "QC Summary Basic": qc_basic,
            "QC Summary Adv":   qc_adv,
        }
        if not loadings_df.empty:
            tables["PCA Loadings"] = loadings_df

        xlsx_path = export_xlsx_report(report_dir, tables, "autogcms_report.xlsx")
        print(f"  Reporte guardado → {xlsx_path}")

        if ok_min and X_proc.shape[1] > 0:
            top_n_list = getattr(config, 'HEATMAP_TOP_N', 50)
            variances  = X_proc.var(axis=0)
            top_feats  = variances.nlargest(min(top_n_list, len(variances))).index.tolist()
            hm_feat_df = feat_meta_df[feat_meta_df["compound_name"].isin(top_feats)].copy()
            hm_feat_path = report_dir / "heatmap_feature_list.xlsx"
            hm_feat_df.to_excel(hm_feat_path, index=False)
            print(f"  Lista de features del heatmap → {hm_feat_path}")

    except Exception:
        print("  ⚠ Error exportando reporte XLSX:")
        traceback.print_exc()

    # ─── Resumen final ────────────────────────────────────────────────────────
    print("\n✔ Pipeline completado.")
    print(f"  Muestras procesadas : {len(sample_results)}")
    print(f"  Features alineadas  : {len(features)}")
    print(f"  Features (proc.)    : {len(kept_cols)}")
    print(f"  Salida              : {output_dir}")


# ── Punto de entrada ──────────────────────────────────────────────────────────
_run_pipeline()
