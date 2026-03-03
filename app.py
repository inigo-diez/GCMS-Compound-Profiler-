"""
GCMS Compound Profiler — Streamlit web app.

Run with:
    streamlit run app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

# ─── Path setup ─────────────────────────────────────────────────────────────
APP_DIR = Path(__file__).parent
sys.path.insert(0, str(APP_DIR))

# ─── Page config ────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="GCMS Compound Profiler",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Helpers ────────────────────────────────────────────────────────────────

def show_fig(fig: plt.Figure, caption: str = "") -> None:
    """Display a matplotlib figure in Streamlit and immediately close it."""
    st.pyplot(fig, use_container_width=True)
    if caption:
        st.caption(caption)
    plt.close(fig)


# ─── App ────────────────────────────────────────────────────────────────────

def main() -> None:
    st.title("🔬 GCMS Compound Profiler")
    st.markdown(
        "Automated GC-MS pipeline: **peak detection → compound identification → "
        "multivariate analysis** (PCA, heatmaps)."
    )

    # ── Sidebar ──────────────────────────────────────────────────────────────
    with st.sidebar:
        st.header("⚙️ Configuration")

        data_dir = st.text_input(
            "📁 Data directory (.CDF files)",
            value=str(Path.home() / "Desktop"),
            help="Absolute path to the folder containing your .CDF files",
        )

        nist_dir = st.text_input(
            "🔬 NIST directory (optional)",
            value="",
            help="Path to the NIST MSSEARCH folder. Leave empty to skip identification.",
        )

        st.divider()
        st.subheader("Peak detection")
        min_prominence = st.number_input("Min prominence", value=5000, min_value=100, step=500)
        snr_threshold  = st.number_input("S/N threshold", value=3.0, min_value=0.5, step=0.5)
        min_distance   = st.number_input("Min distance (scans)", value=10, min_value=1, step=1)

        st.divider()
        st.subheader("Matrix & visualization")
        rt_tolerance  = st.number_input("RT tolerance (min)", value=0.05, min_value=0.01,
                                         step=0.01, format="%.2f")
        heatmap_top_n = st.number_input("Heatmap — top N features", value=50, min_value=5, step=5)

    # ── Detect CDF files ──────────────────────────────────────────────────────
    data_path = Path(data_dir)
    if not data_path.exists():
        st.warning(f"Directory not found: `{data_dir}`")
        st.info("Enter the path to your folder containing `.CDF` files in the sidebar.")
        return

    cdf_files: list[Path] = []
    seen: set[str] = set()
    for f in sorted(list(data_path.glob("*.CDF")) + list(data_path.glob("*.cdf"))):
        key = str(f).lower()
        if key not in seen:
            seen.add(key)
            cdf_files.append(f)

    if not cdf_files:
        st.warning(f"No `.CDF` files found in `{data_dir}`")
        return

    st.success(f"Found **{len(cdf_files)}** `.CDF` file(s) in `{data_dir}`")

    selected_names = st.multiselect(
        "Select samples to process",
        options=[f.name for f in cdf_files],
        default=[f.name for f in cdf_files],
    )
    selected_files = [f for f in cdf_files if f.name in selected_names]

    if not selected_files:
        st.warning("Select at least one sample.")
        return

    if not st.button("▶ Run Pipeline", type="primary"):
        return

    # ── Update config dynamically ─────────────────────────────────────────────
    import config
    config.DATA_DIR            = data_path
    config.PEAK_MIN_PROMINENCE = float(min_prominence)
    config.PEAK_SNR_THRESHOLD  = float(snr_threshold)
    config.PEAK_MIN_DISTANCE   = int(min_distance)
    config.RT_TOLERANCE_MIN    = float(rt_tolerance)
    config.HEATMAP_TOP_N       = int(heatmap_top_n)
    config.OUTPUT_DIR          = APP_DIR / "output"
    config.OUTPUT_DIR.mkdir(exist_ok=True)
    if nist_dir:
        config.NIST_DIR      = Path(nist_dir)
        config.NIST_MSP_PATH = Path(nist_dir)

    # ── Import pipeline modules after config is set ───────────────────────────
    from pipeline import setup_logging
    from preprocessing import preprocess
    from peak_detection import find_and_filter_peaks
    from spectra import extract_all_spectra, export_msp
    from nist_search import load_nist_library, match_all_spectra, build_identification_table
    from matrix_builder import (
        SampleResult, align_peaks, build_data_matrix,
        build_feature_metadata, export_matrix,
    )
    from visualization import (
        plot_tic_preprocessing, plot_peaks, plot_peak_quality,
        plot_tic_overlay, plot_data_matrix_heatmap,
    )
    from qc import build_tic_matrix, pca_tic_qc, plot_pca_tic

    setup_logging()

    # ── Load NIST library ─────────────────────────────────────────────────────
    with st.spinner("Loading reference library..."):
        library = load_nist_library(config.NIST_MSP_PATH if nist_dir else None)

    if library:
        st.info(f"Reference library loaded: **{len(library)}** compounds")
    else:
        st.warning("No NIST library loaded — peaks will be detected but not identified.")

    # ── Per-sample processing ──────────────────────────────────────────────────
    st.header("Per-sample results")

    sample_results: dict[str, SampleResult] = {}
    all_gcms_data: dict = {}
    all_id_tables: list[pd.DataFrame] = []

    progress = st.progress(0, text="Starting...")

    for i, cdf_path in enumerate(selected_files):
        progress.progress(i / len(selected_files), text=f"Processing {cdf_path.name} ...")

        with st.expander(f"📊 {cdf_path.stem}", expanded=(i == 0)):

            # Module 1 — Preprocessing
            with st.spinner("Preprocessing..."):
                data = preprocess(cdf_path)

            # Module 2 — Peak detection
            peaks = find_and_filter_peaks(data)

            # Module 3 — Spectrum extraction
            spectra = extract_all_spectra(data, peaks)
            export_msp(spectra, config.OUTPUT_DIR / f"{data.name}_spectra.msp",
                       sample_name=data.name)

            # Module 4 — NIST matching
            if spectra and library:
                matches  = match_all_spectra(spectra, library, sample_name=data.name)
                id_table = build_identification_table(spectra, matches)
            else:
                matches  = []
                id_table = pd.DataFrame()

            # Save CSVs
            peaks.table.to_csv(
                config.OUTPUT_DIR / f"{data.name}_peaks.csv",
                sep=";", encoding="utf-8-sig", index=False,
            )
            if not id_table.empty:
                id_table.to_csv(
                    config.OUTPUT_DIR / f"{data.name}_identifications.csv",
                    sep=";", encoding="utf-8-sig", index=False,
                )

            # Store results
            sample_results[data.name] = SampleResult(
                name=data.name,
                peak_table=peaks.table,
                spectra=spectra,
                matches=matches,
            )
            all_gcms_data[data.name] = data
            if not id_table.empty:
                id_copy = id_table.copy()
                id_copy.insert(0, "sample", data.name)
                all_id_tables.append(id_copy)

            # ─ Metrics ─
            n_id = (
                id_table.loc[id_table["rank"] == 1, "compound_name"].ne("No match").sum()
                if not id_table.empty else 0
            )
            c1, c2, c3 = st.columns(3)
            c1.metric("Peaks detected",    peaks.n_peaks)
            c2.metric("Spectra extracted", len(spectra))
            c3.metric("Identified",        int(n_id))

            # ─ Plots ─
            col_a, col_b = st.columns(2)
            with col_a:
                show_fig(
                    plot_tic_preprocessing(data, save_path=None),
                    "Preprocessing: raw → smoothed → baseline corrected",
                )
            with col_b:
                show_fig(
                    plot_peaks(data, peaks, save_path=None),
                    "Detected peaks (red triangles = peak apexes)",
                )

            show_fig(
                plot_peak_quality(peaks, save_path=None),
                "Peak QC: intensity distribution, S/N, width, RT vs intensity",
            )

            # ─ Tables ─
            st.subheader("Top 20 peaks by intensity")
            st.dataframe(peaks.table.head(20), use_container_width=True)

            if not id_table.empty:
                st.subheader("Top identifications")
                top_ids = (
                    id_table[id_table["rank"] == 1]
                    [["rt_min", "compound_name", "formula", "match_score"]]
                    .head(20)
                )
                st.dataframe(top_ids, use_container_width=True)

    progress.progress(1.0, text="All samples processed ✓")

    # ── Multivariate analysis (≥ 2 samples) ────────────────────────────────
    if len(sample_results) < 2:
        st.info("Process at least **2 samples** to enable multivariate analysis.")
        return

    st.header("Multivariate analysis")

    with st.spinner("Aligning peaks and building data matrix..."):
        features        = align_peaks(sample_results)
        sample_names    = list(sample_results.keys())
        data_matrix     = build_data_matrix(features, sample_names)
        feature_metadata = build_feature_metadata(features)
        export_matrix(data_matrix, feature_metadata, config.OUTPUT_DIR)

    st.success(
        f"Data matrix: **{data_matrix.shape[0]} samples × {data_matrix.shape[1]} features**"
    )

    # TIC overlay + PCA TIC
    col1, col2 = st.columns(2)
    with col1:
        show_fig(
            plot_tic_overlay(all_gcms_data, save_path=None),
            "TIC overlay — all samples superimposed",
        )
    with col2:
        tic_matrix, _ = build_tic_matrix(all_gcms_data)
        pca_tic = pca_tic_qc(tic_matrix)
        if pca_tic:
            show_fig(
                plot_pca_tic(pca_tic, save_path=None),
                "PCA 1 — Technical QC: batch effects and outlier samples",
            )

    # PCA features + Heatmap
    col3, col4 = st.columns(2)
    with col3:
        if data_matrix.shape[1] >= 2:
            feat_log = pd.DataFrame(
                np.log1p(data_matrix.fillna(0).values),
                index=data_matrix.index,
                columns=data_matrix.columns,
            )
            pca_feat = pca_tic_qc(feat_log)
            if pca_feat:
                show_fig(
                    plot_pca_tic(pca_feat, save_path=None),
                    "PCA 2 — Group discrimination: experimental group separation",
                )
    with col4:
        if not data_matrix.empty:
            show_fig(
                plot_data_matrix_heatmap(data_matrix, save_path=None),
                "Heatmap — log₁₀ intensity per sample and feature",
            )

    # ── Downloads ────────────────────────────────────────────────────────────
    st.header("⬇ Download results")

    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
        st.download_button(
            "Data matrix (.csv)",
            data=data_matrix.to_csv(sep=";").encode("utf-8-sig"),
            file_name="data_matrix.csv",
            mime="text/csv",
        )
    with col_d2:
        st.download_button(
            "Feature metadata (.csv)",
            data=feature_metadata.to_csv(sep=";", index=False).encode("utf-8-sig"),
            file_name="feature_metadata.csv",
            mime="text/csv",
        )
    if all_id_tables:
        all_ids = pd.concat(all_id_tables, ignore_index=True)
        with col_d3:
            st.download_button(
                "All identifications (.csv)",
                data=all_ids.to_csv(sep=";", index=False).encode("utf-8-sig"),
                file_name="identifications.csv",
                mime="text/csv",
            )


if __name__ == "__main__":
    main()
