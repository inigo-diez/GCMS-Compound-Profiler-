# gcms-compound-profiler

---

## Sobre este proyecto

Este software surge como desarrollo complementario en el contexto del Trabajo de Fin de Grado (TFG) orientado al análisis metabolómico mediante cromatografía de gases acoplada a espectrometría de masas (GC-MS). El objetivo principal es **automatizar íntegramente el procesamiento de datos cromatográficos**, eliminando la necesidad de intervención manual en cada etapa del flujo de análisis.

El pipeline aborda tres problemas fundamentales del tratamiento de datos en GC-MS:

1. **Detección automática de picos cromatográficos** sobre el cromatograma de iones totales (TIC), incluyendo el manejo de compuestos co-eluidos mediante deconvolución gaussiana.

2. **Identificación automática de compuestos** mediante comparación espectral contra librerías de referencia en formato NIST, asignando a cada pico el candidato más probable con su puntuación de concordancia (Aún en proceso debido a la nula disponibilidad de las librerias oficial de compuestos)

3. **Análisis multivariante para la discriminación entre grupos experimentales**: a partir de la matriz de features alineadas entre muestras, se generan dos análisis de componentes principales (PCA) y dos mapas de calor (heatmaps) complementarios:

   - **PCA 1 — Control de calidad técnico (TIC binned):** se calcula sobre el perfil completo del cromatograma de iones totales, discretizado en ventanas de tiempo fijas. Su objetivo es detectar problemas instrumentales: muestras atípicas, efectos de lote o variabilidad en la inyección. No analiza compuestos individuales, sino la forma global del cromatograma.

   - **PCA 2 — Discriminación entre grupos (features alineadas):** se calcula sobre la matriz de intensidades de los compuestos detectados y alineados entre todas las muestras, previa normalización y transformación logarítmica. Es el análisis principal del experimento: permite visualizar si los distintos grupos experimentales se separan en el espacio multivariante y qué compuestos contribuyen más a dicha separación (loadings).

   - **Heatmap 1 — Intensidades log-normalizadas:** representa la intensidad de cada compuesto (en escala log₁₀) para cada muestra. Permite identificar visualmente qué metabolitos están presentes o ausentes en cada grupo y con qué abundancia relativa.

   - **Heatmap 2 — Z-score con clustering jerárquico:** aplica una estandarización por z-score sobre la misma matriz y agrupa muestras y compuestos por similitud de perfil mediante clustering jerárquico. Resulta especialmente útil para detectar patrones de co-expresión entre metabolitos y subgrupos de muestras que no son evidentes en la escala absoluta.

El conjunto de estas herramientas permite pasar de los archivos brutos del instrumento a resultados estadísticos interpretables de forma reproducible, trazable y sin intervención manual, con el fin de agilizar y estandarizar el tratamiento de datos en estudios de metabolómica basados en GC-MS.

> **Formato de entrada requerido: `.CDF` (Common Data Format / AIA)**
> El pipeline trabaja exclusivamente con archivos `.CDF`, que es el formato de exportación estándar e independiente del fabricante para instrumentos GC-MS. Estos archivos contienen el cromatograma de iones totales (TIC) y los espectros de masas completos adquiridos a lo largo del análisis. La práctica totalidad de los sistemas GC-MS comerciales (Agilent, Shimadzu, Thermo, Waters, etc.) permiten exportar los datos brutos a este formato desde su software propietario. El pipeline no es compatible con formatos nativos de fabricante (`.D`, `.RAW`, `.qgd`, etc.) ni con archivos de datos procesados.

---

Automated end-to-end pipeline for **GC-MS chromatogram analysis**: automatic peak detection, compound identification via spectral library matching, and multivariate group discrimination (PCA, heatmaps).

---

## Validated on real data

Results obtained processing **24 real GC-MS samples** (fecal metabolomics, TFG dataset):

| Metric | Value |
|---|---|
| Samples processed | 24 |
| Average peaks detected per sample | ~70 (range: 13–98) |
| Aligned features across all samples | 309 |
| Features retained after QC filtering | 45 (blank filter + missing value filter) |

---

## Example outputs

All figures below were generated automatically by the pipeline on a real fecal metabolomics dataset (20 samples).

---

### 1. Pipeline execution — sample loading

The notebook detects all `.CDF` files in the configured folder and begins processing them sequentially.

![Sample loading](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20202724.png)

---

### 2. Per-sample processing log

For each sample the pipeline reports: scans read, peaks detected, spectra extracted (including co-eluting peaks resolved by Gaussian deconvolution), and NIST identifications.

![Processing log 1](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20202812%20(2).png)
![Processing log 2](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20202857%20(3).png)
![Processing log 3](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20202932%20(4).png)
![Processing log 4](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20203005%20(5).png)
![Processing log 5](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20203048%20(6).png)

---

### 3. Final summary

Once all samples are processed, the pipeline prints the total peaks detected per sample and confirms the data is ready for multivariate analysis.

![Summary](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Captura%20de%20pantalla%202026-03-03%20203119%20(7).png)

---

### 4. Chromatogram preprocessing

Three-panel plot showing the signal processing applied to the raw TIC: (1) original signal, (2) Savitzky-Golay smoothing to reduce high-frequency noise, (3) AsLS baseline correction to remove the chemical background drift.

![Chromatogram preprocessing](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Visualizaci%C3%B3n%20cromatograma.png)

---

### 5. Automatic peak detection

The corrected TIC with all automatically detected peaks marked (red triangles). Each triangle indicates the apex of a detected chromatographic peak, labelled with its retention time in minutes.

![Peak detection](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Espectro%20de%20masas%20de%20uno%20de%20los%20cromatogramas.png)

---

### 6. Top 20 peaks by intensity

Table of the 20 most intense peaks for a given sample, including retention time, intensity, signal-to-noise ratio (S/N), and peak width. This table is exported automatically as a `.csv` for each sample.

![Top 20 peaks](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Top%2020%20picos%20de%20mayor%20intensidad%20de%20uno%20de%20los%20cromatogramas%20analizados.png)

---

### 7. Peak quality control

Four-panel QC plot per sample: intensity distribution, S/N ratio distribution (with detection threshold), peak width distribution, and scatter plot of RT vs. intensity. Used to assess the overall quality of the peak detection result.

![Peak QC](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Analaisis%20del%20control%20de%20calidad%20de%20uno%20de%20los%20cromatogramas.png)

---

### 8. Mass spectrum of an individual peak

For each detected peak, the pipeline extracts its full mass spectrum (m/z vs. intensity). This spectrum is then matched against the NIST reference library to identify the compound. The top fragments are reported for manual inspection.

![Mass spectrum](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Espectro%20de%20masas%20de%20uno%20de%20los%20picos.png)

---

### 9. PCA — Technical QC (TIC binned)

PCA computed on the full chromatogram profile of all samples (binned TIC). Detects batch effects, outlier injections, or instrument drift without relying on individual compound identification.

![PCA TIC QC](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/pca_tic_qc.png)

---

### 10. PCA — Group discrimination (aligned features)

PCA computed on the aligned feature matrix (log-normalized compound intensities across all samples). This is the main biological analysis: it reveals whether the experimental groups separate in multivariate space and which compounds drive the discrimination (loadings).

![PCA Features](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/pca_features.png)

---

### 11. Heatmap — Log-normalized intensities

Heatmap of the top variable features (log₁₀ intensity scale) across all samples. Rows are samples, columns are aligned features. Provides an immediate visual overview of which compounds are abundant or absent in each sample.

![Heatmap main](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/heatmap_main.png)

---

### 12. Heatmap — Z-score + hierarchical clustering

Same feature matrix after z-score standardization, with hierarchical clustering applied to both samples and features. Highlights relative differences between samples and reveals co-varying compound groups not visible in the absolute intensity scale.

![Heatmap z-score](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/heatmap_zscore.png)

---

## What does it do?

Given a set of raw GC-MS files (`.CDF` format), the pipeline:

1. **Reads and preprocesses** the chromatogram — smoothing + baseline correction
2. **Detects peaks** automatically in the Total Ion Chromatogram (TIC)
3. **Extracts mass spectra** at each peak, handling co-eluting compounds via Gaussian deconvolution
4. **Identifies compounds** by matching spectra against a reference library (NIST format `.msp`)
5. **Aligns peaks** across all samples and builds a feature × sample data matrix
6. **Generates QC reports and visualizations**: TIC overlays, PCA, heatmaps, match score histograms

---

## What are `.CDF` files?

`.CDF` (Common Data Format, also called AIA/ANDI format) is the **universal export standard for GC-MS instruments**. It stores the full raw data from each analysis run:
- The **Total Ion Chromatogram (TIC)**: signal intensity over time
- The **full mass spectrum** at each scan (mass-to-charge ratio vs. intensity)

Most GC-MS software (Agilent MassHunter, Shimadzu GCMS Solution, Thermo Xcalibur, etc.) can export to `.CDF`. It is instrument-independent — this pipeline works regardless of the GC-MS manufacturer.

---

## Requirements

### Software
- **Python 3.10+**
- **NIST MS Search** (optional but recommended for compound identification)
  - Required for matching spectra against the NIST/EPA/NIH Mass Spectral Library
  - If not available, the pipeline still runs but peaks will not be identified

### Python dependencies

```bash
pip install -r requirements.txt
```

| Package | Use |
|---|---|
| `netCDF4` | Reading `.CDF` files |
| `scipy` | Signal smoothing, baseline correction, peak detection |
| `numpy` | Numerical operations |
| `pandas` | Data tables and export |
| `matplotlib` / `seaborn` | All visualizations |
| `openpyxl` | Excel report export |
| `jupyter` | Running the notebook |
| `matchms` | Fallback spectral matching (if NIST MS Search is unavailable) |

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/inigo-diez/gcms-compound-profiler.git
cd gcms-compound-profiler
pip install -r requirements.txt
```

### 2. Configure paths in `config.py`

Open `config.py` and set:

```python
# Folder containing your .CDF files
DATA_DIR = Path(r"C:\path\to\your\CDF_files")

# Output folder (created automatically)
OUTPUT_DIR = Path(__file__).parent / "output"

# Path to NIST MS Search installation (optional)
NIST_DIR = Path(r"C:\path\to\MSSEARCH")
```

Point `DATA_DIR` to the folder where your `.CDF` files are stored. The pipeline will automatically detect and process all `.CDF` files in that folder.

![Selecting CDF files](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/Seleccionar%20.cdf.png)

### 3. Run the pipeline

**Option A — Streamlit web app (recommended):**

```bash
pip install streamlit
streamlit run app.py
```

Opens automatically in your browser. Select your `.CDF` folder, choose samples, and click **Run Pipeline**. All results appear inline: preprocessing plots, peak detection, identifications, PCA and heatmaps. Download buttons for CSV/Excel at the end.

**Option B — Jupyter Notebook:**

Open `notebooks/run_pipeline_BUENO.ipynb` and execute the cells step by step.

**Option C — Python script:**

```python
from pipeline import run_full_pipeline, setup_logging

setup_logging()
results = run_full_pipeline()
```

---

## Pipeline modules

| Module | File | Description |
|---|---|---|
| 1 — Preprocessing | `preprocessing.py` | Reads `.CDF`, applies Savitzky-Golay smoothing and AsLS baseline correction |
| 2 — Peak detection | `peak_detection.py` | Detects peaks in TIC using prominence, width, distance and S/N filters |
| 3 — Spectrum extraction | `spectra.py` | Extracts mass spectra at peak apexes; detects co-eluting compounds via Gaussian deconvolution; exports `.msp` |
| 4 — Identification | `nist_search.py` | Matches spectra against NIST library via cosine similarity; reports top candidates per peak |
| 5 — Data matrix | `matrix_builder.py` | Aligns peaks across samples by retention time; builds feature × sample intensity matrix |
| 6 — QC & visualization | `visualization.py` / `qc.py` | Generates all plots and the final Excel report |

---

## Outputs

All results are saved in the `output/` folder.

### Per sample
| File | Description |
|---|---|
| `<sample>_peaks.csv` | Peak table: retention time, intensity, S/N, width |
| `<sample>_identifications.csv` | Top compound candidates per peak (name, formula, CAS, match score) |
| `<sample>_spectra.msp` | Extracted mass spectra in NIST format |
| `<sample>_preprocessing.png` | Raw TIC vs. smoothed vs. baseline-corrected |
| `<sample>_peaks.png` | Detected peaks overlaid on TIC |
| `<sample>_peak_quality.png` | S/N and prominence distributions |

### Global (all samples)
| File | Description |
|---|---|
| `autogcms_data_matrix.csv` | Feature × sample intensity matrix (aligned peaks) |
| `autogcms_feature_metadata.csv` | Feature table: mean RT, putative compound ID, detection frequency |
| `autogcms_all_identifications.csv` | Combined identification table for all samples |

### Final report (`output/final_report/`)
| File | Description |
|---|---|
| `tic_overlay.png` | TIC of all samples overlaid — visual QC for run stability |
| `tic_overlay_normalized.png` | Same, normalized to peak area |
| `pca_tic_qc.png` | PCA of binned TICs — detects batch effects and outlier samples |
| `pca_features.png` | PCA of aligned features — group separation and sample clustering |
| `heatmap_features.png` | Heatmap of top variable features across samples |
| `match_scores.png` | Distribution of NIST match scores |
| `sample_summary.png` | Per-sample bar charts: n peaks, n identified, % missing |
| `qc_summary.csv` | QC metrics table (n peaks, n identified, % missing per sample) |
| `sample_stats.csv` | Summary statistics per sample |
| `autogcms_report.xlsx` | Multi-sheet Excel: QC Summary, Identifications, Data Matrix, Feature Metadata, Sample Stats |

---

## Key parameters (`config.py`)

| Parameter | Default | Description |
|---|---|---|
| `PEAK_MIN_PROMINENCE` | 5000 | Minimum peak prominence (adjust to your signal scale) |
| `PEAK_MIN_DISTANCE` | 10 | Minimum distance between peaks (in scans) |
| `PEAK_SNR_THRESHOLD` | 3.0 | Minimum signal-to-noise ratio |
| `NIST_MATCH_THRESHOLD` | 600 | Minimum match score (0–999) for compound identification |
| `RT_TOLERANCE_MIN` | 0.05 | Retention time tolerance for peak alignment across samples (minutes) |
| `HEATMAP_TOP_N` | 50 | Number of most variable features shown in the heatmap |
| `NORMALIZATION_METHOD` | `"tic"` | Normalization before PCA: `"none"`, `"tic"`, `"total_sum"` |
| `SCALING_METHOD` | `"pareto"` | Scaling before PCA: `"none"`, `"pareto"`, `"autoscale"` |
| `ENABLE_BLANK_FILTER` | `True` | Remove features dominated by blank samples |

---

## Project structure

```
gcms-compound-profiler/
├── config.py               # All parameters — edit this file to configure the pipeline
├── preprocessing.py        # Module 1: CDF reading, smoothing, baseline correction
├── peak_detection.py       # Module 2: Automatic peak detection and S/N filtering
├── spectra.py              # Module 3: Spectrum extraction and .msp export
├── nist_search.py          # Module 4: NIST library matching
├── matrix_builder.py       # Module 5: RT alignment and data matrix
├── visualization.py        # Module 6: All plots and QC visualizations
├── qc.py                   # QC metrics and PCA analysis
├── pipeline.py             # Orchestrator: run_single_sample() and run_full_pipeline()
├── nist_api_search.py      # Alternative: NIST WebBook API search
├── nist_binary_reader.py   # Reader for NIST binary library format
├── requirements.txt
└── notebooks/
    └── run_pipeline_BUENO.ipynb   # Interactive notebook with step-by-step execution
```

---

## License

MIT License — free to use, modify and distribute with attribution.
