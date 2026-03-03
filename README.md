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

## Example outputs

| PCA 1 — Technical QC (TIC binned) | PCA 2 — Group discrimination (aligned features) |
|:---:|:---:|
| ![PCA TIC QC](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/pca_tic_qc.png) | ![PCA Features](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/pca_features.png) |

| Heatmap 1 — Log-normalized intensities | Heatmap 2 — Z-score + hierarchical clustering |
|:---:|:---:|
| ![Heatmap Main](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/heatmap_main.png) | ![Heatmap Z-score](https://raw.githubusercontent.com/inigo-diez/GCMS-Compound-Profiler-/main/figures/heatmap_zscore.png) |

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

### 3. Run the pipeline

**Option A — Jupyter Notebook (recommended):**

Open `notebooks/run_pipeline_BUENO.ipynb` and execute the cells step by step.

**Option B — Python script:**

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
