"""
Configuración centralizada del pipeline AutoGCMS.

Modifica los parámetros aquí sin tocar el resto del código.
"""
from pathlib import Path

# ─────────────────────────────────────────────
# Rutas
# ─────────────────────────────────────────────
DATA_DIR = Path(r"C:\Users\LENOVO\Desktop\TFG\CDF")
OUTPUT_DIR = Path(__file__).parent / "output"
NIST_DIR = Path(r"C:\Users\LENOVO\OneDrive\Documentos\MATLAB\MSSEARCH")

# ─────────────────────────────────────────────
# Preprocesado (Módulo 1)
# ─────────────────────────────────────────────
SAVGOL_WINDOW = 15          # Ventana Savitzky-Golay (debe ser impar)
SAVGOL_POLYORDER = 3        # Orden del polinomio
BASELINE_LAMBDA = 1e6       # Suavidad de la corrección de baseline (AsLS)
BASELINE_P = 0.01           # Asimetría de la corrección de baseline

# ─────────────────────────────────────────────
# Detección de picos (Módulo 2)
# ─────────────────────────────────────────────
PEAK_MIN_PROMINENCE = 5000  # Prominencia mínima para detectar un pico
PEAK_MIN_DISTANCE = 10      # Distancia mínima entre picos (en scans)
PEAK_MIN_WIDTH = 3          # Ancho mínimo del pico (en scans)
PEAK_SNR_THRESHOLD = 3.0    # Ratio señal/ruido mínimo
MAX_PEAKS = None            # Máximo de picos a retener (None = sin límite)

# ─────────────────────────────────────────────
# Extracción de espectros (Módulo 3)
# ─────────────────────────────────────────────
APEX_SCANS_AVERAGE = 3      # Número de scans alrededor del apex para promediar
COELUTION_BIC_DELTA = 10    # Delta BIC para decidir co-elución (2 Gauss vs 1 Gauss)
SPECTRA_AVG_SCANS = 2       # Scans a promediar alrededor de cada subpico co-eluido

# ─────────────────────────────────────────────
# Identificación de compuestos (Módulo 4)
# ─────────────────────────────────────────────
# Ruta a la librería NIST exportada en formato .msp
# Opción A: apuntar a un archivo .msp concreto (exportado desde NIST MS Search)
# Opción B: apuntar a un directorio con archivos .MSP (se cargan todos)
# Si el archivo no existe, se buscan .MSP en NIST_DIR automáticamente.
NIST_MSP_PATH = NIST_DIR  # Busca todos los .MSP en el directorio MSSEARCH

NIST_MATCH_THRESHOLD = 600  # Match score mínimo (0-999), equivale a coseno ~0.6
NIST_TOP_N = 5              # Número de candidatos a devolver por pico (aumentado)
COSINE_THRESHOLD = 0.5      # Umbral de similitud coseno (0.0 - 1.0) - reducido para más matches

# ─────────────────────────────────────────────
# Alineamiento y matriz (Módulo 5)
# ─────────────────────────────────────────────
RT_TOLERANCE_MIN = 0.05     # Tolerancia de RT para alineamiento (minutos)

# ─────────────────────────────────────────────
# QC y PCA (Módulo 2 / 6)
# ─────────────────────────────────────────────
TIC_BIN_SIZE_MIN = 0.05     # Tamaño del bin en minutos para TIC binned (PCA QC)
PCA_N_COMPONENTS = 2        # Componentes PCA para QC

# ─────────────────────────────────────────────
# Análisis de Features (Módulo 6 - Nuevos)
# ─────────────────────────────────────────────
ENABLE_FEATURE_PCA = True                   # Activar PCA de features alineadas
ENABLE_FEATURE_HEATMAP = True               # Activar heatmap de features
HEATMAP_TOP_N = 50                          # N features más variables para heatmap

# Preprocesado estadístico:
FEATURE_MAX_MISSING_FRAC = 0.5              # Máx fracción de ceros/NA por feature (0-1)
IMPUTE_METHOD = "min_half"                  # ["min_half", "median"]
NORMALIZATION_METHOD = "tic"                # ["none", "tic", "total_sum"]
LOG_TRANSFORM = "log10"                     # ["none", "log10", "log2"]
SCALING_METHOD = "pareto"                   # ["none", "pareto", "autoscale"]

ENABLE_BLANK_FILTER = True                  # Filtrar features dominadas por blanks
BLANK_FILTER_RULE = "sample_to_blank_ratio" # Regla: "sample_to_blank_ratio" o "blank_intensity_ratio"
BLANK_RATIO_THRESHOLD = 3.0                 # Umbral de ratio muestra/blank

ENABLE_QC_RSD_FILTER = True                 # Filtrar por RSD si hay QCs
QC_RSD_THRESHOLD = 30.0                     # % RSD máximo en QCs para retener feature

# PCA visualización:
PCA_LABEL_POINTS = False                    # Anotar nombres de muestras en PCA
PCA_TOP_LOADINGS_N = 20                     # Top N features en diagrama de loadings

# ─────────────────────────────────────────────
# Controles de activación (Módulo 6)
# ─────────────────────────────────────────────
ENABLE_PCA_TIC      = True   # PCA de TIC binned (QC técnico)
ENABLE_PCA_FEATURES = True   # PCA de features alineadas (análisis principal)
ENABLE_HEATMAP_MAIN  = True  # Heatmap principal (valores log-normalizados)
ENABLE_HEATMAP_ZSCORE = True # Heatmap exploratorio (z-score + clustering)

# Umbrales mínimos de datos para análisis:
MIN_SAMPLES_FOR_PCA      = 3   # Muestras válidas mínimas para calcular PCA
MIN_SAMPLES_FOR_HEATMAP  = 2   # Muestras válidas mínimas para heatmap
MIN_FEATURES_FOR_ANALYSIS = 2  # Features mínimas tras filtrado

# Heatmap opciones:
HEATMAP_MAIN_USE_IDENTIFIED_ONLY = False  # Solo features con ID NIST en heatmap
HEATMAP_CLUSTERING = True                 # Clustering jerárquico (requiere scipy)

# ─────────────────────────────────────────────
# General
# ─────────────────────────────────────────────
LOGGING_LEVEL = "INFO"
