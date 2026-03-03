# MÓDULO 6 - MEJORAS DE ANÁLISIS

## Resumen de Cambios Implementados

### 1. **config.py** - Nuevos parámetros configurables

Se agregaron 14 parámetros nuevos para controlar análisis de features (Módulo 6):

| Parámetro | Valor default | Descripción |
|-----------|----------------|------------|
| `ENABLE_FEATURE_PCA` | `True` | Activar PCA de features alineadas |
| `ENABLE_FEATURE_HEATMAP` | `True` | Activar heatmap de features |
| `HEATMAP_TOP_N` | `50` | N features más variables para heatmap |
| `FEATURE_MAX_MISSING_FRAC` | `0.5` | Máx fracción de NaN/ceros por feature |
| `IMPUTE_METHOD` | `"min_half"` | Método imputación: ["min_half", "median"] |
| `NORMALIZATION_METHOD` | `"tic"` | Normalización: ["none", "tic", "total_sum"] |
| `LOG_TRANSFORM` | `"log10"` | Transform: ["none", "log10", "log2"] |
| `SCALING_METHOD` | `"pareto"` | Escalado: ["none", "pareto", "autoscale"] |
| `ENABLE_BLANK_FILTER` | `True` | Filtrar features dominadas por blanks |
| `BLANK_RATIO_THRESHOLD` | `3.0` | Ratio min sample/blank para retener |
| `ENABLE_QC_RSD_FILTER` | `True` | Filtrar features por RSD si hay QCs |
| `QC_RSD_THRESHOLD` | `30.0` | RSD máximo (%) en QCs |
| `PCA_LABEL_POINTS` | `False` | Anotar nombres en PCA |
| `PCA_TOP_LOADINGS_N` | `20` | Top N features en diagrama loadings |

---

### 2. **qc.py** - Nuevas funciones para Módulo 6

#### **Metadatos e Inferencia**
- `infer_sample_metadata(sample_names) -> DataFrame`
  - Infiere sample_type (Sample/QC/Blank/Control) desde nombres
  - Crea tabla metadata para usar en análisis

#### **Preprocesado Estadístico (Pipeline configurable)**
- `preprocess_feature_matrix(...) -> (DataFrame, dict)`
  - 7 pasos en orden: blank filter → missing value filter → imputation → normalization → log transform → scaling → QC RSD filter
  - Retorna matriz procesada + estadísticas de cada paso
  - Totalmente configurable vía config.py

#### **PCA de Features (ANÁLISIS PRINCIPAL)**
- `pca_features(matrix) -> dict`
  - PCA sobre matriz de features preprocesada
  - Devuelve scores, loadings, varianza explicada, nombres
  
- `plot_pca_features(pca_result, sample_metadata) -> Figure`
  - Scatter PC1 vs PC2
  - Coloreado por tipo de muestra (QC/blank/sample)
  - Guarda en PNG

- `plot_pca_loadings(pca_result, feature_metadata, top_n) -> Figure`
  - Biplot de loadings
  - Muestra top N features que contribuyen más
  - Incluye anotaciones con nombres de features

#### **Heatmap de Features (ANÁLISIS PRINCIPAL)**
- `plot_heatmap_features(matrix, top_n, sample_metadata) -> Figure`
  - Heatmap de top N features más variables
  - Z-score por columna para visualización
  - Clustering jerárquico automático si scipy disponible
  - Etiquetas legibles de muestras y features
  - Guarda en PNG

#### **QC Summary**
- `compute_advanced_qc_summary(...) -> DataFrame`
  - Compara antes/después: n_features, missing values, sum intensities
  - Por cada muestra
  - Exporta a XLSX

---

### 3. **run_pipeline.ipynb (Celda 13)** - Nuevo esquema Módulo 6

La celda del pipeline completo (antigua cel #VSC-45ec3898) fue completamente reescrita con estructura 7 pasos:

**PASO 1: Procesamiento muestras** (módulos 1-4)
- Itera sobre todos los .cdf
- Acumula results en `sample_results` dict

**PASO 2: Alineamiento picos** (módulo 5)
- `align_peaks()` agrupa picos por RT
- Construye matriz raw (muestras × features)

**PASO 3: QC técnico**
- PCA de TIC binned (marco para futuro uso)
- Requiere GCMSData objects (mejora futura)

**PASO 4: Metadatos y preprocesado**
- `infer_sample_metadata()` del nombre del archivo
- `preprocess_feature_matrix()` con pipeline 7 pasos
- Genera X_processed matriz lista para análisis

**PASO 5: PCA de features** (ANÁLISIS PRINCIPAL)
- `pca_features()` sobre X_processed
- Genera 2 plots: scores + loadings
- Exporta varianza explicada CSV

**PASO 6: Heatmap de features** (ANÁLISIS PRINCIPAL)
- `plot_heatmap_features()` con top N features
- Clustering automático
- Exporta lista de features usados XLSX

**PASO 7: QC summary y exports**
- `compute_advanced_qc_summary()` -> XLSX
- Preprocessing summary -> TXT
- Feature matrices (raw + processed) -> XLSX
- Sample metadata -> XLSX
- Combined identifications -> CSV

---

### 4. **Archivo de Salida (final_report/)**

La ejecución genera estos archivos nuevos:

```
final_report/
├── pca_features_scores.png          # PCA scores (PC1 vs PC2)
├── pca_features_loadings.png        # PCA biplot loadings
├── pca_explained_variance.csv       # % varianza por PC
├── heatmap_topN_features.png        # Heatmap con clustering
├── heatmap_feature_list.xlsx        # Lista de top N features
├── qc_summary.xlsx                  # QC antes/después preproceso
├── preprocessing_summary.txt        # Resumen pipeline estadístico
├── feature_matrix_raw.xlsx          # Matriz original (muestras × features)
├── feature_matrix_processed.xlsx    # Matriz postprocesada
└── sample_metadata_used.xlsx        # Metadata sample_type/group
```

---

## Separación Clara: PCA de QC vs. Análisis Químico

### ❌ PCA de TIC (QC técnico SOLAMENTE)
- Usa TIC binned normalizado
- PASO 3 en el pipeline
- Propósito: validar consistencia técnica entre muestras
- NO para análisis químico principal

### ✅ PCA de Features (ANÁLISIS PRINCIPAL)
- Usa matriz de features alineadas preprocesada
- PASO 5 en el pipeline
- Propósito: descubrimiento de patrones biológicos/químicos
- Este es el análisis principal

### ✅ Heatmap de Features (ANÁLISIS PRINCIPAL)
- Usa matriz de features alineadas preprocesada
- PASO 6 en el pipeline
- Propósito: exploración visual de perfiles metabolómicos
- Clustering jerárquico automático

---

## Pipeline Preprocesado (7 pasos)

Ejecutado en `preprocess_feature_matrix()`:

1. **Blank filtering** (si `ENABLE_BLANK_FILTER=True`)
   - Elimina features donde sample/blank ratio < BLANK_RATIO_THRESHOLD (default 3.0)

2. **Missing value filter**
   - Elimina features con > FEATURE_MAX_MISSING_FRAC (default 50%) de NaN/ceros

3. **Imputación**
   - `"min_half"`: reemplaza NaN con min(feature)/2
   - `"median"`: reemplaza con mediana

4. **Normalización**
   - `"none"`: sin cambios
   - `"tic"`: divide cada muestra por su suma (proporciones)
   - `"total_sum"`: normaliza por suma total del dataset

5. **Log transform**
   - `"none"`: sin cambios
   - `"log10"`: log₁₀(x + 1)
   - `"log2"`: log₂(x + 1)

6. **Scaling**
   - `"none"`: sin cambios
   - `"pareto"`: (x - mean) / √std
   - `"autoscale"`: (x - mean) / std

7. **QC RSD filter** (si `ENABLE_QC_RSD_FILTER=True` y hay QCs)
   - Elimina features con RSD > QC_RSD_THRESHOLD (default 30%) en muestras QC

---

## Metadata Inference

Función `infer_sample_metadata()` heurística:

```
Si nombre contiene:
  "qc"     → sample_type = "QC"
  "control" → sample_type = "QC"
  "blank"   → sample_type = "Blank"
  "blanco"  → sample_type = "Blank"
  
Else → sample_type = "Sample"
```

⚠️ TODO: Implementar lectura de metadata formal (CSV, JSON, etc) si existe

---

## Notas de Implementación

### ✅ Completado
- Todas las funciones de preprocesado estadístico
- PCA de features con loadings y varianza explicada
- Heatmap con clustering jerárquico automático
- QC summary before/after
- Exports matrices raw/processed
- Configurabilidad vía config.py
- Manejo de edge cases (no QCs, no blanks, matriz vacía)
- Logs claros de cada paso

### ⚠️ Mejoras Futuras
- Pasar GCMSData objects al pipeline para activar PCA QC de TIC
- Metadata formal desde archivo externo
- Anotación de identificaciones NIST en PCA loadings
- Validación de NIST hits en top loadings
- Export a formato compatible con software metabolómica (mzTab, etc)
- Heatmap con valores de p-value si contrastes disponibles

---

## Ejecución (Célula 13)

```python
# La célula corre todo automáticamente:
# 1. Lee todos los .cdf
# 2. Procesa módulos 1-4
# 3. Alinea picos (módulo 5)
# 4. Preprocesa estadísticamente
# 5. Genera PCA y heatmap
# 6. Exporta report completo

# Salida esperada: final_report/ con 10+ archivos
```

---

## Archivos Modificados

| Archivo | Cambios |
|---------|---------|
| `config.py` | +14 parámetros nuevos |
| `qc.py` | +8 funciones nuevas (~250 líneas) |
| `run_pipeline.ipynb` | Cell 13: reescrita completamente (7 pasos) |
| `matrix_builder.py` | Sin cambios (compatible) |
| Otros módulos | Sin cambios |

---

## Compatible con Versión Anterior

✅ Los módulos 1-5 NO fueron modificados
✅ Las funciones nuevas son aditivas (backward compatible)
✅ Config.py solo añade parámetros (valores default funcionales)
✅ Puede ejecutarse parcialmente (solo paso 1, por ejemplo)

