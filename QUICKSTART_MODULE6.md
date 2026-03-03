# GUÍA RÁPIDA - MÓDULO 6 MEJORADO

## Inicio Rápido

### 1. Ejecutar Pipeline Completo

En el notebook `run_pipeline.ipynb`:

```
Celda 1-5: Setup e inspección de datos (igual que antes)
Celda 6-12: Procesamiento de una muestra individual (igual que antes)
Celda 13: ⭐ NUEVO - Pipeline Módulo 6 completo
```

**Simplemente ejecuta la celda 13** y obtendrás:
- PCA de features preprocesadas
- Heatmap con clustering automático  
- Matrices raw/processed exportadas
- QC summary y reporte final

### 2. Personalizar en config.py

Edita los parámetros antes de ejecutar celda 13:

```python
# ACTIVAR/DESACTIVAR ANÁLISIS
ENABLE_FEATURE_PCA = True          # Quiero PCA
ENABLE_FEATURE_HEATMAP = True      # Quiero heatmap

# CANTIDAD DE FEATURES
HEATMAP_TOP_N = 50                 # Mostrar top 50 features más variables
PCA_TOP_LOADINGS_N = 20            # Mostrar top 20 en loadings

# PREPROCESADO (selecciona UNA de las opciones por línea)
NORMALIZATION_METHOD = "tic"       # ["none", "tic", "total_sum"]
LOG_TRANSFORM = "log10"             # ["none", "log10", "log2"]
SCALING_METHOD = "pareto"           # ["none", "pareto", "autoscale"]

# FILTROS
ENABLE_BLANK_FILTER = True          # Filtrar por ratio muestra/blank
BLANK_RATIO_THRESHOLD = 3.0         # ratio mínimo

ENABLE_QC_RSD_FILTER = True         # Filtrar por variabilidad en QCs
QC_RSD_THRESHOLD = 30.0             # RSD máximo (%)
```

### 3. Salidas en `final_report/`

Después de ejecutar celda 13, aparecen estos archivos:

| Archivo | Contenido |
|---------|-----------|
| `pca_features_scores.png` | Gráfico PC1 vs PC2 (coloreado por tipo muestra) |
| `pca_loadings_top_features.png` | Biplot con top features que contribuyen |
| `pca_explained_variance.csv` | % varianza explicada por cada PC |
| `heatmap_topN_features.png` | Heatmap con clustering jerárquico |
| `heatmap_feature_list.xlsx` | Lista de top features del heatmap |
| `qc_summary.xlsx` | Antes/después: n_features, missing %, sumas |
| `feature_matrix_raw.xlsx` | Matriz muestras × features (original) |
| `feature_matrix_processed.xlsx` | Matriz muestras × features (postprocesada) |
| `sample_metadata_used.xlsx` | Tipo de muestra inferido (Sample/QC/Blank) |
| `preprocessing_summary.txt` | Detalles del pipeline aplicado |

---

## Ejemplos de Diferencias de Configuración

### Opción A: Metabolómica conservadora  
```python
NORMALIZATION_METHOD = "tic"
LOG_TRANSFORM = "log10"
SCALING_METHOD = "autoscale"    # Estandarización clásica
ENABLE_BLANK_FILTER = True
BLANK_RATIO_THRESHOLD = 5.0     # Más estricto con blanks
QC_RSD_THRESHOLD = 25.0         # Más estricto con QCs
```

### Opción B: Descubrimiento exploratorio
```python
NORMALIZATION_METHOD = "none"   # Sin normalización
LOG_TRANSFORM = "none"           # Sin log
SCALING_METHOD = "pareto"        # Menos agresivo
ENABLE_BLANK_FILTER = False      # Mantener todos
ENABLE_QC_RSD_FILTER = False     # Confiar en todas
```

### Opción C: Recomendada (standard)
```python
# (Los valores default de config.py)
NORMALIZATION_METHOD = "tic"
LOG_TRANSFORM = "log10"
SCALING_METHOD = "pareto"
ENABLE_BLANK_FILTER = True
BLANK_RATIO_THRESHOLD = 3.0
ENABLE_QC_RSD_FILTER = True
QC_RSD_THRESHOLD = 30.0
```

---

## Interpretación de Resultados

### PCA Scores (pca_features_scores.png)
- **X axis (PC1)**: Primera dirección de máxima varianza
- **Y axis (PC2)**: Segunda dirección  
- **Colores**: 
  - 🔵 Blue = Samples (muestras reales)
  - 🟠 Orange = QC (control de calidad)
  - 🔴 Red = Blank (blancos)
- **Interpretación**: Clusters en el espacio indican muestras similares

### PCA Loadings (pca_loadings_top_features.png)
- **Flechas**: Dirección y magnitud de contribución de cada feature a PC1/PC2
- **Features cerca de PC1/PC2**: Contribuyen más a esa componente
- **Anotaciones**: Nombres de los top 20 features
- **Interpretación**: Identifica qué features diferencian las muestras

### Heatmap (heatmap_topN_features.png)
- **Filas**: Muestras
- **Columnas**: Top 50 features más variables
- **Colores**: Rojo (valores altos), Azul (valores bajos)
- **Dendrogramas**: Clustering jerárquico automático
- **Interpretación**: Identifica perfiles metabolómicos similares

### QC Summary (qc_summary.xlsx)
Compara antes vs después del preprocesado:
- `n_features_raw`: cuántos features (columnas) detectados
- `n_features_processed`: cuántos quedan después de filtros
- `pct_missing_raw`: % de valores faltantes originales
- **Interpretación**: Ver cuántos datos se perdieron en filtrado

---

## Cases de Uso

### Caso 1: Quiero ver qué features son importantes
→ Mira `pca_loadings_top_features.png`
→ Los que están más alejados del origen son los más importantes

### Caso 2: Quiero agrupar muestras similares
→ Mira `heatmap_topN_features.png`
→ El dendrograma en el eje Y muestra el clustering

### Caso 3: Tengo muy pocas features después de preproceso
→ Edita config.py: reduce `FEATURE_MAX_MISSING_FRAC`, `BLANK_RATIO_THRESHOLD`, o `QC_RSD_THRESHOLD`

### Caso 4: El heatmap tiene demasiadas features
→ Edita `HEATMAP_TOP_N` (e.g., 50 → 30)

### Caso 5: Quiero validar calidad técnica
→ Mira `qc_summary.xlsx`: % missing, n_features antes/después

---

## Debug: Si algo falla

Mira el output en la celda 13:

```
PASO 4: ...feature_metadata.xlsx con info de cada feature
  Metadata inferida:
    - Samples: 2
    - QCs: 0
    - Blanks: 0
```

Si ves:
- **Samples: 0** → Revisa nombres de archivo, podrían tener caracteres extraños
- **Features finales: 0** → Todos los features fueron filtrados (ajusta umbrales en config.py)
- **No se pudo completar PCA** → Pocas muestras (< 2) o features (< 1)

---

## Próximos Pasos Recomendados

1. **Validar identificaciones NIST** en `combined_identifications.csv`
   - ¿Los compuestos identificados tienen sentido?
   - ¿Son real hits o demos?

2. **Comparar resultados con software externo**
   - Exporta `feature_matrix_processed.xlsx` a Python/R/MetaboAnalyst
   - Reproduce PCA/heatmap para validar

3. **Mejorar metadata**
   - El sistema ahora infiere sample_type de nombres
   - TODO: Crear `sample_metadata.csv` formal si tienes info de grupos/tratamientos

4. **Personalizar análisis biológico**
   - Implementar contrastes entre grupos (si tienes repetidos)
   - Calcular p-values, fold-changes
   - Identificación de features biomarker

---

## Archivos Modificados Resumen

```
config.py               (+14 parámetros nuevos)
qc.py                   (+8 funciones nuevas)
run_pipeline.ipynb      (Celda 13 completamente reescrita)
```

✅ Módulos 1-5: **SIN CAMBIOS** (compatible 100%)

---

## Soporte

Si necesitas ajustar el análisis:

1. **Cambiar número de PCA components**: `config.PCA_N_COMPONENTS`
2. **Cambiar número de top loadings**: `config.PCA_TOP_LOADINGS_N`
3. **Cambiar número de features en heatmap**: `config.HEATMAP_TOP_N`
4. **Cambiar filtros estadísticos**: Ver tabla en config.py section "Análisis de Features"

