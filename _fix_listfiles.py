"""Actualiza la celda cell-list-files en el notebook con una versión robusta."""
import json
from pathlib import Path

nb_path = Path(r"c:\Users\LENOVO\Desktop\Claude-Proyectos\autogcms\notebooks\run_pipeline.ipynb")

new_src = """\
# Importar si no está disponible (p. ej. si se corre esta celda sola)
try:
    list_cdf_files
except NameError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent) if '__file__' in dir() else '.')
    try:
        from pipeline import list_cdf_files
        import config
    except Exception as _e:
        print(f"⚠ No se pudo importar pipeline: {_e}")
        print("Ejecuta primero la celda de imports (▶ la celda con 'import config...')")
        raise SystemExit(0)

import config as _cfg
_search_dir = _cfg.DATA_DIR
print(f"Buscando archivos .cdf en: {_search_dir}")

cdf_files = list_cdf_files()

if not cdf_files:
    print(f"  ⚠ No se encontraron archivos .cdf en {_search_dir}")
    print("     Comprueba config.DATA_DIR en config.py o reinicia el kernel y repite.")
    raise SystemExit(0)

print(f"Archivos .cdf disponibles ({len(cdf_files)}):\\n")
for i, f in enumerate(cdf_files):
    print(f"  [{i}] {f.name}")

# === SELECCIONAR AQUI ===
FILE_IDX = 0
# ========================

if FILE_IDX >= len(cdf_files):
    print(f"\\n⚠ FILE_IDX={FILE_IDX} fuera de rango (hay {len(cdf_files)} archivos, índices 0-{len(cdf_files)-1}).")
    raise SystemExit(0)

test_file = cdf_files[FILE_IDX]
print(f"\\nSeleccionado: {test_file.name}")
"""

lines = new_src.splitlines(keepends=True)

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

cell_found = False
for cell in nb["cells"]:
    if cell.get("id") == "cell-list-files":
        cell["source"] = lines
        cell["outputs"] = []
        cell["execution_count"] = None
        cell_found = True
        print(f"Celda cell-list-files actualizada ({len(lines)} líneas).")
        break

if not cell_found:
    print("ERROR: celda cell-list-files no encontrada.")
    exit(1)

with open(nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("Notebook guardado.")
