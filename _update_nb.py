"""Inyecta el contenido de _pipeline_src.py en la celda 3958fffb del notebook."""
import json
from pathlib import Path

nb_path  = Path(r"c:\Users\LENOVO\Desktop\Claude-Proyectos\autogcms\notebooks\run_pipeline.ipynb")
src_path = Path(r"c:\Users\LENOVO\Desktop\Claude-Proyectos\autogcms\_pipeline_src.py")

with open(src_path, "r", encoding="utf-8") as f:
    raw_src = f.read()

# Convertir a lista de líneas con \n al final (formato JSON del notebook)
lines = raw_src.splitlines(keepends=True)
# Asegurarse de que la última línea no tenga \n duplicado
# (el notebook format espera exactamente lo mismo que se vería en el editor)

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

cell_found = False
for cell in nb["cells"]:
    if cell.get("id") == "3958fffb":
        cell["source"] = lines
        cell["cell_type"] = "code"
        cell["outputs"] = []
        cell["execution_count"] = None
        cell_found = True
        print(f"Celda 3958fffb actualizada ({len(lines)} líneas, {len(raw_src)} caracteres).")
        break

if not cell_found:
    print("ERROR: celda 3958fffb no encontrada.")
    exit(1)

with open(nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("Notebook guardado.")

# Verificación
with open(nb_path, "r", encoding="utf-8") as f:
    nb2 = json.load(f)
for cell in nb2["cells"]:
    if cell.get("id") == "3958fffb":
        src_check = "".join(cell["source"])
        print(f"Verificación OK: {len(src_check)} caracteres en la celda.")
        break
