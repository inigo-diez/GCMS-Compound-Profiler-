"""Actualiza la celda selector b068f5f9: Browse añade directo a Selected."""
import json
from pathlib import Path

nb_path = Path(r"c:\Users\LENOVO\Desktop\Claude-Proyectos\autogcms\notebooks\run_pipeline.ipynb")

new_src = """\
# SELECTOR DE ARCHIVOS .CDF
# Selecciona los archivos → pulsa "Confirmar" → el pipeline arranca en la celda siguiente
from IPython.display import display, HTML, clear_output
from ipywidgets import SelectMultiple, Button, VBox, HBox, Label, Output
import IPython, pipeline, config
from pathlib import Path
import tkinter as tk
from tkinter import filedialog

original_list_cdf  = pipeline.list_cdf_files
selected_cdf_paths = []
out_run = Output()

# ── Widgets ──────────────────────────────────────────────────────────────────
existing_list = SelectMultiple(rows=8, description='', layout={'width': '340px'})
selected_list = SelectMultiple(rows=8, description='', layout={'width': '340px'})
browse_btn  = Button(description='+ Buscar archivos .cdf', button_style='info',    layout={'width': '200px'})
add_btn     = Button(description='Añadir →',               button_style='primary', layout={'width': '110px'})
remove_btn  = Button(description='← Quitar',               button_style='warning', layout={'width': '110px'})
confirm_btn = Button(description='▶ Confirmar y analizar', button_style='success', layout={'width': '200px'})
reset_btn   = Button(description='↺ Reiniciar',            button_style='',        layout={'width': '120px'})
status_lbl  = Label(value='')

def _refresh_existing():
    try:    paths = list(original_list_cdf())
    except: paths = []
    existing_list.options = [(p.name, str(p)) for p in paths]
_refresh_existing()

def _on_browse(b):
    \"\"\"Abre el explorador y añade los archivos elegidos DIRECTAMENTE a Selected.\"\"\"
    try:
        root = tk.Tk(); root.withdraw(); root.wm_attributes('-topmost', True)
        paths = filedialog.askopenfilenames(
            title='Selecciona archivos .cdf para analizar',
            filetypes=[('CDF files', '*.cdf *.CDF'), ('Todos', '*.*')])
        root.destroy()
    except Exception as e:
        status_lbl.value = f'Error explorador: {e}'; return

    if not paths:
        status_lbl.value = 'No se seleccionó ningún archivo.'; return

    # Añadir al panel Selected directamente
    curr_sel = dict(selected_list.options)
    added = 0
    for p in paths:
        src = Path(p)
        if str(src) not in curr_sel.values():
            curr_sel[src.name] = str(src)
            added += 1

    selected_list.options = list(curr_sel.items())
    status_lbl.value = f'{added} archivo(s) añadido(s) → {len(curr_sel)} en Selected.'

def _on_add(b):
    \"\"\"Mueve los items seleccionados en Existing → Selected.\"\"\"
    vals = list(existing_list.value)
    if not vals: status_lbl.value = 'Selecciona primero en la lista izquierda.'; return
    curr = dict(selected_list.options)
    for v in vals:
        if v not in curr.values(): curr[Path(v).name] = v
    selected_list.options = list(curr.items())
    status_lbl.value = f'{len(selected_list.options)} archivo(s) en Selected.'

def _on_remove(b):
    to_remove = set(selected_list.value)
    selected_list.options = [(n, v) for n, v in selected_list.options if v not in to_remove]
    status_lbl.value = f'{len(selected_list.options)} archivo(s) en Selected.'

def _on_reset(b):
    global selected_cdf_paths
    selected_cdf_paths = []
    pipeline.list_cdf_files = original_list_cdf
    IPython.get_ipython().user_ns.pop('selected_cdf_paths', None)
    selected_list.options = []
    for w in [browse_btn, add_btn, confirm_btn, existing_list, selected_list]: w.disabled = False
    with out_run: clear_output()
    _refresh_existing()
    status_lbl.value = 'Reiniciado.'

def _on_confirm(b):
    global selected_cdf_paths
    opts = list(selected_list.options)
    if not opts: status_lbl.value = '⚠ Añade archivos a Selected primero.'; return

    selected_cdf_paths = [Path(v) for _, v in opts]
    pipeline.list_cdf_files = lambda data_dir=None: selected_cdf_paths
    # Exponer en el namespace del kernel para que la celda siguiente lo lea
    IPython.get_ipython().user_ns['selected_cdf_paths'] = selected_cdf_paths

    for w in [browse_btn, add_btn, confirm_btn, existing_list, selected_list]: w.disabled = True
    status_lbl.value = f'⏳ {len(selected_cdf_paths)} archivo(s) — ejecutando pipeline…'

    with out_run:
        clear_output(wait=True)
        ok = _exec_cell('3958fffb')
        if ok:
            status_lbl.value = f'✔ Pipeline completado — {len(selected_cdf_paths)} muestra(s).'
        else:
            status_lbl.value = '❌ Error en el pipeline (ver output).'

def _exec_cell(cell_id):
    import json, os
    ip = IPython.get_ipython()
    ns = ip.user_ns

    nb_path = ns.get('__vsc_ipynb_file__', None)
    if not nb_path:
        hits = sorted(Path(os.getcwd()).glob('*.ipynb'))
        nb_path = str(hits[0]) if hits else None

    if not nb_path or not Path(nb_path).exists():
        print(f'⚠ No se localizó el notebook.')
        print(f'  Ejecuta manualmente la celda del pipeline (id: {cell_id}).')
        status_lbl.value = '⚠ Ejecuta la celda del pipeline manualmente.'
        return False

    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)

    for cell in nb['cells']:
        if cell.get('id') == cell_id and cell['cell_type'] == 'code':
            src = ''.join(cell['source']).strip()
            if not src:
                return True
            result = ip.run_cell(src)
            return not (result.error_before_exec or result.error_in_exec)

    print(f'⚠ Celda "{cell_id}" no encontrada en el notebook.')
    return False

browse_btn.on_click(_on_browse); add_btn.on_click(_on_add); remove_btn.on_click(_on_remove)
confirm_btn.on_click(_on_confirm); reset_btn.on_click(_on_reset)

display(HTML('<h3>📂 Seleccionar archivos .cdf para el análisis</h3>'))
display(HTML('<i>Usa "+ Buscar archivos .cdf" para seleccionar archivos desde cualquier carpeta. '
             'Los archivos de la carpeta por defecto aparecen en el panel izquierdo.</i>'))
display(browse_btn)
display(HBox([
    VBox([Label('📂 Carpeta por defecto (opcional)'), existing_list]),
    VBox([add_btn, remove_btn], layout={'justify_content': 'center', 'margin': '40px 12px'}),
    VBox([Label('✅ A analizar (Selected)'), selected_list]),
]))
display(HBox([confirm_btn, reset_btn]))
display(status_lbl)
display(HTML('<hr style="margin-top:12px"><b>Output del pipeline:</b>'))
display(out_run)
"""

lines = new_src.splitlines(keepends=True)

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

cell_found = False
for cell in nb["cells"]:
    if cell.get("id") == "b068f5f9":
        cell["source"] = lines
        cell["outputs"] = []
        cell["execution_count"] = None
        cell_found = True
        print(f"Celda b068f5f9 actualizada ({len(lines)} líneas).")
        break

if not cell_found:
    print("ERROR: celda b068f5f9 no encontrada.")
    exit(1)

with open(nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("Notebook guardado.")

# Verificar
import json as _json
with open(nb_path, "r", encoding="utf-8") as f:
    nb2 = _json.load(f)
print(f"JSON válido: {len(nb2['cells'])} celdas.")
