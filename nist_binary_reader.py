"""
Lector de base de datos binaria NIST MAINLIB (.db, .in6, .sdb, .dsk files).

NIST MAINLIB es una base de datos binaria con estructura:
- MAINLIB.db: contiene datos espectrales comprimidos
- MAINLIB.in6: índice de 6 bytes por entrada (acceso rápido)
- MAINLIB.sdb: sinónimos/datos extendidos

Este módulo intenta leer directamente sin necesidad de exportar a .msp.

Formato general observado en NIST MAINLIB:
- Registro de índice de 6 bytes: [position(4 bytes), length(2 bytes)]
- Cada espectro: encabezado + picos

Notas:
- El formato exacto es parcialmente propietario
- Este es un aproximación basada en reverse-engineering
- Mejor solución: exportar desde NIST MS Search cuando sea posible
"""

from __future__ import annotations
from pathlib import Path
import struct
import logging
from typing import Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)


def read_nist_mainlib_index(mainlib_path: Path) -> list[Tuple[int, int]]:
    """
    Lee archivo .in6 (índice de NIST MAINLIB).
    
    Cada entrada ocupa 6 bytes:
    - 4 bytes: posición en el archivo .db
    - 2 bytes: longitud del registro
    """
    index_file = mainlib_path / "MAINLIB.in6"
    
    if not index_file.exists():
        logger.warning(f"Archivo índice no encontrado: {index_file}")
        return []
    
    entries = []
    try:
        with open(index_file, "rb") as f:
            while True:
                chunk = f.read(6)
                if len(chunk) < 6:
                    break
                
                # Formato: little-endian
                pos, length = struct.unpack("<IH", chunk)
                
                # Validar valores razonables
                if pos > 0 and length > 0 and length < 100000:
                    entries.append((pos, length))
    
    except Exception as e:
        logger.error(f"Error leyendo índice: {e}")
        return []
    
    logger.info(f"índice NIST leído: {len(entries)} entradas")
    return entries


def _parse_nist_record(data: bytes) -> dict:
    """
    Intenta parsear un registro binario NIST.
    
    Estructura aproximada (reverse-engineered):
    - Header variable
    - Campo ID (NIST#)
    - Nombre del compuesto
    - Fórmula, MW, etc.
    - Número de picos
    - Lista de picos (m/z, intensidad)
    
    NOTA: Este parser es parcial y puede requerir ajustes
    """
    spectrum = {
        "name": "Unknown",
        "formula": "",
        "mw": 0,
        "mz": np.array([]),
        "intensities": np.array([]),
        "nist_id": None
    }
    
    if len(data) < 10:
        return spectrum
    
    try:
        # Intentar extraer texto legible como nombres
        text_attempts = []
        for encoding in ["utf-8", "latin-1", "ascii"]:
            try:
                text = data.decode(encoding, errors="ignore")
                text_attempts.append(text)
            except:
                pass
        
        # Buscar campos conocidos
        for text in text_attempts:
            lines = text.split("\n")
            for line in lines:
                line = line.strip()
                
                if ":" in line:
                    key, _, value = line.partition(":")
                    key = key.strip().upper()
                    value = value.strip()
                    
                    if key == "NAME":
                        spectrum["name"] = value
                    elif key in ("MW", "MOLWT"):
                        try:
                            spectrum["mw"] = float(value)
                        except:
                            pass
                    elif key == "FORMULA":
                        spectrum["formula"] = value
                    elif key == "NIST#" or key == "ID":
                        spectrum["nist_id"] = value
        
        # Si no se encontró nombre, usar primeras palabras legibles
        if spectrum["name"] == "Unknown":
            for text in text_attempts:
                words = text.split()
                valid_words = [w for w in words if len(w) > 2 and w.isalnum()]
                if valid_words:
                    spectrum["name"] = " ".join(valid_words[:3])
                    break
    
    except Exception as e:
        logger.debug(f"Error parseando registro: {e}")
    
    return spectrum


def read_nist_mainlib_records(
    mainlib_path: Path,
    max_records: Optional[int] = None
) -> list[dict]:
    """
    Lee espectros directamente desde MAINLIB.db usando índice.
    
    Args:
        mainlib_path: directorio que contiene MAINLIB.db + MAINLIB.in6
        max_records: máximo de registros a leer (None = todos)
    
    Returns:
        Lista de espectros parsed
    """
    db_file = mainlib_path / "MAINLIB.db"
    
    if not db_file.exists():
        logger.error(f"Base de datos no encontrada: {db_file}")
        return []
    
    # Leer índice
    index = read_nist_mainlib_index(mainlib_path)
    
    if not index:
        logger.error("No se pudo leer índice NIST")
        return []
    
    spectra = []
    
    try:
        with open(db_file, "rb") as f:
            for record_num, (pos, length) in enumerate(index):
                if max_records and record_num >= max_records:
                    break
                
                # Leer registro desde posición
                f.seek(pos)
                data = f.read(length)
                
                if len(data) > 0:
                    spectrum = _parse_nist_record(data)
                    spectra.append(spectrum)
                
                # Log progreso cada 1000 registros
                if (record_num + 1) % 1000 == 0:
                    logger.info(f"  Procesados {record_num + 1} espectros...")
    
    except Exception as e:
        logger.error(f"Error leyendo MAINLIB.db: {e}")
    
    logger.info(f"Leyendo MAINLIB completado: {len(spectra)} espectros")
    return spectra


def try_load_nist_mainlib_binary(
    nist_dir: Path,
    max_attempts: int = 5000
) -> list[dict]:
    """
    Intenta cargar MAINLIB desde forma binaria sin exportar.
    
    Este es un mejor esfuerzo. Si falla, caerá a .msp tradicional.
    """
    logger.info("Intentando leer MAINLIB binario directamente...")
    
    mainlib_dir = nist_dir / "MAINLIB"
    
    if not mainlib_dir.exists():
        logger.warning(f"Directorio MAINLIB no encontrado: {mainlib_dir}")
        return []
    
    # Verificar que tenemos archivos necesarios
    db_file = mainlib_dir / "MAINLIB.db"
    in6_file = mainlib_dir / "MAINLIB.in6"
    
    if not (db_file.exists() and in6_file.exists()):
        logger.warning(f"Archivos MAINLIB.db o MAINLIB.in6 no encontrados")
        return []
    
    file_size_mb = db_file.stat().st_size / (1024 * 1024)
    logger.info(f"📊 MAINLIB.db: {file_size_mb:.1f} MB")
    
    # Leer con límite para no sobrecargar (primeros 5000 registros)
    spectra = read_nist_mainlib_records(mainlib_dir, max_records=max_attempts)
    
    if spectra:
        logger.info(f"✅ Se cargaron {len(spectra)} espectros de MAINLIB binario")
    
    return spectra
