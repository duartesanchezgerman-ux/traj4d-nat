"""Lanzador para VS Code: abre este archivo y pulsa el botón ▶ (Run Python File).
Ejecuta scripts/1_descargar_era5.sh, que es el que hace el trabajo (una sola versión de la lógica).
También se puede lanzar desde la terminal de VS Code:  bash scripts/1_descargar_era5.sh
"""
import subprocess
import sys
from pathlib import Path

script = Path(__file__).with_name("1_descargar_era5.sh")
sys.exit(subprocess.run(["bash", str(script)], cwd=script.parent.parent).returncode)
