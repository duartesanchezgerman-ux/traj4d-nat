#!/bin/bash
# =====================================================================================================
# PASO 1 — Descargar (o reanudar) ERA5 de TODOS los días de estudio del artículo A1.
#
# Uso, desde la carpeta del proyecto en la Terminal del Mac:
#     bash scripts/1_descargar_era5.sh
#
# - Se puede lanzar cuantas veces haga falta: los días que ya están en datos/era5 se saltan.
# - Si ya hay otra descarga en marcha, NO lanza otra (duplicaría peticiones en la cola de Copernicus).
# - Mantiene el Mac despierto mientras descarga (caffeinate). La tapa puede ir cerrada solo si el Mac
#   está enchufado a la corriente.
# - La lista de días se guarda en docs/dias_estudio_era5.txt (va al repositorio: así el estudio es
#   reproducible). La primera vez se reconstruye a partir de tu historial de la Terminal (los días que
#   ya pediste con descargar_era5.py) más los 10 días de estudio iniciales.
# - Al terminar comprueba qué días faltan y guarda la huella SHA-256 de cada archivo en docs/ERA5_SHA256.txt.
# =====================================================================================================
set -u
cd "$(dirname "$0")/.." || exit 1
LISTA=docs/dias_estudio_era5.txt
LOG=resultados/descarga_era5.log
mkdir -p resultados datos/era5 docs

echo "== Paso 1: descarga de ERA5 ($(date '+%d/%m/%Y %H:%M')) =="

if pgrep -f "descargar_era5.py" >/dev/null 2>&1; then
  echo "Ya hay una descarga de ERA5 en marcha en otra terminal. No lanzo otra."
  echo "Espera a que termine (o ciérrala) y vuelve a ejecutar este script."
  exit 1
fi
if [ ! -f "$HOME/.cdsapirc" ]; then
  echo "Falta el archivo ~/.cdsapirc con tu clave de Copernicus (ver la cabecera de scripts/datos/descargar_era5.py)."
  exit 1
fi
if [ -f .venv/bin/activate ]; then source .venv/bin/activate; fi

# ---- 1. Lista de días de estudio ------------------------------------------------------------------
if [ ! -s "$LISTA" ]; then
  python3 - "$LISTA" <<'PY'
import re, sys
from datetime import date, timedelta
from pathlib import Path
salida = Path(sys.argv[1])
# Los 10 días de estudio iniciales (Fase 2) siempre forman parte del estudio.
iniciales = ["2023-07-08", "2023-07-09", "2023-07-15", "2023-07-16", "2023-11-01",
             "2023-11-03", "2023-11-04", "2023-11-05", "2023-11-09", "2024-12-25"]
DIA_DE_PRUEBA = "2024-01-15"          # solo se usó para probar la descarga; no es día de estudio
fecha = re.compile(r"\b(20\d\d-\d\d-\d\d)\b")
pedidos = set()
for h in (Path.home() / ".zsh_history", Path.home() / ".bash_history"):
    if not h.exists():
        continue
    for linea in h.read_bytes().decode("utf-8", "ignore").splitlines():
        if "descargar_era5.py" not in linea:
            continue
        if "--dias" in linea:
            pedidos.update(fecha.findall(linea.split("--dias", 1)[1].split("--area")[0]))
        if "--desde" in linea and "--hasta" in linea:
            a = fecha.search(linea.split("--desde", 1)[1]); b = fecha.search(linea.split("--hasta", 1)[1])
            if a and b:
                d, fin = date.fromisoformat(a.group(1)), date.fromisoformat(b.group(1))
                while d <= fin:
                    pedidos.add(d.isoformat()); d += timedelta(days=1)
dias = sorted((set(iniciales) | pedidos) - {DIA_DE_PRUEBA})
cab = ["# Días de estudio de ERA5 del artículo A1 (uno por línea). Generado por scripts/1_descargar_era5.sh",
       f"# a partir del historial de la Terminal + los 10 días iniciales. Total: {len(dias)} días."]
salida.write_text("\n".join(cab + dias) + "\n", encoding="utf-8")
print(f"Lista de días creada: {salida} ({len(dias)} días)")
if len(dias) != 70:
    print(f"ATENCIÓN: esperaba 70 días y salen {len(dias)}. Revisa {salida} antes de seguir "
          "(puedes editarla a mano: un día por línea) y avisa a Claude.")
PY
fi
DIAS=$(grep -v '^#' "$LISTA" | tr '\n' ' ')
N=$(echo $DIAS | wc -w | tr -d ' ')
echo "Días de estudio: $N (lista en $LISTA)"
if [ "$N" -eq 0 ]; then echo "La lista está vacía. Avisa a Claude."; exit 1; fi

# ---- 2. Descarga (una petición por mes, en serie, reanudable) --------------------------------------
echo "Descargando lo que falte. La cola de Copernicus puede tardar horas; no cierres esta ventana."
caffeinate -i python scripts/datos/descargar_era5.py --agrupar-por-mes --dias $DIAS --area 70 -60 30 0 2>&1 | tee -a "$LOG"

# ---- 3. Comprobación final ---------------------------------------------------------------------------
FALTAN=""
for d in $DIAS; do
  [ -f "datos/era5/era5_pl_$d.nc" ] || FALTAN="$FALTAN $d"
done
shasum -a 256 datos/era5/era5_pl_*.nc > docs/ERA5_SHA256.txt
echo "Huellas SHA-256 guardadas en docs/ERA5_SHA256.txt"
if [ -n "$FALTAN" ]; then
  echo "FALTAN $(echo $FALTAN | wc -w | tr -d ' ') días:$FALTAN"
  echo "Vuelve a lanzar este mismo script (reanuda lo que falte)."
  exit 1
fi
echo "OK: están los $N días. Siguiente paso:  bash scripts/2_revision_completa.sh"
