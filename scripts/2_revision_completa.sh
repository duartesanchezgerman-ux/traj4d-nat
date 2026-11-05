#!/bin/bash
# =====================================================================================================
# PASO 2 — Revisión completa con todos los días de ERA5, en el propio Mac (sin pasar los datos a Claude).
#
# Uso, cuando el paso 1 haya terminado sin días pendientes:
#     bash scripts/2_revision_completa.sh
#
# Qué hace, en orden (si algo falla, se para y dice dónde):
#   0. Comprueba que no hay una descarga en marcha y que están todos los días de docs/dias_estudio_era5.txt.
#   1. Comprueba cada archivo de ERA5: se abre, tiene 24 horas y cubre el dominio 30–70°N, 60°W–0°.
#   2. Guarda una copia de los resultados actuales (10 días) en resultados/archivo_10dias/ para comparar.
#   3. Pasa las pruebas automáticas. 3b. La primera vez, reprocesa el volcado ADS-C con la hora en UTC.
#   4. Rehace el análisis: validación de ERA5, comparación de rutas, estadística, estructura/elección y
#      figuras en inglés y en español.
#   5. Escribe resultados/RESUMEN_COMPARATIVO.txt (10 días frente a todos) y hace un commit con fecha y hora.
# Tiempo estimado en un Mac M4: 1–3 horas. Mantiene el Mac despierto mientras trabaja.
# Al terminar, dile a Claude "revisión hecha": leerá los resultados de la carpeta y actualizará el artículo.
# =====================================================================================================
set -u -o pipefail
cd "$(dirname "$0")/.." || exit 1
LISTA=docs/dias_estudio_era5.txt
R=resultados
export PYTHONPATH=src
export DIAS_ESTUDIO="$LISTA"            # los scripts de análisis solo usan los días de la lista
INICIO=$(date +%s)

paso() {  # paso "nombre" comando...   -> ejecuta, guarda la salida y para si falla
  local nombre="$1"; shift
  echo ""; echo "== $nombre  ($(date '+%H:%M'))"
  if ! "$@"; then
    echo "FALLO en: $nombre. No sigo. Copia este mensaje a Claude."
    exit 1
  fi
}

echo "== Paso 2: revisión completa ($(date '+%d/%m/%Y %H:%M')) =="
if pgrep -f "descargar_era5.py" >/dev/null 2>&1; then
  echo "Hay una descarga de ERA5 en marcha. Espera a que termine."; exit 1
fi
[ -s "$LISTA" ] || { echo "No existe $LISTA. Ejecuta antes: bash scripts/1_descargar_era5.sh"; exit 1; }
if [ -f .venv/bin/activate ]; then source .venv/bin/activate; fi
caffeinate -i -w $$ &            # el Mac no se duerme mientras dure este script

# ---- Dependencias (instala lo que falte, p. ej. matplotlib para las figuras) ------------------------
paso "Dependencias de Python" python -m pip install -q -r requirements.txt

# ---- 0 y 1. Datos completos y sanos --------------------------------------------------------------------
paso "Comprobación de los archivos de ERA5" python3 - "$LISTA" <<'PY'
import sys
from pathlib import Path
import numpy as np
import xarray as xr
dias = [l.split("#")[0].strip() for l in open(sys.argv[1], encoding="utf-8")]
dias = [d for d in dias if d]
malos, faltan = [], []
for d in dias:
    f = Path(f"datos/era5/era5_pl_{d}.nc")
    if not f.exists():
        faltan.append(d); continue
    try:
        with xr.open_dataset(f) as ds:
            dim = "valid_time" if "valid_time" in ds.dims else "time"
            n = ds.sizes[dim]
            lat, lon = ds["latitude"].values, ds["longitude"].values
            ok = (n == 24 and lat.min() <= 30.0 and lat.max() >= 70.0 and lon.min() <= -60.0 and lon.max() >= 0.0
                  and all(v in ds for v in ("u", "v", "t")) and np.isfinite(ds["u"].values).all())
            if not ok:
                malos.append(f"{d} (horas {n}, lat {lat.min()}–{lat.max()}, lon {lon.min()}–{lon.max()})")
    except Exception as e:  # noqa: BLE001
        malos.append(f"{d} (no se abre: {e})")
print(f"Días en la lista: {len(dias)} | presentes y correctos: {len(dias) - len(faltan) - len(malos)}")
if faltan:
    print("FALTAN:", " ".join(faltan), "\n-> ejecuta otra vez: bash scripts/1_descargar_era5.sh")
if malos:
    print("ARCHIVOS CON PROBLEMAS:\n  " + "\n  ".join(malos))
sys.exit(1 if (faltan or malos) else 0)
PY

# ---- 2. Copia de los resultados con 10 días (solo la primera vez) -------------------------------------
if [ ! -d "$R/archivo_10dias" ]; then
  mkdir -p "$R/archivo_10dias" datos/procesado/archivo_10dias
  cp "$R"/resultados_f2_validacion_era5.txt "$R"/resultados_f2_comparacion.txt "$R"/resultados_f3_estadistica.txt \
     "$R"/resultados_f3_estructura.txt "$R/archivo_10dias/" 2>/dev/null
  cp datos/procesado/comparacion_tramos.csv datos/procesado/estructura_tramos.csv datos/procesado/validacion_era5_adsc.csv \
     datos/procesado/archivo_10dias/ 2>/dev/null
  echo "Resultados de 10 días guardados en $R/archivo_10dias/"
fi

# ---- 3. Pruebas ---------------------------------------------------------------------------------------
paso "Pruebas automáticas" python -m pytest -q tests

# ---- 3b. Hora del volcado ADS-C: de hora local de Europa central a UTC (corrección del 3-11-2026) ------
# El Timestamp del volcado de OpenSky viene en CET/CEST, no en UTC. Se reprocesa una sola vez y se comprueba
# que las tablas salen idénticas (sha256) a las que usó Claude con los 10 días.
if [ ! -f datos/procesado/HORA_UTC ]; then
  mkdir -p datos/procesado/_hora_local_antiguo
  cp datos/procesado/*.csv datos/procesado/*.csv.gz datos/procesado/_hora_local_antiguo/ 2>/dev/null
  paso "Reprocesado del volcado ADS-C con hora UTC" python scripts/datos/procesar_adsc.py \
       datos/opensky_adsc/adsc_decoded.txt --salida datos/procesado
  paso "Comprobación sha256 de las tablas ADS-C" python3 - <<'PY'
import gzip, hashlib, sys
ok = True
for linea in open("docs/ADSC_UTC_SHA256.txt", encoding="utf-8"):
    h, nombre = linea.split()[0], linea.split()[1]
    with gzip.open(f"datos/procesado/{nombre}.gz", "rb") as f:
        hh = hashlib.sha256(f.read()).hexdigest()
    print(nombre, "OK" if hh == h else f"DISTINTO ({hh})"); ok &= hh == h
sys.exit(0 if ok else 1)
PY
  date > datos/procesado/HORA_UTC
fi
paso "Cobertura y segmentos" bash -c "python scripts/f2_cobertura.py | tee $R/resultados_f2_cobertura.txt"

# ---- 4. Análisis completo -------------------------------------------------------------------------------
paso "Validación de ERA5 frente al viento medido" bash -c "python scripts/f2_validacion_era5.py | tee $R/resultados_f2_validacion_era5.txt"
paso "Comparación ruta volada / óptima (la parte larga)" bash -c "python scripts/f2_comparacion.py | tee $R/resultados_f2_comparacion.txt"
paso "Estadística (bootstrap por días)" bash -c "python scripts/f3_estadistica.py | tee $R/resultados_f3_estadistica.txt"
paso "Estructura frente a elección de puntos de paso" bash -c "python scripts/f3_estructura.py | tee $R/resultados_f3_estructura.txt"
paso "Figuras en inglés" python scripts/f3_figuras.py --idioma en
paso "Figuras en español" python scripts/f3_figuras.py --idioma es

# ---- 5. Resumen comparativo y commit -------------------------------------------------------------------
python3 - <<'PY' | tee "$R/RESUMEN_COMPARATIVO.txt"
from pathlib import Path
claves = ("Informes comparados", "error vectorial RMS", "error vectorial mediano", "sesgo de temperatura",
          "sesgo del módulo con viento", "Tramos candidatos", "exceso de combustible, media",
          "exceso de combustible, mediana", "fracción de tramos con exceso", "ganancia del óptimo",
          "exceso medio en kg", "exceso medio B789", "exceso medio A359", "estructural (restringido",
          "de elección (volada", "volada / ruta planificada", "El punto planificado coincide",
          "fracción estructural", "0,70·MTOW")
archivos = ("resultados_f2_validacion_era5.txt", "resultados_f2_comparacion.txt",
            "resultados_f3_estadistica.txt", "resultados_f3_estructura.txt")
print("RESUMEN COMPARATIVO: 10 días (archivo) frente a todos los días de estudio (nuevo)\n")
for a in archivos:
    viejo, nuevo = Path("resultados/archivo_10dias") / a, Path("resultados") / a
    print("=" * 100); print(a)
    for etiqueta, f in (("10 días", viejo), ("TODOS  ", nuevo)):
        if not f.exists():
            continue
        for l in f.read_text(encoding="utf-8").splitlines():
            if any(c in l for c in claves):
                print(f"  [{etiqueta}] {l.strip()}")
PY

MIN=$(( ($(date +%s) - INICIO) / 60 ))
git add docs/dias_estudio_era5.txt docs/ERA5_SHA256.txt "$R" articulo/en/figures articulo/es/figures 2>/dev/null
git -c user.name="Germán Duarte Sánchez" -c user.email="duartesanchezgerman@gmail.com" commit -q -m "Resultados con $(grep -vc '^#' "$LISTA") días de ERA5 ($(date '+%d/%m/%Y %H:%M'); revisión completa en $MIN min)" \
  && echo "Commit hecho: $(git log --oneline -1)"
echo ""
echo "LISTO en $MIN minutos. Dile a Claude: \"revisión hecha\". El resumen está en $R/RESUMEN_COMPARATIVO.txt"
