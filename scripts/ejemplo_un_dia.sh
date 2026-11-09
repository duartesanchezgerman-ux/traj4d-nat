#!/bin/bash
# Quick reproducibility check on ONE study day, without downloading the whole ERA5 set.
#
#   bash scripts/ejemplo_un_dia.sh               # default day: 2023-07-09 (65 oceanic segments)
#   bash scripts/ejemplo_un_dia.sh 2023-10-01    # any of the 70 days in docs/dias_estudio_era5.txt
#
# Needs: the processed ADS-C tables (see "Data" in README.md; ~1 min) and ~/.cdsapirc for one day of
# ERA5 (~85 MB). Runs the ERA5 validation and the flown-vs-optimal comparison for that day, compares the
# results with the full 70-day run and draws a one-day version of Figures 2 and 5. A few minutes in total.
# Results go to resultados/ejemplo_<day>/; the files of the full run are not touched.
set -eu
cd "$(dirname "$0")/.." || exit 1
DIA=${1:-2023-07-09}
OUT=resultados/ejemplo_$DIA
ERA5=datos/era5/era5_pl_$DIA.nc

for f in datos/procesado/adsc_atlantico.csv.gz datos/procesado/aeronaves_adsc.csv; do
  [ -f "$f" ] || { echo "Missing $f: run the two commands of the 'Data' section of README.md first."; exit 1; }
done
grep -q "^$DIA" docs/dias_estudio_era5.txt || { echo "$DIA is not a study day (docs/dias_estudio_era5.txt)."; exit 1; }

mkdir -p datos/era5 "$OUT"
if [ ! -f "$ERA5" ]; then
  echo "== Downloading ERA5 for $DIA (one file, ~85 MB; the Copernicus queue may take a few minutes)"
  python scripts/datos/descargar_era5.py --dia "$DIA" --area 70 -60 30 0
fi
REF=$(grep "$DIA" docs/ERA5_SHA256.txt | cut -d' ' -f1)
SHA=$( (sha256sum "$ERA5" 2>/dev/null || shasum -a 256 "$ERA5") | cut -d' ' -f1)
if [ "$SHA" = "$REF" ]; then echo "ERA5 file identical to the one used in the article (SHA-256)."
else echo "Note: the ERA5 file differs byte-wise from the article's (Copernicus may re-encode files); the numbers below are the real test."; fi

echo "$DIA" > "$OUT/dias.txt"
echo "== ERA5 vs on-board wind"
DIAS_ESTUDIO="$OUT/dias.txt" SALIDA="$OUT" PYTHONPATH=src python scripts/f2_validacion_era5.py > "$OUT/resultados_validacion_era5.txt"
echo "== Flown vs wind-optimal paths"
DIAS_ESTUDIO="$OUT/dias.txt" SALIDA="$OUT" PYTHONPATH=src python scripts/f2_comparacion.py > "$OUT/resultados_comparacion.txt"
PYTHONPATH=src python scripts/ejemplo_un_dia.py "$DIA" "$OUT"
