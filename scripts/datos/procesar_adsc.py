"""Convierte el volcado ADS-C de OpenSky en tablas CSV comprimidas. Solo usa la librería estándar.

Uso:
    python scripts/datos/procesar_adsc.py datos/opensky_adsc/adsc_decoded.txt --salida datos/procesado

Genera:
    adsc_posiciones.csv.gz   todos los mensajes con posición (una fila por mensaje)
    adsc_atlantico.csv.gz    subconjunto de la caja del Atlántico Norte: 25–75°N, 75°W–5°W
    adsc_resumen.json        recuento de mensajes, centros ATSU y anomalías de lectura
"""
import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from traj4d.adsc import COLUMNAS, Estadisticas, abrir, leer_mensajes  # noqa: E402

CAJA_NAT = dict(lat_min=25.0, lat_max=75.0, lon_min=-75.0, lon_max=-5.0)


def en_caja(f):
    return (CAJA_NAT["lat_min"] <= f["lat"] <= CAJA_NAT["lat_max"]
            and CAJA_NAT["lon_min"] <= f["lon"] <= CAJA_NAT["lon_max"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entrada", nargs="+")
    ap.add_argument("--salida", default="datos/procesado")
    a = ap.parse_args()
    out = Path(a.salida)
    out.mkdir(parents=True, exist_ok=True)

    est = Estadisticas()
    centros, centros_caja = Counter(), Counter()
    n_todo = n_caja = 0
    with gzip.open(out / "adsc_posiciones.csv.gz", "wt", newline="", compresslevel=6) as f1, \
         gzip.open(out / "adsc_atlantico.csv.gz", "wt", newline="", compresslevel=6) as f2:
        w1 = csv.DictWriter(f1, fieldnames=COLUMNAS)
        w2 = csv.DictWriter(f2, fieldnames=COLUMNAS)
        w1.writeheader()
        w2.writeheader()
        for ruta in a.entrada:
            with abrir(ruta) as fh:
                for fila in leer_mensajes(fh, est):
                    if fila["t"] is None or fila["lon"] is None:
                        est.sumar("posicion_incompleta")
                        continue
                    w1.writerow(fila)
                    n_todo += 1
                    centros[fila["centro"]] += 1
                    if en_caja(fila):
                        w2.writerow(fila)
                        n_caja += 1
                        centros_caja[fila["centro"]] += 1

    resumen = dict(estadisticas=dict(est), filas_posiciones=n_todo, filas_atlantico=n_caja,
                   caja_atlantico=CAJA_NAT, centros=dict(centros.most_common()),
                   centros_en_caja=dict(centros_caja.most_common()))
    (out / "adsc_resumen.json").write_text(json.dumps(resumen, indent=1, ensure_ascii=False))
    print(json.dumps({k: v for k, v in resumen.items() if k not in ("centros",)}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
