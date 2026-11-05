"""Extrae de la base de aeronaves de OpenSky (108 MB) solo los aviones que aparecen en el ADS-C.
Solo usa la librería estándar.

Uso:
    python scripts/datos/extraer_aeronaves.py datos/opensky_meta/aircraft-database-complete-2025-08.csv \
        datos/procesado/adsc_posiciones.csv.gz datos/procesado/aeronaves_adsc.csv
"""
import csv
import gzip
import sys

CAMPOS = ["icao24", "registration", "typecode", "model", "manufacturerName", "operator", "operatorIcao",
          "operatorCallsign", "owner", "built", "timestamp"]


def main(bd, adsc, salida):
    icaos, regs = set(), set()
    with gzip.open(adsc, "rt") as f:
        for r in csv.DictReader(f):
            if r["icao24"]:
                icaos.add(r["icao24"])
            if r["matricula"]:
                regs.add(r["matricula"].upper().replace("-", ""))
    n = 0
    # La base de OpenSky entrecomilla los campos con comillas simples.
    with open(bd, newline="", encoding="utf-8", errors="replace") as f, open(salida, "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=CAMPOS)
        w.writeheader()
        for r in csv.DictReader(f, quotechar="'"):
            ic = r.get("icao24", "").lower()
            rg = (r.get("registration") or "").upper().replace("-", "")
            if ic in icaos or (rg and rg in regs):
                w.writerow({k: r.get(k, "") for k in CAMPOS})
                n += 1
    print(f"icao24 en ADS-C: {len(icaos)} | matrículas: {len(regs)} | filas extraídas: {n}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
