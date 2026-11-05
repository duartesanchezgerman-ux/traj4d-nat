"""Comprobación de la hora del volcado ADS-C con trazas ADS-B (adsb.lol, hora UTC).

El campo Timestamp del volcado de OpenSky viene en hora local de Europa central, no en UTC; traj4d.adsc._ts
lo convierte a UTC. Este script lo comprueba: para cada informe ADS-C (ya en UTC) busca, en las trazas ADS-B
del mismo avión (archivos-día D−1, D y D+1), el punto ADS-B más cercano en el espacio. Si está a menos de
8 km, la diferencia de tiempos t_ADS-B − t_ADS-C debe ser ≈ 0 s. También da la diferencia que habría con la
hora local tomada como UTC (−3600 s en invierno, −7200 s en verano) y cuenta los informes corregidos por la
vuelta de hora (traj4d.vuelos.corregir_vuelta_de_hora).

Necesita datos/adsblol/<día>/<icao24>.json.gz (scripts/datos/adsblol_extraer.py); usa los días completos
que haya (los que tienen el archivo _COMPLETO).

Uso:  PYTHONPATH=src python scripts/f2_hora_adsc.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from traj4d.cruces import leer_traza
from traj4d.ruta import distancia
from traj4d.vuelos import corregir_vuelta_de_hora

BASE = Path("datos/adsblol")
DIST_MAX_KM = 8.0

pos = pd.read_csv("datos/procesado/adsc_posiciones.csv.gz", usecols=["t", "icao24", "matricula", "lat", "lon"],
                  dtype={"icao24": "string", "matricula": "string"})
pos = corregir_vuelta_de_hora(pos)
n_corr, n_tot = int(pos.t_corregida_1h.sum()), len(pos)
caja = pos.lat.between(25, 75) & pos.lon.between(-75, -5)
n_corr_caja = int((pos.t_corregida_1h & caja).sum())
minuto = ((pos.t[pos.t_corregida_1h] + 3600) // 60 % 60).astype(int)

pos = pos[pos.icao24.notna()].copy()
pos["dia"] = pd.to_datetime(pos.t, unit="s").dt.strftime("%Y-%m-%d")
dias = sorted(p.name for p in BASE.iterdir() if p.name[:1] == "2" and (p / "_COMPLETO").exists())
LISTA = Path("docs/dias_hora_adsc.txt")       # días fijados para el artículo (si existe la lista)
if LISTA.exists():
    fijos = [l.split("#")[0].strip() for l in LISTA.read_text(encoding="utf-8").splitlines()]
    fijos = [x for x in fijos if x]
    faltan = sorted(set(fijos) - set(dias))
    if faltan:
        raise SystemExit(f"Faltan trazas ADS-B completas de {len(faltan)} días de {LISTA}: {' '.join(faltan)}")
    dias = fijos
cerca = lambda d, k: (pd.Timestamp(d) + pd.Timedelta(days=k)).strftime("%Y-%m-%d")

filas = []
for dia in dias:
    for hx, g in pos[pos.dia == dia].groupby("icao24"):
        fs = [BASE / cerca(dia, k) / f"{hx}.json.gz" for k in (-1, 0, 1)]
        fs = [f for f in fs if f.exists()]
        if not fs:
            continue
        b = pd.concat([leer_traza(f) for f in fs]).drop_duplicates("t").sort_values("t")
        b = b[b.alt_ft.notna()]
        for r in g.itertuples():
            w = b[(b.t > r.t - 4 * 3600) & (b.t < r.t + 4 * 3600)]
            if len(w) < 2:
                continue
            d = distancia(r.lat, r.lon, w.lat.to_numpy(), w.lon.to_numpy()) / 1e3
            i = int(np.argmin(d))
            if d[i] <= DIST_MAX_KM:
                filas.append(dict(dia=dia, icao24=hx, t_adsc=r.t, t_adsb=float(w.t.iat[i]), dist_km=float(d[i])))

df = pd.DataFrame(filas)
local = pd.to_datetime(df.t_adsc, unit="s", utc=True).dt.tz_convert("Europe/Zurich")
df["desfase_utc_s"] = df.t_adsb - df.t_adsc
df["desfase_si_local_s"] = df.desfase_utc_s - local.map(lambda x: x.utcoffset().total_seconds())
df["mes"] = df.dia.str[:7]
Path("datos/procesado").mkdir(exist_ok=True)
df.to_csv("datos/procesado/hora_adsc_vs_adsb.csv", index=False)

print("=" * 100)
print("HORA DEL VOLCADO ADS-C FRENTE A ADS-B (adsb.lol, UTC)")
print(f"Días ADS-B completos: {len(dias)} | coincidencias a < {DIST_MAX_KM:.0f} km: {len(df)} informes, "
      f"{df.icao24.nunique()} aviones, {df.dia.nunique()} días ({', '.join(sorted(df.mes.unique()))})")
ok = df.desfase_utc_s.abs() < 600
q = df.desfase_utc_s[ok].quantile([.05, .5, .95])
print(f"t_ADS-B − t_ADS-C con la hora convertida a UTC: mediana {q[.5]:+.1f} s | p5 {q[.05]:+.1f} s | "
      f"p95 {q[.95]:+.1f} s | fuera de ±600 s: {int((~ok).sum())}")
print("Con la hora local tomada como UTC (por mes, mediana en s):")
print(df.groupby("mes").desfase_si_local_s.agg(["size", "median"]).round(0).to_string())
print(f"\nVuelta de hora: {n_corr} informes corregidos de {n_tot} ({100 * n_corr / n_tot:.2f} %); "
      f"en la caja 25–75°N, 75–5°W: {n_corr_caja}")
print("  minuto de generación de los corregidos:", minuto.value_counts().sort_index().to_dict())
