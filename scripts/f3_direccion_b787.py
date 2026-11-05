"""Defecto de la dirección del viento en el ADS-C (sobre todo Boeing 787): ¿de qué depende?

Con 10 días el defecto afectaba al 4 % de los informes meteorológicos de la familia 787; con 70 días, al 13,5 %.
Este script lo desglosa para saber si es un rasgo de ciertos aviones, operadores, periodos o versiones:
  * diferencia de dirección con ERA5 (histograma) y cociente de velocidades;
  * fracción afectada por mes, por tipo (B788/B789/B78X), por operador y por avión (icao24);
  * si un avión afectado lo está siempre o solo a ratos (persistencia).
Mismo criterio que f2_validacion_era5.py: |Δdir| > 90°, viento observado > 10 m/s, cociente 0,6–1,5.

Uso:  PYTHONPATH=src python scripts/f3_direccion_b787.py | tee resultados/resultados_f3_direccion_b787.txt
Salida adicional: datos/procesado/direccion_viento_todos.csv (todos los informes, con la marca).
"""
import glob
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from traj4d import atmosfera as atm
from traj4d.viento import CampoERA5
from traj4d.vuelos import cargar_informes

warnings.filterwarnings("ignore")
KT = atm.KT

df = cargar_informes("datos/procesado/adsc_atlantico.csv.gz")
df = df[df.viento_kt.notna() & df.temp_c.notna() & df.alt_ft.notna()].copy()
wd = np.radians(df.viento_dir_deg)
df["u_obs"], df["v_obs"] = -df.viento_kt * KT * np.sin(wd), -df.viento_kt * KT * np.cos(wd)
df["fl"] = df.alt_ft / 100.0
df["dia"] = df.fecha.dt.strftime("%Y-%m-%d")
archivos = {Path(f).stem.replace("era5_pl_", ""): f for f in glob.glob("datos/era5/era5_pl_*.nc")}
if os.environ.get("DIAS_ESTUDIO"):
    lista = {l.split("#")[0].strip() for l in open(os.environ["DIAS_ESTUDIO"], encoding="utf-8")} - {""}
    archivos = {k: v for k, v in archivos.items() if k in lista}
partes = []
for dia, f in sorted(archivos.items()):
    sub = df[df.dia == dia]
    if sub.empty:
        continue
    u, v, _ = CampoERA5(f).muestra_lote(sub.t, sub.fl, sub.lat, sub.lon)
    partes.append(sub.assign(u_era=u, v_era=v)[np.isfinite(u)])
d = pd.concat(partes)
d["vel_obs"], d["vel_era"] = np.hypot(d.u_obs, d.v_obs), np.hypot(d.u_era, d.v_era)
d["ddir"] = ((np.degrees(np.arctan2(-d.u_obs, -d.v_obs)) - np.degrees(np.arctan2(-d.u_era, -d.v_era)) + 180) % 360) - 180
d["dir_obs"] = np.degrees(np.arctan2(-d.u_obs, -d.v_obs)) % 360
d["dir_era"] = np.degrees(np.arctan2(-d.u_era, -d.v_era)) % 360
d["corrupta"] = (d.ddir.abs() > 90) & (d.vel_obs > 10) & (d.vel_obs / d.vel_era).between(0.6, 1.5)
aer = pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str).drop_duplicates("icao24").set_index("icao24")
d["tipo"] = d.icao24.map(aer.typecode)
d["operador"] = d.icao24.map(aer.operatorIcao).fillna("?")
d["mes"] = d.fecha.dt.strftime("%Y-%m")
d[["t", "icao24", "tipo", "operador", "lat", "lon", "fl", "vel_obs", "vel_era", "dir_obs", "dir_era", "ddir", "corrupta"]].to_csv(
    "datos/procesado/direccion_viento_todos.csv", index=False)

pd.set_option("display.width", 200)
b = d[d.tipo.str[:3] == "B78"]
print("=" * 100)
print(f"Informes meteorológicos comparados: {len(d)} en {d.dia.nunique()} días | marcados: {int(d.corrupta.sum())} "
      f"({100 * d.corrupta.mean():.2f} %) | familia 787: {len(b)} informes, {100 * b.corrupta.mean():.1f} % marcados")
c = d[d.corrupta]
print("\n1. Diferencia de dirección de los marcados (|Δdir|, grados): "
      + " | ".join(f"p{q}: {np.percentile(c.ddir.abs(), q):.0f}" for q in (5, 25, 50, 75, 95)))
print("   cociente de velocidades obs/ERA5 de los marcados: mediana "
      f"{(c.vel_obs / c.vel_era).median():.3f} | p5–p95 {(c.vel_obs / c.vel_era).quantile(.05):.2f}–"
      f"{(c.vel_obs / c.vel_era).quantile(.95):.2f}")
print("   signo de Δdir en los marcados: " + str(np.sign(c.ddir).value_counts().to_dict()))
print("\n2. Familia 787 por tipo:")
print(b.groupby("tipo").corrupta.agg(["size", "mean"]).assign(mean=lambda x: (100 * x["mean"]).round(1)).to_string())
print("\n3. Familia 787 por mes (fracción marcada, %):")
print(b.groupby("mes").corrupta.agg(["size", "mean"]).assign(mean=lambda x: (100 * x["mean"]).round(1)).to_string())
print("\n4. Familia 787 por operador (≥ 100 informes):")
op = b.groupby("operador").corrupta.agg(["size", "mean"])
print(op[op["size"] >= 100].assign(mean=lambda x: (100 * x["mean"]).round(1)).sort_values("mean", ascending=False).to_string())
print("\n5. Familia 787 por avión: fracción de informes marcados de cada avión (aviones con ≥ 20 informes)")
av = b.groupby("icao24").corrupta.agg(["size", "mean"])
av = av[av["size"] >= 20]
print(f"   aviones: {len(av)} | con 0 %: {(av['mean'] == 0).sum()} | con > 0 y < 50 %: "
      f"{((av['mean'] > 0) & (av['mean'] < .5)).sum()} | con ≥ 50 %: {(av['mean'] >= .5).sum()} | "
      f"con ≥ 90 %: {(av['mean'] >= .9).sum()}")
print("   (si el defecto fuera de ciertos aviones, la fracción por avión sería casi siempre 0 % o casi 100 %)")
print("\n6. Persistencia por avión y día: ¿un avión afectado lo está todo el día?")
ad = b.groupby(["icao24", "dia"]).corrupta.agg(["size", "mean"])
ad = ad[ad["size"] >= 5]
print(f"   pares avión-día con ≥ 5 informes: {len(ad)} | 0 %: {(ad['mean'] == 0).sum()} | intermedio: "
      f"{((ad['mean'] > 0) & (ad['mean'] < 1)).sum()} | 100 %: {(ad['mean'] == 1).sum()}")
print("\n7. Familia 787: fracción marcada según la dirección real del viento (ERA5, de dónde sopla) y Δdir medio")
b2 = b[b.vel_obs > 10]
print(b2.groupby(pd.cut(b2.dir_era, np.arange(0, 361, 30)), observed=True).agg(
    n=("corrupta", "size"), marcada=("corrupta", "mean"),
    ddir_marcados=("ddir", lambda x: x[b2.loc[x.index, "corrupta"]].median())).round(3).to_string())
print("\n8. Fuera de la familia 787 (marcados por tipo):")
print(c[c.tipo.str[:3] != "B78"].tipo.value_counts().to_string())
