"""Fase 2: validación de ERA5 frente al viento y la temperatura MEDIDOS por los aviones (ADS-C, Tag 16).
Ejecutar:  PYTHONPATH=src python scripts/f2_validacion_era5.py | tee resultados/resultados_f2_validacion_era5.txt

Para cada informe con meteo de los días de estudio dentro del dominio de ERA5 descargado:
interpolación de ERA5 en (t, nivel de presión del FL, lat, lon) y comparación de u, v y T.
La altitud del ADS-C es altitud presión (referida a 1013,25 hPa): equivale directamente a nuestro FL.

Defecto de origen detectado (2 nov 2026): un ~1,4 % de los informes (el 4 % de los del Boeing 787, casi
ningún otro tipo) tiene la DIRECCIÓN del viento desplazada entre 110° y 130°, con el MÓDULO correcto.
Se clasifican como "dirección corrupta" si |Δdir| > 90°, el viento observado es > 10 m/s y el cociente
de módulos está entre 0,6 y 1,5 (hasta el 5-11-2026, 0,7–1,3: con 70 días se vio que dejaba pasar 29 informes
    girados del 787 con viento fuerte, en los que ERA5 infravalora algo el módulo). Se informa con y sin ellos. Nota: el criterio usa ERA5 como
referencia, pero ningún reanálisis confunde un viento de 15–40 m/s girándolo 120° con el módulo
exacto: el patrón es el de un fallo de codificación, no un error meteorológico.

Advertencia de independencia: ERA5 asimila observaciones de aviones (AMDAR y parte de los
informes ADS-C de algunos operadores). Por eso el error medido aquí es una cota INFERIOR del error
de una previsión, que es lo que realmente tendría un piloto en vuelo.
"""
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from traj4d import atmosfera as atm
from traj4d.vuelos import cargar_informes
from traj4d.viento import CampoERA5

def filtrar_dias_estudio(dias):
    """Si la variable de entorno DIAS_ESTUDIO apunta a una lista de días (uno por línea, '#' = comentario),
    solo se usan esos días. Sin la variable se usan todos los archivos de datos/era5 (comportamiento anterior)."""
    import os
    ruta = os.environ.get("DIAS_ESTUDIO")
    if not ruta:
        return dias
    lista = {l.split("#")[0].strip() for l in open(ruta, encoding="utf-8")} - {""}
    return {k: v for k, v in dias.items() if k in lista}


warnings.filterwarnings("ignore")
KT = atm.KT

df = cargar_informes("datos/procesado/adsc_atlantico.csv.gz")
df = df[df.viento_kt.notna() & df.temp_c.notna() & df.alt_ft.notna()].copy()
wd = np.radians(df.viento_dir_deg)
df["u_obs"] = -df.viento_kt * KT * np.sin(wd)
df["v_obs"] = -df.viento_kt * KT * np.cos(wd)
df["T_obs"] = df.temp_c + 273.15
df["fl"] = df.alt_ft / 100.0
df["dia"] = df.fecha.dt.strftime("%Y-%m-%d")

partes = []
_archivos = filtrar_dias_estudio({Path(f).stem.replace("era5_pl_", ""): f for f in glob.glob("datos/era5/era5_pl_*.nc")})
for f in sorted(_archivos.values()):
    dia = Path(f).stem.replace("era5_pl_", "")
    sub = df[df.dia == dia]
    if sub.empty:
        continue
    c = CampoERA5(f)
    u, v, dT = c.muestra_lote(sub.t, sub.fl, sub.lat, sub.lon)
    sub = sub.assign(u_era=u, v_era=v, T_era=dT + atm.temperatura(atm.fl_a_m(sub.fl)))
    partes.append(sub[np.isfinite(u)])
d = pd.concat(partes)
d["du"], d["dv"], d["dT"] = d.u_era - d.u_obs, d.v_era - d.v_obs, d.T_era - d.T_obs
d["err_vec"] = np.hypot(d.du, d.dv)
d["vel_obs"], d["vel_era"] = np.hypot(d.u_obs, d.v_obs), np.hypot(d.u_era, d.v_era)
ddir = ((np.degrees(np.arctan2(-d.u_obs, -d.v_obs)) - np.degrees(np.arctan2(-d.u_era, -d.v_era)) + 180) % 360) - 180
d["dir_corrupta"] = (np.abs(ddir) > 90) & (d.vel_obs > 10) & (d.vel_obs / d.vel_era).between(0.6, 1.5)
aer = pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str).drop_duplicates("icao24").set_index("icao24")
d["tipo"] = d.icao24.map(aer.typecode)
d["estacion"] = np.where(d.fecha.dt.month.isin([11, 12, 1, 2]), "invierno", "verano")


def tabla(g):
    return pd.Series({
        "n": len(g),
        "sesgo_u": g.du.mean(), "sesgo_v": g.dv.mean(),
        "rmse_u": np.sqrt((g.du**2).mean()), "rmse_v": np.sqrt((g.dv**2).mean()),
        "error_vectorial_rms": np.sqrt((g.err_vec**2).mean()), "error_vectorial_p50": g.err_vec.median(),
        "error_vectorial_p95": g.err_vec.quantile(.95),
        "sesgo_modulo": (g.vel_era - g.vel_obs).mean(), "viento_medio_obs": g.vel_obs.mean(),
        "sesgo_T": g.dT.mean(), "rmse_T": np.sqrt((g.dT**2).mean()),
    })


print("=" * 90)
print("ERA5 frente a viento y temperatura medidos por aviones (ADS-C Tag 16). Unidades: m/s y K")
print(f"Informes comparados: {len(d)} | días: {d.dia.nunique()} | FL {d.fl.min():.0f}–{d.fl.max():.0f} | "
      f"lat {d.lat.min():.1f}–{d.lat.max():.1f} | lon {d.lon.min():.1f}–{d.lon.max():.1f}")
pd.set_option("display.width", 200)
fam = d.tipo.str[:3]
print(f"\nDirección corrupta: {int(d.dir_corrupta.sum())} informes ({100 * d.dir_corrupta.mean():.2f} %) | "
      f"por familia: {fam[d.dir_corrupta].value_counts().to_dict()} | "
      f"fracción en la familia B78: {100 * d.dir_corrupta[fam == 'B78'].mean():.1f} %")
print("\nGLOBAL (bruto | sin dirección corrupta)")
print(pd.concat([tabla(d), tabla(d[~d.dir_corrupta])], axis=1, keys=["bruto", "filtrado"]).round(2).to_string())
d = d[~d.dir_corrupta]
print("\n--- En adelante, sin los informes con dirección corrupta ---")
print("\nPOR ESTACIÓN"); print(d.groupby("estacion").apply(tabla).round(2).to_string())
print("\nPOR DÍA"); print(d.groupby("dia").apply(tabla)[["n", "error_vectorial_rms", "sesgo_modulo", "viento_medio_obs", "sesgo_T"]].round(2).to_string())
d["clase_viento"] = pd.cut(d.vel_obs, [0, 20, 40, 60, 120], labels=["<20", "20–40", "40–60", ">60"])
print("\nPOR INTENSIDAD DEL VIENTO OBSERVADO [m/s]")
print(d.groupby("clase_viento").apply(tabla)[["n", "error_vectorial_rms", "error_vectorial_p95", "sesgo_modulo"]].round(2).to_string())
d[["t", "icao24", "tipo", "lat", "lon", "fl", "u_obs", "v_obs", "T_obs", "u_era", "v_era", "T_era"]].to_csv(
    "datos/procesado/validacion_era5_adsc.csv", index=False)
