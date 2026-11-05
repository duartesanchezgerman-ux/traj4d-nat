"""Fase 3: intervalos de confianza con bootstrap por DÍAS (remuestreo por bloques).
Ejecutar:  PYTHONPATH=src python scripts/f3_estadistica.py | tee resultados/resultados_f3_estadistica.txt

Por qué por días: los tramos y los informes de un mismo día comparten la situación meteorológica (el
mismo chorro y el mismo error de ERA5), así que no son independientes. Remuestrear tramos individuales
estrecharía artificialmente los intervalos. Se remuestrean días completos con reemplazo (B = 5000) y se
dan intervalos de percentil al 95 %.
Entradas: datos/procesado/comparacion_tramos.csv y datos/procesado/validacion_era5_adsc.csv
(generadas por f2_comparacion.py y f2_validacion_era5.py; se actualizan solas al añadir días).
"""
import numpy as np
import pandas as pd

B, SEMILLA = 5000, 20261102
rng = np.random.default_rng(SEMILLA)


def bootstrap_por_dia(df, col_dia, estadistico):
    dias = df[col_dia].unique()
    grupos = {d: g for d, g in df.groupby(col_dia)}
    est = np.empty(B)
    for b in range(B):
        muestra = pd.concat([grupos[d] for d in rng.choice(dias, size=len(dias), replace=True)])
        est[b] = estadistico(muestra)
    validas = np.isfinite(est)
    if not validas.all():
        print(f"    (aviso: {int((~validas).sum())} de {B} remuestras sin datos para este estadístico; se excluyen)")
    return estadistico(df), np.nanpercentile(est, 2.5), np.nanpercentile(est, 97.5), len(dias)


def linea(nombre, r, unidad):
    v, lo, hi, nd = r
    print(f"  {nombre:52s} {v:8.3f} {unidad}  IC95 [{lo:7.3f}, {hi:7.3f}]  ({nd} días)")


t = pd.read_csv("datos/procesado/comparacion_tramos.csv")
print("=" * 100)
print(f"RUTA VOLADA FRENTE A ÓPTIMO LATERAL — {len(t)} tramos")
linea("exceso de combustible, media", bootstrap_por_dia(t, "dia", lambda d: d.exceso_vs_opt_pct.mean()), "%")
linea("exceso de combustible, mediana", bootstrap_por_dia(t, "dia", lambda d: d.exceso_vs_opt_pct.median()), "%")
linea("fracción de tramos con exceso > 0,5 %", bootstrap_por_dia(t, "dia", lambda d: 100 * (d.exceso_vs_opt_pct > .5).mean()), "%")
linea("ganancia del óptimo frente a la ortodrómica, media", bootstrap_por_dia(t, "dia", lambda d: d.ganancia_opt_vs_gc_pct.mean()), "%")
linea("exceso medio en kg por tramo", bootstrap_por_dia(t, "dia", lambda d: (d.comb_volada - d.comb_opt).mean()), "kg")
for tipo in ("B789", "A359"):
    s = t[t.tipo == tipo]
    if s.dia.nunique() >= 3:
        linea(f"exceso medio {tipo} (n={len(s)})", bootstrap_por_dia(s, "dia", lambda d: d.exceso_vs_opt_pct.mean()), "%")

v = pd.read_csv("datos/procesado/validacion_era5_adsc.csv")
v["dia"] = pd.to_datetime(v.t, unit="s").dt.strftime("%Y-%m-%d")
v["e2"] = (v.u_era - v.u_obs) ** 2 + (v.v_era - v.v_obs) ** 2
print("=" * 100)
print(f"ERA5 FRENTE A VIENTO MEDIDO (sin dirección corrupta) — {len(v)} informes")
linea("error vectorial RMS", bootstrap_por_dia(v, "dia", lambda d: np.sqrt(d.e2.mean())), "m/s")
linea("error vectorial mediano", bootstrap_por_dia(v, "dia", lambda d: np.sqrt(d.e2).median()), "m/s")
linea("sesgo de temperatura (ERA5 − avión)", bootstrap_por_dia(v, "dia", lambda d: (d.T_era - d.T_obs).mean()), "K")
vel = lambda d: (np.hypot(d.u_era, d.v_era) - np.hypot(d.u_obs, d.v_obs))[np.hypot(d.u_obs, d.v_obs) > 60].mean()
linea("sesgo del módulo con viento > 60 m/s", bootstrap_por_dia(v, "dia", vel), "m/s")
print(f"\nB = {B} remuestreos por días, semilla {SEMILLA}.")
