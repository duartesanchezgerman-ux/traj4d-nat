"""Fase 2: primera comparación entre la ruta volada y la ruta óptima, tramo a tramo (días de estudio).
Ejecutar:  PYTHONPATH=src python scripts/f2_comparacion.py | tee resultados/resultados_f2_comparacion.txt

Unidad de análisis (decisión del 2 nov 2026): el tramo oceánico observado, con ≥15° de longitud y
≥4 informes. Para cada tramo:
  * Ruta volada: la poligonal que une los informes ADS-C en orden temporal.
  * Ruta óptima: programación dinámica lateral entre el primer y el último informe, a la misma hora
    de inicio, mismo FL (mediana del tramo) y mismo Mach (mediana medida del vuelo; si falta, la
    mediana medida de su tipo).
  * Referencia: la ortodrómica entre los mismos extremos.
  * Las tres rutas se evalúan con el mismo integrador RK4, el viento 4D de ERA5 y el modelo
    Poll–Schumann del tipo de avión.
Hipótesis declaradas: masa = 0,80·MTOW, con sensibilidad a 0,70 y 0,90 en el resumen; crucero a FL
constante (se descartan los tramos con más de 4000 ft de variación); CI irrelevante, porque a Mach fijo
el combustible y el tiempo varían en el mismo sentido.
"""
import glob
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params

from traj4d.crucero import integrar_ruta, perfil_constante
from traj4d.optimizador import optimizar_lateral
from traj4d.rendimiento import PollSchumannModel
from traj4d.ruta import Ruta, distancia, puntos_distintos
from traj4d.viento import CampoERA5
from traj4d.vuelos import asignar_tipo, cargar_informes, resumen_vuelos, segmentar_vuelos

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
PARAMS = load_aircraft_engine_params()
DOM = dict(lat=(30.5, 69.5), lon=(-59.5, -0.5))          # dominio de ERA5 descargado, con margen

df = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                  pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
df = segmentar_vuelos(df[df.icao24.notna()])
v = resumen_vuelos(df)
mach_tipo = df.dropna(subset=["mach"]).groupby("typecode")["mach"].median()
dias = {Path(f).stem.replace("era5_pl_", ""): f for f in glob.glob("datos/era5/era5_pl_*.nc")}
dias = filtrar_dias_estudio(dias)
v["dia"] = v.fecha.astype(str)
sel = v[(v.span_lon >= 15) & (v.n >= 4) & v.typecode.isin(PARAMS.keys()) & v.dia.isin(dias)
        & (v.lat_min >= DOM["lat"][0]) & (v.lat_max <= DOM["lat"][1])
        & (v.lon_min >= DOM["lon"][0]) & (v.lon_max <= DOM["lon"][1])]

# Un solo día de ERA5 en memoria cada vez (con 70 días no caben todos) y tramos ordenados por día.
_campo = {}


def campo_de(dia):
    if dia not in _campo:
        _campo.clear()
        _campo[dia] = CampoERA5(dias[dia])
    return _campo[dia]


sel = sel.sort_values("dia", kind="stable")
modelos, filas = {}, []
descartes = {"fl_variable": 0, "fl_fuera_era5": 0, "fuera_tiempo": 0, "no_es_cruce": 0, "rejilla_fuera_dominio": 0}
for vid, s in sel.iterrows():
    r = df[df.vuelo_id == vid].sort_values("t")
    fl_rango = (r.alt_ft.max() - r.alt_ft.min()) / 100
    fl = float(np.round(r.alt_ft.median() / 100))
    if fl_rango > 40:
        descartes["fl_variable"] += 1; continue
    if not 270 <= fl <= 440:
        descartes["fl_fuera_era5"] += 1; continue
    tipo = s.typecode
    mach = float(r.mach.median()) if r.mach.notna().any() else float(mach_tipo.get(tipo, 0.84))
    origen_mach = "medido" if r.mach.notna().any() else "tipo"
    modelos.setdefault(tipo, PollSchumannModel(tipo))
    mod, campo = modelos[tipo], campo_de(s.dia)
    t0 = float(r.t.iloc[0])
    if r.t.iloc[-1] > campo.limites["t"][1] - 7200:      # el vuelo debe acabar antes de las 22 UTC del día del archivo
        descartes["fuera_tiempo"] += 1; continue
    lat, lon = r.lat.to_numpy(), r.lon.to_numpy()
    keep = puntos_distintos(lat, lon)
    volada = Ruta(lat[keep], lon[keep])
    gc = Ruta([lat[0], lat[-1]], [lon[0], lon[-1]])
    # Control de calidad: un tramo de cruce avanza de forma monótona en longitud (retrocesos <= 0,5°)
    # y su distancia no supera en más de un 5 % la ortodrómica. Si no, son idas y vueltas, esperas o
    # dos vuelos del mismo avión unidos por la segmentación.
    dlon = np.diff(lon) * (1 if lon[-1] > lon[0] else -1)
    if dlon.min() < -0.5 or volada.longitud > 1.05 * gc.longitud:
        descartes["no_es_cruce"] += 1; continue
    prof = perfil_constante(fl, mach)
    try:
        res = {}
        masas = (0.70, 0.80, 0.90) if int(vid) % 5 == 0 else (0.80,)    # sensibilidad en 1 de cada 5 tramos (por id, no depende del orden)
        for frac in masas:
            m0 = frac * PARAMS[tipo].amass_mtow
            n_et = max(8, int(round(gc.longitud / 100e3)))
            dp = optimizar_lateral(mod, campo, lat[0], lon[0], lat[-1], lon[-1], fl, mach, m0, t0=t0,
                                   n_etapas=n_et, K=20, delta_km=15.0, salto_max=3)
            ev = {k: integrar_ruta(mod, m0, rr, prof, campo, t0=t0, ds=20e3)
                  for k, rr in (("volada", volada), ("optima", dp.ruta), ("gc", gc))}
            res[frac] = (ev, dp)
    except ValueError as e:
        if "out of bounds" not in str(e):
            raise
        descartes["rejilla_fuera_dominio"] += 1; continue
    ev, dp = res[0.80]
    filas.append(dict(
        vuelo_id=vid, dia=s.dia, tipo=tipo, indicativo=s.indicativo, sentido=s.sentido, fl=fl, mach=mach,
        origen_mach=origen_mach, n_informes=len(r), dist_volada_km=volada.longitud / 1e3,
        dist_gc_km=gc.longitud / 1e3, dist_opt_km=dp.ruta.longitud / 1e3,
        desv_max_opt_km=float(np.abs(dp.desplazamientos_km).max()),
        comb_volada=ev["volada"].combustible, comb_opt=ev["optima"].combustible, comb_gc=ev["gc"].combustible,
        t_volada_min=ev["volada"].t[-1] / 60, t_opt_min=ev["optima"].t[-1] / 60,
        exceso_vs_opt_pct=100 * (ev["volada"].combustible / ev["optima"].combustible - 1),
        exceso_vs_gc_pct=100 * (ev["volada"].combustible / ev["gc"].combustible - 1),
        ganancia_opt_vs_gc_pct=100 * (1 - ev["optima"].combustible / ev["gc"].combustible),
        **{f"exceso_vs_opt_m{int(100 * k)}": (100 * (res[k][0]["volada"].combustible / res[k][0]["optima"].combustible - 1)
                                              if k in res else np.nan) for k in (0.70, 0.90)},
    ))

t = pd.DataFrame(filas)
t.to_csv(Path(os.environ.get("SALIDA", "datos/procesado")) / "comparacion_tramos.csv", index=False)  # SALIDA: ejemplo de un día
pd.set_option("display.width", 200)
print("=" * 90)
print(f"Tramos candidatos en los días de estudio: {len(sel)} | evaluados: {len(t)} | descartes: {descartes}")
print(f"Mach medido en el propio vuelo: {(t.origen_mach == 'medido').sum()} | Mach de su tipo: {(t.origen_mach == 'tipo').sum()}")
print("\nCOMPROBACIÓN: el óptimo no debe ser peor que la volada ni que la ortodrómica (salvo la resolución de la rejilla)")
print(f"  óptimo > volada en {(t.exceso_vs_opt_pct < -0.05).sum()} tramos (mínimo {t.exceso_vs_opt_pct.min():.3f} %) | "
      f"óptimo > ortodrómica en {(t.ganancia_opt_vs_gc_pct < -0.05).sum()} tramos")
print("  -> Suelo de ruido numérico (rejilla de 15 km + aproximación de la DP): los excesos en |x| < 0,2 % no son distinguibles de 0.")
print("\nDISTRIBUCIÓN del exceso de la volada frente al óptimo [%]:")
print("  " + " | ".join(f"p{int(100 * q)} {t.exceso_vs_opt_pct.quantile(q):.3f}" for q in (.1, .25, .5, .75, .9, .95, .99)))
print(f"  fracción > 0,5 %: {100 * (t.exceso_vs_opt_pct > .5).mean():.1f} % | > 1 %: {100 * (t.exceso_vs_opt_pct > 1).mean():.1f} %")


def resumen(g):
    return pd.Series({"n": len(g), "exceso_vs_opt_media": g.exceso_vs_opt_pct.mean(),
                      "exceso_vs_opt_p50": g.exceso_vs_opt_pct.median(),
                      "exceso_vs_opt_p90": g.exceso_vs_opt_pct.quantile(.9),
                      "exceso_vs_gc_media": g.exceso_vs_gc_pct.mean(),
                      "ganancia_opt_vs_gc_media": g.ganancia_opt_vs_gc_pct.mean(),
                      "kg_exceso_media": (g.comb_volada - g.comb_opt).mean(),
                      "dist_km_media": g.dist_volada_km.mean()})


print("\nGLOBAL"); print(resumen(t).round(3).to_string())
t["estacion"] = np.where(pd.to_datetime(t.dia).dt.month.isin([11, 12, 1, 2]), "invierno", "verano")
print("\nPOR ESTACIÓN Y SENTIDO"); print(t.groupby(["estacion", "sentido"]).apply(resumen).round(3).to_string())
print("\nTIPOS PRINCIPALES"); print(t[t.tipo.isin(["B789", "A359", "B77W", "B772", "A333"])].groupby("tipo").apply(resumen).round(3).to_string())
sub = t.dropna(subset=["exceso_vs_opt_m70"])
print(f"\nSENSIBILIDAD A LA MASA (submuestra de {len(sub)} tramos; exceso medio de la volada frente al óptimo, %)")
print(f"  0,70·MTOW: {sub.exceso_vs_opt_m70.mean():.3f} | 0,80·MTOW: {sub.exceso_vs_opt_pct.mean():.3f} | 0,90·MTOW: {sub.exceso_vs_opt_m90.mean():.3f}")
