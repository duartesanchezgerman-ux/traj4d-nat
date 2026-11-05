"""Fase 3: ¿cuánta de la ineficiencia lateral es ESTRUCTURAL?
Ejecutar (después de f2_comparacion.py):
    PYTHONPATH=src python scripts/f3_estructura.py | tee resultados/resultados_f3_estructura.txt

Hipótesis (esquema del artículo, 2 nov 2026): en el espacio oceánico del Atlántico Norte las rutas se planifican
con puntos de paso definidos por latitud y longitud en cada meridiano múltiplo de 10° (NAT Doc 007), con la
latitud en grados enteros o medios grados. Entre puntos de paso el avión vuela el arco de círculo máximo.
El óptimo de la programación dinámica (f2_comparacion) no tiene esa restricción. Descomponemos:

    exceso total    = volada / óptimo libre − 1
    estructural     = óptimo restringido / óptimo libre − 1     (coste de la rejilla de puntos de paso)
    de elección     = volada / óptimo restringido − 1           (lo que se podía ganar dentro de la rejilla)

Óptimo restringido: ruta por puntos en cada meridiano múltiplo de 10° interior al tramo (a más de 0,25° de los
extremos), con latitud en pasos de 0,5°, entre los mismos extremos que la volada. Se busca por descenso
coordenado (±0,5° por meridiano hasta que no mejora), con dos arranques: la latitud de cruce de la volada y la
de la ortodrómica, redondeadas a 0,5°. Se evalúa con el mismo integrador RK4, viento ERA5 4D, modelo, FL, Mach,
masa (0,80·MTOW), hora de inicio y paso (20 km) que f2_comparacion. Como el óptimo libre sale de una rejilla
discreta (15 km), el óptimo restringido puede salir algo mejor que él: eso es el suelo de ruido (±0,2 %).

Además se comprueba con el Tag 13 del ADS-C (ruta prevista: siguiente punto y el posterior) qué granularidad
tienen los puntos de paso planificados en los meridianos múltiplos de 10°.
"""
import glob
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params

from traj4d.crucero import integrar_ruta, perfil_constante
from traj4d.rendimiento import PollSchumannModel
from traj4d.ruta import Ruta, distancia, puntos_distintos
from traj4d.viento import CampoERA5
from traj4d.vuelos import asignar_tipo, cargar_informes, segmentar_vuelos

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
TOL = 0.05            # tolerancia [°] para decidir si un punto está en un meridiano o en una latitud de 0,5°
PASO = 0.5            # granularidad de latitud de la rejilla de puntos de paso

df = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                  pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
df = segmentar_vuelos(df[df.icao24.notna()])
comp = pd.read_csv("datos/procesado/comparacion_tramos.csv")
dias = {Path(f).stem.replace("era5_pl_", ""): f for f in glob.glob("datos/era5/era5_pl_*.nc")}
dias = filtrar_dias_estudio(dias)


def clase_lat(lat):
    if abs(lat - round(lat)) < TOL:
        return "entero"
    if abs(2 * lat - round(2 * lat)) < 2 * TOL:
        return "medio"
    return "otro"


# --- 1. Granularidad de los puntos de paso PLANIFICADOS (Tag 13) en todo el NAT -------------------------------
wp = pd.concat([df[["wp1_lat", "wp1_lon"]].set_axis(["lat", "lon"], axis=1),
                df[["wp2_lat", "wp2_lon"]].set_axis(["lat", "lon"], axis=1)]).dropna()
wp = wp[(wp.lat.abs() <= 90) & (wp.lon.abs() <= 180) & wp.lon.between(-60, -10) & wp.lat.between(30, 70)]
en_meridiano = (wp.lon - 10 * np.round(wp.lon / 10)).abs() < TOL
wm = wp[en_meridiano].round(2).drop_duplicates()
cl = wm.lat.map(clase_lat).value_counts(normalize=True) * 100
print("=" * 100)
print(f"Puntos de paso planificados (Tag 13) distintos entre 10°W y 60°W: {len(wp.round(2).drop_duplicates())}"
      f" | en un meridiano múltiplo de 10°: {len(wm)}")
print("  latitud de los que están en un meridiano múltiplo de 10° [%]: "
      + " | ".join(f"{k} {v:.1f}" for k, v in cl.items()))

# --- 2. Óptimo restringido por tramo (se omite con --resumen: solo se rehace el resumen desde el CSV) ---------
SOLO_RESUMEN = "--resumen" in sys.argv
modelos, filas = {}, []
cache = {}
_campo = {}


def campo_de(dia):
    """Un solo día de ERA5 en memoria cada vez; los tramos se recorren ordenados por día."""
    if dia not in _campo:
        _campo.clear()
        _campo[dia] = CampoERA5(dias[dia])
    return _campo[dia]


def coste(mod, m0, lat_pts, lon_pts, prof, campo, t0, clave):
    if clave not in cache:
        cache[clave] = integrar_ruta(mod, m0, Ruta(lat_pts, lon_pts), prof, campo, t0=t0, ds=20e3).combustible
    return cache[clave]


for _, c in ([] if SOLO_RESUMEN else comp.sort_values("dia", kind="stable").iterrows()):
    r = df[df.vuelo_id == c.vuelo_id].sort_values("t")
    lat, lon = r.lat.to_numpy(), r.lon.to_numpy()
    keep = puntos_distintos(lat, lon)
    lat, lon = lat[keep], lon[keep]
    sgn = 1 if lon[-1] > lon[0] else -1
    lo, hi = sorted((lon[0], lon[-1]))
    meridianos = [m for m in range(-60, 1, 10) if lo + 0.25 < m < hi - 0.25]
    meridianos = sorted(meridianos, key=lambda m: sgn * m)
    if not meridianos:
        continue
    tipo = c.tipo
    modelos.setdefault(tipo, PollSchumannModel(tipo))
    mod, campo = modelos[tipo], campo_de(c.dia)
    m0 = 0.80 * PARAMS[tipo].amass_mtow
    t0 = float(r.t.iloc[0])
    prof = perfil_constante(c.fl, c.mach)
    cache.clear()

    # latitud de cruce de la volada en cada meridiano (interpolación en longitud a lo largo de la poligonal)
    orden = np.argsort(sgn * lon)
    lat_cruce_vol = [float(np.interp(sgn * m, sgn * lon[orden], lat[orden])) for m in meridianos]
    lat_cruce_gc = [float(np.interp(sgn * m, [sgn * lon[0], sgn * lon[-1]], [lat[0], lat[-1]])) for m in meridianos]
    # Nota: la interpolación lineal de la ortodrómica en longitud es aproximada; solo es un punto de arranque.

    # puntos planificados en esos meridianos (Tag 13, del propio tramo)
    pl = pd.concat([r[["wp1_lat", "wp1_lon"]].set_axis(["lat", "lon"], axis=1),
                    r[["wp2_lat", "wp2_lon"]].set_axis(["lat", "lon"], axis=1)]).dropna()
    plan = {}
    for m in meridianos:
        q = pl[(pl.lon - m).abs() < TOL]
        if len(q):
            plan[m] = float(q.lat.round(2).mode().iloc[0])

    def evalua(lats):
        lats = tuple(round(x * 2) / 2 for x in lats)
        return coste(mod, m0, [lat[0], *lats, lat[-1]], [lon[0], *meridianos, lon[-1]], prof, campo, t0, lats), lats

    try:
        mejor = (np.inf, None)
        for arranque in (lat_cruce_vol, lat_cruce_gc):
            f, x = evalua(arranque)
            mejora = True
            while mejora:
                mejora = False
                for i in range(len(x)):
                    for d in (-PASO, PASO):
                        y = list(x); y[i] += d
                        if not 30.5 <= y[i] <= 69.5:
                            continue
                        fy, ty = evalua(y)
                        if fy < f - 1e-6:
                            f, x, mejora = fy, ty, True
            if f < mejor[0]:
                mejor = (f, x)
        f_plan = np.nan
        if len(plan) == len(meridianos):
            f_plan = coste(mod, m0, [lat[0], *plan.values(), lat[-1]], [lon[0], *meridianos, lon[-1]],
                           prof, campo, t0, ("plan",))
    except ValueError as e:
        if "out of bounds" not in str(e):
            raise
        continue
    f_rest, x_rest = mejor
    filas.append(dict(
        vuelo_id=c.vuelo_id, dia=c.dia, tipo=tipo, sentido=c.sentido, n_meridianos=len(meridianos),
        comb_volada=c.comb_volada, comb_opt=c.comb_opt, comb_rest=f_rest, comb_plan=f_plan,
        exceso_total_pct=c.exceso_vs_opt_pct,
        estructural_pct=100 * (f_rest / c.comb_opt - 1),
        eleccion_pct=100 * (c.comb_volada / f_rest - 1),
        volada_vs_plan_pct=100 * (c.comb_volada / f_plan - 1) if np.isfinite(f_plan) else np.nan,
        lat_plan=";".join(f"{v:.2f}" for v in plan.values()),
        clases_plan=";".join(clase_lat(v) for v in plan.values()),
        lat_rest=";".join(f"{v:.1f}" for v in x_rest),
        dif_plan_rest_grados=(float(np.max(np.abs(np.array(list(plan.values())) - np.array(x_rest))))
                              if len(plan) == len(meridianos) else np.nan),
    ))

if SOLO_RESUMEN:
    t = pd.read_csv("datos/procesado/estructura_tramos.csv", keep_default_na=False, na_values=[""])
    for col in ("lat_plan", "clases_plan", "lat_rest"):
        t[col] = t[col].astype(str).replace("nan", "")
else:
    t = pd.DataFrame(filas)
    t.to_csv("datos/procesado/estructura_tramos.csv", index=False)
pd.set_option("display.width", 200)
print("=" * 100)
print(f"Tramos de f2_comparacion: {len(comp)} | con al menos un meridiano múltiplo de 10° interior: {len(t)}"
      f" | con todos los puntos planificados identificados (Tag 13): {t.comb_plan.notna().sum()}")
print(f"  meridianos interiores por tramo: {t.n_meridianos.value_counts().sort_index().to_dict()}")
clases = pd.Series(";".join(t.clases_plan[t.clases_plan != ""]).split(";")).value_counts(normalize=True) * 100
print("  latitud de los puntos planificados de estos tramos [%]: " + " | ".join(f"{k} {v:.1f}" for k, v in clases.items()))
print("\nDESCOMPOSICIÓN DEL EXCESO (media y mediana, %)")
for col, nom in (("exceso_total_pct", "total (volada / óptimo libre)"),
                 ("estructural_pct", "estructural (restringido / óptimo libre)"),
                 ("eleccion_pct", "de elección (volada / restringido)"),
                 ("volada_vs_plan_pct", "volada / ruta planificada (Tag 13)")):
    x = t[col].dropna()
    print(f"  {nom:45s} media {x.mean():7.3f} | mediana {x.median():7.3f} | p10 {x.quantile(.1):7.3f}"
          f" | p90 {x.quantile(.9):7.3f} | n {len(x)}")
def dif_plan_rest(fila):
    """Máxima diferencia [°] entre el punto planificado (redondeado a la rejilla de 0,5°) y el del óptimo restringido."""
    if not fila.lat_plan or not fila.lat_rest:
        return np.nan
    p = np.round(2 * np.array([float(v) for v in fila.lat_plan.split(";")])) / 2
    q = np.array([float(v) for v in fila.lat_rest.split(";")])
    return float(np.max(np.abs(p - q))) if len(p) == len(q) else np.nan


ok = t.apply(dif_plan_rest, axis=1).dropna()
print(f"\nEl punto planificado coincide con el del óptimo restringido en todos los meridianos: {100 * (ok == 0).mean():.1f} %"
      f" | difiere en 0,5°: {100 * (ok == 0.5).mean():.1f} % | en ≥1°: {100 * (ok >= 1).mean():.1f} %  (n = {len(ok)})")
t["dif_plan_rest"] = t.apply(dif_plan_rest, axis=1)
t["clase_dif"] = pd.cut(t.dif_plan_rest, [-0.1, 0.1, 0.6, 10], labels=["igual", "0,5°", ">=1°"])
print("\nSEGÚN LA DIFERENCIA ENTRE EL PUNTO PLANIFICADO Y EL DEL ÓPTIMO RESTRINGIDO (máxima en el tramo)")
print(t.groupby("clase_dif", observed=True)[["exceso_total_pct", "estructural_pct", "eleccion_pct"]]
      .agg(["mean", "median", "count"]).round(3).to_string())
s_ = t.dropna(subset=["dif_plan_rest"])
print(f"  parte de la suma del exceso de elección que aportan los tramos con >=1°: "
      f"{100 * s_[s_.dif_plan_rest >= 1].eleccion_pct.sum() / s_.eleccion_pct.sum():.0f} %")
rng = np.random.default_rng(20261102)
dias_u = t.dia.unique()
grupos = {d: t[t.dia == d] for d in dias_u}
boot = {k: [] for k in ("exceso_total_pct", "estructural_pct", "eleccion_pct")}
for _ in range(5000):
    m = pd.concat([grupos[d] for d in rng.choice(dias_u, len(dias_u), replace=True)])
    for k in boot:
        boot[k].append(m[k].mean())
print(f"\nIC95 bootstrap por días (B = 5000, semilla 20261102), media:")
for k, v in boot.items():
    print(f"  {k:20s} {t[k].mean():.3f} %  [{np.percentile(v, 2.5):.3f}, {np.percentile(v, 97.5):.3f}]")
print(f"  fracción estructural de la media del exceso: {100 * t.estructural_pct.mean() / t.exceso_total_pct.mean():.0f} %")
print("\nPOR SENTIDO")
print(t.groupby("sentido")[["exceso_total_pct", "estructural_pct", "eleccion_pct"]].agg(["mean", "median", "count"]).round(3).to_string())
