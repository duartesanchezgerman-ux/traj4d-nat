"""Auditoría independiente de los cálculos del artículo A1.

Recalcula por caminos distintos a los del código principal (otra implementación, otra semilla u otro dato)
cada pieza de la cadena y comprueba que da lo mismo. No sustituye a las pruebas automáticas (tests/): las
completa con comprobaciones sobre los datos reales.

  A. Atmósfera estándar: presiones de la tabla ISA (ICAO Doc 7488) en FL300–FL410.
  B. Convenio del viento: con el signo contrario, ERA5 y los aviones no casan.
  C. Interpolación de ERA5: xarray (lineal en tiempo, log p, lat, lon) frente a la del código (scipy).
  D. Hora: duración observada de cada tramo (ADS-C) frente a la que da el modelo con el viento de ERA5.
  E. Consumo: combustible por hora de cada tipo frente a valores públicos típicos de crucero.
  F. Integrador: integración de Euler de 1 km escrita aquí frente al RK4 de 20 km del código, en 25 tramos.
  G. Optimizador: el óptimo no es peor que la ortodrómica; los casos "volada mejor que óptimo" son ruido.
  H. Descomposición: (1+E_total) = (1+E_estructural)·(1+E_elección) en cada tramo.
  I. Bootstrap: otra implementación y otra semilla dan los mismos intervalos.
  J. Cifras del manuscrito: cada número clave de main.tex (inglés) recalculado desde las tablas por tramo.

Uso:  PYTHONPATH=src python scripts/auditoria_a1.py | tee resultados/resultados_auditoria_a1.txt
Necesita los archivos de ERA5 de los días de estudio y datos/procesado/ de f2_comparacion y f3_estructura.
"""
from __future__ import annotations

import glob
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
KT, FT = 0.514444, 0.3048
FALLOS = []


def ok(nombre, cond, detalle=""):
    print(f"  [{'OK' if cond else 'FALLO'}] {nombre}" + (f" — {detalle}" if detalle else ""))
    if not cond:
        FALLOS.append(nombre)


def cab(t):
    print("\n" + "=" * 100 + f"\n{t}\n" + "=" * 100)


comp = pd.read_csv("datos/procesado/comparacion_tramos.csv")
estr = pd.read_csv("datos/procesado/estructura_tramos.csv")
val = pd.read_csv("datos/procesado/validacion_era5_adsc.csv")
dias = {Path(f).stem.replace("era5_pl_", ""): f for f in glob.glob("datos/era5/era5_pl_*.nc")}
dias = {d: f for d, f in dias.items() if d in set(comp.dia)}

# ------------------------------------------------------------------------------------------------ A
cab("A. ATMÓSFERA ESTÁNDAR (presión en hPa; tabla ISA de referencia)")
from traj4d import atmosfera as atm  # noqa: E402
tabla = {300: 300.90, 350: 238.42, 370: 216.63, 390: 196.77, 410: 178.74}   # ISA (ICAO Doc 7488), hPa
for fl, p_ref in tabla.items():
    p = float(atm.presion(atm.fl_a_m(fl))) / 100
    ok(f"FL{fl}: {p:.2f} hPa (tabla {p_ref:.2f})", abs(p - p_ref) < 0.05)

# ------------------------------------------------------------------------------------------------ B
cab("B. CONVENIO DEL VIENTO (dirección DE la que sopla: u = −V·sin, v = −V·cos)")
d = val.copy()
vel = np.hypot(d.u_obs, d.v_obs)
rms_ok = np.sqrt(((d.u_era - d.u_obs) ** 2 + (d.v_era - d.v_obs) ** 2).mean())
rms_mal = np.sqrt(((d.u_era + d.u_obs) ** 2 + (d.v_era + d.v_obs) ** 2).mean())
print(f"  error vectorial RMS (informes válidos): convenio usado {rms_ok:.2f} m/s | "
      f"convenio contrario {rms_mal:.2f} m/s | viento medio {vel.mean():.1f} m/s")
ok("el convenio usado es el correcto", rms_ok < 0.3 * rms_mal)

# ------------------------------------------------------------------------------------------------ C
cab("C. INTERPOLACIÓN DE ERA5: xarray frente al código (300 informes al azar)")
import xarray as xr  # noqa: E402
rng = np.random.default_rng(7)
m = val.sample(300, random_state=7)
m["dia"] = pd.to_datetime(m.t, unit="s").dt.strftime("%Y-%m-%d")
dif = []
for dia, g in m.groupby("dia"):
    if dia not in dias:
        continue
    ds = xr.open_dataset(dias[dia])
    ds = ds.rename({k: v for k, v in {"valid_time": "time", "pressure_level": "level"}.items() if k in ds.coords})
    ds = ds.assign_coords(logp=("level", np.log(ds.level.values * 100.0))).swap_dims(level="logp")
    p = atm.presion(atm.fl_a_m(g.fl.values))
    pts = dict(time=xr.DataArray(pd.to_datetime(g.t.values, unit="s"), dims="k"),
               logp=xr.DataArray(np.log(p), dims="k"),
               latitude=xr.DataArray(g.lat.values, dims="k"), longitude=xr.DataArray(g.lon.values, dims="k"))
    x = ds[["u", "v", "t"]].interp(**pts)
    dif.append(np.c_[x.u.values - g.u_era.values, x.v.values - g.v_era.values, x.t.values - g.T_era.values])
dif = np.vstack(dif)
mx = np.nanmax(np.abs(dif), axis=0)
print(f"  diferencia máxima xarray − código: u {mx[0]:.2e} m/s | v {mx[1]:.2e} m/s | T {mx[2]:.2e} K (n = {len(dif)})")
ok("misma interpolación", (mx < 1e-3).all())

# ------------------------------------------------------------------------------------------------ D
cab("D. HORA, POSICIONES Y VIENTO: duraciones y velocidades observadas frente al modelo")
import pyproj  # noqa: E402
from traj4d.ruta import distancia as _dist_esfera  # noqa: E402
from traj4d.vuelos import asignar_tipo, cargar_informes, segmentar_vuelos  # noqa: E402
inf = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                   pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
inf = segmentar_vuelos(inf[inf.icao24.notna()])
geod = pyproj.Geod(ellps="WGS84")
filas_d = []
for vid, g in inf[inf.vuelo_id.isin(comp.vuelo_id)].groupby("vuelo_id"):
    g = g.sort_values("t")
    la, lo = g.lat.to_numpy(), g.lon.to_numpy()
    filas_d.append(dict(vuelo_id=vid, t_obs_min=(g.t.iat[-1] - g.t.iat[0]) / 60,
                        d_elip_km=geod.line_length(lo, la) / 1e3,
                        d_esf_km=float(_dist_esfera(la[:-1], lo[:-1], la[1:], lo[1:]).sum()) / 1e3,
                        gs_tag14=g.gs_kt.median()))
c = comp.merge(pd.DataFrame(filas_d), on="vuelo_id")
c["f_elip"] = c.d_elip_km / c.d_esf_km
c["gs_pos"] = c.d_elip_km / (c.t_obs_min / 60) / 1.852
x = c.dropna(subset=["gs_tag14"])
print(f"  1) Velocidad sobre el suelo: medida a bordo (etiqueta 14) − distancia WGS84/tiempo entre informes: "
      f"mediana {(x.gs_tag14 - x.gs_pos).median():+.2f} kt (n = {len(x)})")
ok("la hora y las posiciones de los informes son coherentes con la velocidad medida a bordo",
   abs((x.gs_tag14 - x.gs_pos).median()) < 1.5)
print(f"  2) Distancias del código (esfera, R = 6371 km) frente a WGS84: las de WGS84 son un "
      f"{100 * (c.f_elip.median() - 1):.2f} % más largas (p5–p95 {100 * (c.f_elip.quantile(.05) - 1):.2f}–"
      f"{100 * (c.f_elip.quantile(.95) - 1):.2f} %)")
c["err_t_pct"] = 100 * (c.t_volada_min * c.f_elip / c.t_obs_min - 1)
for origen, g in c.groupby("origen_mach"):
    print(f"  3) Mach {origen:6s} n={len(g):3d} | duración del modelo (corregida a WGS84) − observada: mediana "
          f"{g.err_t_pct.median():+.2f} % | p5–p95 {g.err_t_pct.quantile(.05):+.2f} a {g.err_t_pct.quantile(.95):+.2f} %")
gm = c[c.origen_mach == "medido"]
ok("con Mach medido, la duración del modelo cuadra con la observada (|mediana| < 0,5 %)", abs(gm.err_t_pct.median()) < 0.5)
print(f"     casos con |error| > 5 %: {int((c.err_t_pct.abs() > 5).sum())} (informes con hora errónea no detectable "
      "al final del tramo, ver Métodos)")
print("  El factor esfera/WGS84 afecta casi por igual a la ruta volada y a la óptima: su efecto en el exceso relativo se")
print("  estimó en +0,0002 puntos porcentuales de media (p5–p95: −0,003 a +0,004) con las latitudes del óptimo restringido.")

# ------------------------------------------------------------------------------------------------ E
cab("E. CONSUMO HORARIO EN CRUCERO frente a valores públicos típicos (kg/h, a 0,8·MTOW)")
tipicos = {"B789": (5000, 6500), "A359": (5300, 7000), "B77W": (6800, 9000), "B772": (6000, 8000),
           "A333": (5200, 7000), "B788": (4600, 6000), "B78X": (5300, 7000), "A35K": (6200, 8200)}
c["kg_h"] = c.comb_volada / c.t_volada_min * 60
for tipo, (a, b) in tipicos.items():
    g = c[c.tipo == tipo]
    if len(g):
        v = g.kg_h.median()
        ok(f"{tipo}: {v:,.0f} kg/h (n={len(g)}; rango típico {a:,}–{b:,})", a <= v <= b)

# ------------------------------------------------------------------------------------------------ F
cab("F. INTEGRADOR: Euler de 1 km (escrito aquí) frente al RK4 de 20 km del código, 25 tramos")
from traj4d.rendimiento import PollSchumannModel  # noqa: E402
from traj4d.viento import CampoERA5  # noqa: E402
from traj4d.ruta import distancia  # noqa: E402
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params  # noqa: E402
P = load_aircraft_engine_params()
cf = c[c.dia.isin(dias)].sort_values("vuelo_id")          # solo días con ERA5 en esta máquina
muestra = cf.iloc[::max(1, len(cf) // 25)].head(25)
filas_f, modelos, campo_dia = [], {}, None
for r in muestra.sort_values("dia").itertuples():
    if campo_dia != r.dia:
        campo, campo_dia = CampoERA5(dias[r.dia]), r.dia
    rep = inf[inf.vuelo_id == r.vuelo_id].sort_values("t")
    lat, lon = rep.lat.to_numpy(), rep.lon.to_numpy()
    sel_ = [0] + [i for i in range(1, len(lat) - 1)
                  if float(distancia(lat[i], lon[i], lat[i - 1], lon[i - 1])) > 5e3
                  and float(distancia(lat[i], lon[i], lat[-1], lon[-1])) > 5e3] + [len(lat) - 1]
    lat, lon = lat[sel_], lon[sel_]
    mod = modelos.setdefault(r.tipo, PollSchumannModel(r.tipo))
    m_, t_ = 0.80 * P[r.tipo].amass_mtow, float(rep.t.iloc[0])
    m0 = m_
    for i in range(len(lat) - 1):                    # tramo a tramo entre informes, en ortodrómica
        L = float(distancia(lat[i], lon[i], lat[i + 1], lon[i + 1]))
        n = max(1, int(np.ceil(L / 1000.0)))
        for k in range(n):
            f = (k + 0.5) / n
            la = lat[i] + f * (lat[i + 1] - lat[i]); lo = lon[i] + f * (lon[i + 1] - lon[i])
            # rumbo de la derrota local (aprox. plana, suficiente en 1 km)
            chi = np.arctan2((lon[i + 1] - lon[i]) * np.cos(np.radians(la)), lat[i + 1] - lat[i])
            u, v, dT = campo.muestra(t_, r.fl, la, lo)
            tas = float(atm.mach_a_tas(r.mach, atm.fl_a_m(r.fl), dT))
            wl, wt = u * np.sin(chi) + v * np.cos(chi), u * np.cos(chi) - v * np.sin(chi)
            vg = wl + np.sqrt(tas ** 2 - wt ** 2)
            dt = (L / n) / vg
            m_ -= float(mod.ff(m_, r.mach, r.fl, dT)) * dt
            t_ += dt
    filas_f.append((r.vuelo_id, m0 - m_, r.comb_volada))
F = pd.DataFrame(filas_f, columns=["vuelo_id", "euler", "codigo"])
F["dif_pct"] = 100 * (F.euler / F.codigo - 1)
print(f"  diferencia de combustible Euler − código: media {F.dif_pct.mean():+.3f} % | máx |·| {F.dif_pct.abs().max():.3f} %")
ok("los dos integradores coinciden (< 0,3 %; la ruta del código interpola en ortodrómica por tramos)",
   F.dif_pct.abs().max() < 0.3)

# ------------------------------------------------------------------------------------------------ G
cab("G. OPTIMIZADOR")
ok("el óptimo nunca es peor que la ortodrómica (tolerancia 0,05 %)", (comp.ganancia_opt_vs_gc_pct > -0.05).all(),
   f"mínimo {comp.ganancia_opt_vs_gc_pct.min():+.3f} %")
neg = comp[comp.exceso_vs_opt_pct < 0]
print(f"  tramos con la volada mejor que el óptimo: {len(neg)} de {len(comp)}; "
      f"mediana {neg.exceso_vs_opt_pct.median():+.3f} %, mínimo {neg.exceso_vs_opt_pct.min():+.3f} %")
print(f"  de ellos, por debajo de −0,2 % (suelo de ruido declarado): {(neg.exceso_vs_opt_pct < -0.2).sum()}")
ok("casi todos los casos negativos están dentro del suelo de ruido (±0,2 %)",
   (neg.exceso_vs_opt_pct < -0.2).sum() <= 0.02 * len(comp))
ok("la ruta volada nunca es más corta que la ortodrómica", (comp.dist_volada_km >= comp.dist_gc_km - 0.5).all())

# ------------------------------------------------------------------------------------------------ H
cab("H. DESCOMPOSICIÓN ESTRUCTURA × ELECCIÓN")
lhs = 1 + estr.exceso_total_pct / 100
rhs = (1 + estr.estructural_pct / 100) * (1 + estr.eleccion_pct / 100)
ok("(1+E_total) = (1+E_estr)·(1+E_elec) en todos los tramos", np.allclose(lhs, rhs, atol=1e-9),
   f"máx. error {np.abs(lhs - rhs).max():.1e}")
x = estr.merge(comp[["vuelo_id", "exceso_vs_opt_pct"]], on="vuelo_id")
ok("el exceso total de f3_estructura es el mismo que el de f2_comparacion",
   np.allclose(x.exceso_total_pct, x.exceso_vs_opt_pct, atol=1e-9))
ok("el óptimo restringido nunca es mejor que el libre (tolerancia 0,2 %)", (estr.estructural_pct > -0.2).all(),
   f"mínimo {estr.estructural_pct.min():+.3f} %")

# ------------------------------------------------------------------------------------------------ I
cab("I. BOOTSTRAP POR DÍAS: implementación y semilla distintas (B = 20 000, semilla 12345)")
r2 = np.random.default_rng(12345)
grupos = [g.exceso_vs_opt_pct.to_numpy() for _, g in comp.groupby("dia")]
est = []
for _ in range(20000):
    idx = r2.integers(0, len(grupos), len(grupos))
    est.append(np.concatenate([grupos[i] for i in idx]).mean())
lo, hi = np.percentile(est, [2.5, 97.5])
import re as _re
_m = _re.search(r"exceso de combustible, media\s+([\d.]+) %\s+IC95 \[\s*([\d.]+),\s*([\d.]+)\]",
                Path("resultados/resultados_f3_estadistica.txt").read_text(encoding="utf-8"))
ref_lo, ref_hi = float(_m.group(2)), float(_m.group(3))
print(f"  exceso medio {comp.exceso_vs_opt_pct.mean():.3f} %  IC95 [{lo:.3f}, {hi:.3f}]  "
      f"(f3_estadistica, B = 5000, semilla 20261102: [{ref_lo:.3f}, {ref_hi:.3f}])")
ok("los dos bootstraps dan el mismo intervalo a ±0,01", abs(lo - ref_lo) < 0.01 and abs(hi - ref_hi) < 0.01)

# ------------------------------------------------------------------------------------------------ J
cab("J. CIFRAS DEL MANUSCRITO (inglés): cada frase se reconstruye con los datos y se busca literal en main.tex")
tex = Path("articulo/en/main.tex").read_text(encoding="utf-8")
c["estacion"] = np.where(pd.to_datetime(c.dia).dt.month.isin([11, 12, 1, 2]), "W", "S")
kg = c.comb_volada - c.comb_opt
vv = val.copy()          # la tabla de validación solo guarda los informes válidos
e = np.hypot(vv.u_era - vv.u_obs, vv.v_era - vv.v_obs)
vel = np.hypot(vv.u_obs, vv.v_obs)


def _dif(fila):          # otra implementación: diferencia máxima con el plan llevado a la rejilla de 0,5°
    if not isinstance(fila.lat_plan, str) or not isinstance(fila.lat_rest, str):
        return np.nan
    a = [round(float(v) * 2) / 2 for v in fila.lat_plan.split(";")]
    b = [float(v) for v in fila.lat_rest.split(";")]
    return max(abs(x - y) for x, y in zip(a, b)) if len(a) == len(b) else np.nan


pc = estr.assign(dif=estr.apply(_dif, axis=1)).dropna(subset=["dif"])
cl = np.select([pc.dif < 0.25, pc.dif < 0.75], ["igual", "medio"], "uno")
m1 = lambda x: f"{x:,}"                       # 1,885
f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:.3f}"
pct0 = lambda x: f"{100 * x:.0f}"


def fila(g, ci=None):
    return (f"{m1(len(g))} & {f2(g.exceso_vs_opt_pct.mean())}" + (f" ({ci})" if ci else "")
            + f" & {f2(g.exceso_vs_opt_pct.median())} & {(g.comb_volada - g.comb_opt).mean():.0f}")


def grupo(est, sen):
    return c[(c.estacion == est) & (c.sentido == sen)]


sens = c.dropna(subset=["exceso_vs_opt_m70"])
clases = pd.cut(vel, [0, 20, 40, 60, 120], include_lowest=True)   # la clase baja incluye la calma (0 m/s)
tab_v = {k: e[clases == k] for k in clases.cat.categories}
afirmaciones = [
    f"{m1(len(c))} were evaluated",
    f"burns {f2(c.exceso_vs_opt_pct.mean())}\\% more fuel than the wind-optimal path",
    f"or {kg.mean():.0f}~kg",
    f"on segments of {c.dist_volada_km.mean():,.0f}~km on average",
    f"the median is {f2(c.exceso_vs_opt_pct.median())}\\%",
    f"{pct0((c.exceso_vs_opt_pct > .5).mean())}\\% of the segments exceed 0.5\\% and {pct0((c.exceso_vs_opt_pct > 1).mean())}\\% exceed 1\\%",
    f"The optimum is on average {f2(c.ganancia_opt_vs_gc_pct.mean())}\\%",
    f"Winter, eastbound & {fila(grupo('W', 'E'))} \\\\",
    f"Winter, westbound & {fila(grupo('W', 'O'))} \\\\",
    f"Summer, eastbound & {fila(grupo('S', 'E'))} \\\\",
    f"Summer, westbound & {fila(grupo('S', 'O'))} \\\\",
    f"Boeing 787-9 & {m1((c.tipo == 'B789').sum())} & {f2(c[c.tipo == 'B789'].exceso_vs_opt_pct.mean())}",
    f"Airbus A350-900 & {m1((c.tipo == 'A359').sum())} & {f2(c[c.tipo == 'A359'].exceso_vs_opt_pct.mean())}",
    f"subsample ({m1(len(sens))} segments, one in five), the mean excess is {f3(sens.exceso_vs_opt_m70.mean())}\\%, "
    f"{f3(sens.exceso_vs_opt_pct.mean())}\\% and {f3(sens.exceso_vs_opt_m90.mean())}\\%",
    f"Mach number was measured on the flight itself in {m1((c.origen_mach == 'medido').sum())} segments",
    f"for only {pct0((c.origen_mach == 'medido').mean())}\\% of the segments",
    f"comparison of {m1(len(vv))} valid reports",
    f"The vector root-mean-square error is {f2(np.sqrt((e ** 2).mean()))}~m/s",
    f"the median vector error {f2(e.median())}~m/s",
    f"the 95th percentile {f2(e.quantile(.95))}~m/s",
    f"ERA5 is {f2(-(vv.T_era - vv.T_obs).mean())}~K colder than the aircraft",
    f"costs on average only {f3(estr.estructural_pct.mean())}\\%",
    f"that is, {100 * estr.estructural_pct.mean() / estr.exceso_total_pct.mean():.0f}\\% of the total excess",
    f"Of the {m1(len(pc))} segments whose planned waypoints were identified",
    f"Identical on every meridian & {m1((cl == 'igual').sum())}",
    f"Differs by 0.5$^\\circ$      & {m1((cl == 'medio').sum())}",
    f"Differs by $\\geq$1$^\\circ$  & {m1((cl == 'uno').sum())}",
    f"from the {pct0((cl != 'igual').mean())}\\% of segments with identified planned waypoints",
]
for k, x in tab_v.items():            # tabla de ERA5 por clase de viento: n y RMS
    afirmaciones.append(f"& {m1(len(x))}")
    afirmaciones.append(f"& {f2(np.sqrt((x ** 2).mean()))} &")
for texto in afirmaciones:
    ok(texto[:90], texto in tex, "" if texto in tex else "NO APARECE en main.tex")

cab("RESUMEN")
print(f"  Comprobaciones con fallo: {len(FALLOS)}" + (": " + "; ".join(FALLOS) if FALLOS else " (ninguna)"))
