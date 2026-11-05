"""Fase 3: figuras del artículo 1 (PDF vectorial + PNG 300 ppp), en inglés o en español.
Ejecutar:  PYTHONPATH=src python scripts/f3_figuras.py --idioma en   -> articulo/en/figures/
           PYTHONPATH=src python scripts/f3_figuras.py --idioma es   -> articulo/es/figures/
(Los nombres de archivo son los mismos en los dos idiomas; solo cambian los textos y el separador decimal.)

Estilo: paleta categórica validada (máx. 3 series por figura, validada en todos los pares para
daltonismo), rampa secuencial azul de un solo tono para magnitudes, texto en tinta neutra, rejilla
recesiva, marcas finas, leyenda + etiquetas directas selectivas. Se regeneran solas con los datos actualizados.
"""
import argparse
import warnings
from pathlib import Path

import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params

from traj4d.crucero import optimo_local
from traj4d.optimizador import optimizar_lateral
from traj4d.rendimiento import PollSchumannModel
from traj4d.ruta import Ruta, distancia
from traj4d.viento import CampoERA5
from traj4d.vuelos import asignar_tipo, cargar_informes, segmentar_vuelos

warnings.filterwarnings("ignore")
_ap = argparse.ArgumentParser(); _ap.add_argument("--idioma", choices=("en", "es"), default="en")
IDIOMA = _ap.parse_args().idioma
OUT = Path(f"articulo/{IDIOMA}/figures"); OUT.mkdir(parents=True, exist_ok=True)
ES = IDIOMA == "es"


def t(en, es):
    """Texto en el idioma elegido."""
    return es if ES else en


def num(x, nd=2):
    """Número con nd decimales y el separador decimal del idioma."""
    s = f"{x:.{nd}f}"
    return s.replace(".", ",") if ES else s


def miles(x):
    """Entero con separador de miles del idioma (coma en inglés, punto en español)."""
    s = f"{x:,}"
    return s.replace(",", ".") if ES else s


def decimal_ejes(fig):
    if not ES:
        return
    from matplotlib.ticker import ScalarFormatter

    class Coma(ScalarFormatter):
        """ScalarFormatter con coma decimal (conserva la lógica de posiciones de matplotlib)."""
        def __call__(self, x, pos=None):
            return super().__call__(x, pos).replace(".", ",")

    for ax in fig.axes:
        for eje in (ax.xaxis, ax.yaxis):
            if type(eje.get_major_formatter()) is ScalarFormatter:
                eje.set_major_formatter(Coma())
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"          # categóricos 1–3 (validados todos los pares)
TINTA, TINTA2, REJILLA, SUPERFICIE = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
SEQ = LinearSegmentedColormap.from_list("azul", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "text.color": TINTA, "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
    "legend.frameon": False, "legend.fontsize": 7, "figure.facecolor": SUPERFICIE, "savefig.facecolor": SUPERFICIE,
})
# Costas: Natural Earth 1:110m (dominio público), guardadas en el repositorio para no depender de geopandas.
from matplotlib.collections import PolyCollection
_NE = json.load(open(Path(__file__).with_name("datos") / "naturalearth_lowres.geojson", encoding="utf-8"))
_ANILLOS = []
for _f in _NE["features"]:
    _g = _f["geometry"]
    for _poli in (_g["coordinates"] if _g["type"] == "MultiPolygon" else [_g["coordinates"]]):
        _ANILLOS.append(np.asarray(_poli[0]))      # anillo exterior (a esta escala no hay huecos relevantes)


def mundo(ax, color, edgecolor, linewidth, zorder):
    ax.add_collection(PolyCollection(_ANILLOS, facecolors=color, edgecolors=edgecolor, linewidths=linewidth, zorder=zorder))


def guardar(fig, nombre):
    decimal_ejes(fig)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(OUT / f"{nombre}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print("OK", nombre)


def costas(ax, lon=(-80, 5), lat=(25, 72)):
    mundo(ax, color="#f0efec", edgecolor="#a3a29c", linewidth=0.4, zorder=1)
    ax.set_xlim(*lon); ax.set_ylim(*lat); ax.set_aspect(1 / np.cos(np.radians(np.mean(lat))))
    ax.set_xlabel(t("Longitude [°]", "Longitud [°]")); ax.set_ylabel(t("Latitude [°]", "Latitud [°]"))


df = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                  pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
comp = pd.read_csv("datos/procesado/comparacion_tramos.csv")
val = pd.read_csv("datos/procesado/validacion_era5_adsc.csv")

# ---------------------------------------------------------------- Fig 1: cobertura
fig, ax = plt.subplots(figsize=(7, 4.2))
costas(ax)
hb = ax.hexbin(df.lon, df.lat, gridsize=70, bins="log", cmap=SEQ, mincnt=1, linewidths=0, zorder=2,
               extent=(-75, -5, 25, 75))
ax.add_patch(plt.Rectangle((-60, 30), 60, 40, fill=False, ls="--", lw=1.0, ec=TINTA, zorder=3))
ax.text(-59, 69, t("ERA5 domain (study days)", "Dominio de ERA5 (días de estudio)"), fontsize=7, color=TINTA, va="top", zorder=3)
cb = fig.colorbar(hb, ax=ax, shrink=0.8, pad=0.02); cb.set_label(t("ADS-C position reports per cell (log)", "Informes de posición ADS-C por celda (log)"))
cb.outline.set_visible(False)
ax.set_title(t(f"ADS-C position reports in the North Atlantic box (n = {miles(len(df))}), Jul 2023 – Jan 2025",
             f"Informes de posición ADS-C en la caja del Atlántico Norte (n = {miles(len(df))}), jul 2023 – ene 2025"), loc="left")
guardar(fig, "fig1_cobertura")

# ---------------------------------------------------------------- Fig 2: error de ERA5 vs viento
val["err"] = np.hypot(val.u_era - val.u_obs, val.v_era - val.v_obs)
val["vel"] = np.hypot(val.u_obs, val.v_obs)
bordes = np.arange(0, 95, 10)
val["bin"] = pd.cut(val.vel, bordes)
g = val.groupby("bin").agg(n=("err", "size"), rms=("err", lambda e: np.sqrt((e**2).mean())), p50=("err", "median"),
                           p95=("err", lambda e: e.quantile(.95)))
g = g[g.n >= 30]
x = np.array([i.mid for i in g.index])
fig, ax = plt.subplots(figsize=(3.5, 2.6))
ax.plot(x, g.rms, color=S1, marker="o", ms=4, label="RMS")
ax.plot(x, g.p50, color=S2, marker="o", ms=4, label=t("Median", "Mediana"))
ax.fill_between(x, g.p50, g.p95, color=S2, alpha=0.12, lw=0, label=t("Median – p95", "Mediana – p95"))
ax.text(x[-1] + 2, g.rms.iloc[-1], "RMS", color=TINTA, fontsize=7, va="center")
ax.text(x[-1] + 2, g.p50.iloc[-1], t("Median", "Mediana"), color=TINTA, fontsize=7, va="center")
for xi, n in zip(x, g.n):
    ax.text(xi, -0.9, f"{n}", ha="center", fontsize=5.5, color=TINTA2)
ax.set_ylim(-1.5, None)
ax.set_xlabel(t("Observed wind speed (ADS-C) [m/s]", "Viento observado (ADS-C) [m/s]")); ax.set_ylabel(t("ERA5 vector wind error [m/s]", "Error vectorial del viento de ERA5 [m/s]"))
ax.set_title(t(f"ERA5 vs in-situ wind, {miles(len(val))} reports", f"ERA5 frente al viento medido, {miles(len(val))} informes"), loc="left")
ax.legend(loc="upper left")
guardar(fig, "fig2_error_era5")

# ---------------------------------------------------------------- Fig 3: Mach real vs modelo
P = load_aircraft_engine_params()
m = df[df.mach.notna() & df.alt_ft.between(33500, 41500)].copy()
m["fl"] = (m.alt_ft / 1000).round() * 10
fig, axs = plt.subplots(1, 2, figsize=(7, 2.7), sharey=True)
for ax, tipo in zip(axs, ("B789", "A359")):
    s = m[m.typecode == tipo]
    q = s.groupby("fl").mach.quantile([.25, .5, .75]).unstack()
    q = q[s.groupby("fl").size() >= 100]
    mod = PollSchumannModel(tipo)
    for ci, col in ((0, S1), (50, S2), (100, S3)):
        econ = [optimo_local(mod, None, 0.75 * P[tipo].amass_mtow, 0.0, ci, fl_fijo=fl)[1] for fl in q.index]
        ax.plot(q.index, econ, color=col, lw=1.5, label=t(f"Model ECON, CI {ci} kg/min", f"ECON del modelo, CI {ci} kg/min"))
        ax.text(q.index[-1] + 2, econ[-1], f"CI {ci}", fontsize=6.5, color=TINTA, va="center")
    ax.errorbar(q.index, q[0.5], yerr=[q[0.5] - q[0.25], q[0.75] - q[0.5]], fmt="o", ms=4, color=TINTA,
                ecolor=TINTA2, elinewidth=1, capsize=0, label=t("Observed (median, IQR)", "Observado (mediana, RIC)"), zorder=5)
    ax.set_title(f"{tipo}: n = {miles(len(s))}", loc="left"); ax.set_xlabel(t("Flight level", "Nivel de vuelo"))
    ax.set_xlim(q.index.min() - 5, q.index.max() + 14)
axs[0].set_ylabel(t("Mach number", "Número de Mach"))
axs[0].legend(loc="lower right", fontsize=6)
guardar(fig, "fig3_mach")

# ---------------------------------------------------------------- Fig 4: validación del modelo
vm = pd.read_csv("datos/procesado/validacion_modelo_ps.csv")
fig, axs = plt.subplots(1, 2, figsize=(7, 2.6), sharey=True)
for ax, tipo in zip(axs, ("B789", "A359")):
    s = vm[vm.tipo == tipo]
    ax.axhline(0, color=TINTA2, lw=0.8)
    ax.axhspan(-5, 5, color="#f0efec", zorder=0)
    ax.plot(s.dist_km, s.PS_vs_EEA, color=S1, marker="o", ms=4, label=t("Poll–Schumann vs EEA/EUROCONTROL", "Poll–Schumann frente a EEA/EUROCONTROL"))
    ax.set_title(tipo, loc="left"); ax.set_xlabel(t("Great-circle distance [km]", "Distancia ortodrómica [km]"))
axs[0].set_ylabel(t("Trip-fuel difference [%]", "Diferencia de combustible de viaje [%]"))
axs[0].text(2900, 4.2, t("±5 %", "±5 %"), fontsize=6.5, color=TINTA2)
axs[1].legend(loc="lower right", fontsize=6)
guardar(fig, "fig4_validacion_modelo")

# ---------------------------------------------------------------- Fig 5: exceso volado vs óptimo
x = comp.exceso_vs_opt_pct.clip(-1, 3)
fig, ax = plt.subplots(figsize=(3.5, 2.6))
ax.axvspan(-0.2, 0.2, color="#f0efec", zorder=0)
ax.text(0, ax.get_ylim()[1] if False else 0, "", fontsize=1)
ax.hist(x, bins=np.arange(-1, 3.05, 0.1), color=S1, edgecolor=SUPERFICIE, linewidth=0.6, zorder=2)
yl = ax.get_ylim()[1]
for v, et in ((comp.exceso_vs_opt_pct.median(), t("median", "mediana")), (comp.exceso_vs_opt_pct.mean(), t("mean", "media"))):
    media = et in ("mean", "media")
    ax.axvline(v, color=TINTA, lw=1, ls="--" if media else "-", zorder=3)
    ax.text(v + 0.04 if media else v + 0.28, yl * (0.80 if media else 0.97), f"{et} {num(v)} %",
            fontsize=6.5, color=TINTA)
ax.text(-0.95, yl * 0.95, t("numerical\nnoise floor\n±0.2 %", "suelo de ruido\nnumérico\n±0,2 %"), fontsize=6, color=TINTA2, va="top")
ax.set_xlabel(t("Flown fuel excess over lateral optimum [%] (clipped to [−1, 3])", "Exceso de combustible de la volada sobre el óptimo lateral [%] (recortado a [−1, 3])"))
ax.set_ylabel(t("Segments", "Tramos"))
ax.set_title(t(f"{len(comp)} oceanic segments, {comp.dia.nunique()} days", f"{len(comp)} tramos oceánicos, {comp.dia.nunique()} días"), loc="left")
guardar(fig, "fig5_exceso")

# ---------------------------------------------------------------- Fig 6: ejemplo de tramo
cand = comp[(comp.tipo == "B789") | (comp.tipo == "A359")].sort_values("exceso_vs_opt_pct")
ej = cand.iloc[int(0.9 * len(cand))]                    # percentil 90: un caso con margen visible pero típico
dseg = segmentar_vuelos(df[df.icao24.notna()])
r = dseg[dseg.vuelo_id == ej.vuelo_id].sort_values("t")
campo = CampoERA5(f"datos/era5/era5_pl_{ej.dia}.nc")
t0 = float(r.t.iloc[0])
mod = PollSchumannModel(ej.tipo)
dp = optimizar_lateral(mod, campo, r.lat.iloc[0], r.lon.iloc[0], r.lat.iloc[-1], r.lon.iloc[-1], ej.fl, ej.mach,
                       0.8 * P[ej.tipo].amass_mtow, t0=t0, n_etapas=max(8, int(round(distancia(r.lat.iloc[0], r.lon.iloc[0], r.lat.iloc[-1], r.lon.iloc[-1]) / 100e3))),
                       K=20, delta_km=15.0, salto_max=3)
la0, la1 = min(r.lat.min(), dp.ruta.lat.min()) - 4, max(r.lat.max(), dp.ruta.lat.max()) + 4
lo0, lo1 = r.lon.min() - 4, r.lon.max() + 4
LA, LO = np.meshgrid(np.arange(max(la0, 30.5), min(la1, 69.5), 0.25), np.arange(max(lo0, -59.5), min(lo1, -0.5), 0.25), indexing="ij")
tmid = float(r.t.median())
u, v, _ = campo.muestra_lote(np.full(LA.size, tmid), np.full(LA.size, ej.fl), LA.ravel(), LO.ravel())
U, V = u.reshape(LA.shape), v.reshape(LA.shape)
fig, ax = plt.subplots(figsize=(7, 3.8))
pc = ax.pcolormesh(LO, LA, np.hypot(U, V), cmap=SEQ, shading="auto", zorder=0, alpha=0.7, linewidth=0)
ax.grid(False)
mundo(ax, color="none", edgecolor="#a3a29c", linewidth=0.4, zorder=1)
k = 12
ax.quiver(LO[::k, ::k], LA[::k, ::k], U[::k, ::k], V[::k, ::k], color="#ffffff", alpha=0.8, scale=900, width=0.0018, zorder=2)
ax.plot(r.lon, r.lat, color=S2, marker="o", ms=3.5, lw=1.5, label=t("Flown (ADS-C reports)", "Volada (informes ADS-C)"), zorder=4)
_gc = Ruta([r.lat.iloc[0], r.lat.iloc[-1]], [r.lon.iloc[0], r.lon.iloc[-1]])
_pts = np.array([_gc.punto(f * _gc.longitud)[:2] for f in np.linspace(0, 1, 40)])
ax.plot(_pts[:, 1], _pts[:, 0], color="#ffffff", lw=1.0, ls=(0, (3, 2)), label="Great circle (drawn on top)", zorder=6)
ax.plot(dp.ruta.lon, dp.ruta.lat, color="#ffffff", lw=3.2, zorder=4)
ax.plot(dp.ruta.lon, dp.ruta.lat, color=TINTA, lw=1.5, label="DP lateral optimum", zorder=5)
ax.set_xlim(LO.min(), LO.max()); ax.set_ylim(LA.min(), LA.max())
ax.set_aspect(1 / np.cos(np.radians(LA.mean())))
cb = fig.colorbar(pc, ax=ax, shrink=0.85, pad=0.02); cb.set_label(t(f"ERA5 wind speed at FL{ej.fl:.0f} [m/s]", f"Viento de ERA5 en FL{ej.fl:.0f} [m/s]")); cb.outline.set_visible(False)
ax.set_xlabel(t("Longitude [°]", "Longitud [°]")); ax.set_ylabel(t("Latitude [°]", "Latitud [°]"))
ax.set_title(t(f"Example segment ({ej.tipo}, {ej.dia}, FL{ej.fl:.0f}, M{ej.mach:.3f})\nflown excess {num(ej.exceso_vs_opt_pct)} % "
               f"({ej.comb_volada - ej.comb_opt:.0f} kg)",
               f"Tramo de ejemplo ({ej.tipo}, {ej.dia}, FL{ej.fl:.0f}, M{num(ej.mach, 3)})\nexceso de la volada {num(ej.exceso_vs_opt_pct)} % "
               f"({ej.comb_volada - ej.comb_opt:.0f} kg)"), loc="left")
from matplotlib.lines import Line2D
manejadores = [Line2D([], [], color=S2, marker="o", ms=3.5, lw=1.5, label=t("Flown (ADS-C reports)", "Volada (informes ADS-C)")),
               Line2D([], [], color=TINTA, lw=1.5, label=t("DP lateral optimum (unconstrained)", "Óptimo lateral DP (sin restricción)")),
               Line2D([], [], color=TINTA2, lw=1.0, ls=(0, (3, 2)), label=t("Great circle (white dashes, on top)", "Ortodrómica (trazos blancos, encima)"))]
_est = Path("datos/procesado/estructura_tramos.csv")
if _est.exists():
    e = pd.read_csv(_est)
    e = e[e.vuelo_id == ej.vuelo_id]
    if len(e) and isinstance(e.lat_rest.iloc[0], str):
        lats = [float(v) for v in e.lat_rest.iloc[0].split(";")]
        sgn = 1 if r.lon.iloc[-1] > r.lon.iloc[0] else -1
        lo, hi = sorted((r.lon.iloc[0], r.lon.iloc[-1]))
        mer = sorted([mm for mm in range(-60, 1, 10) if lo + 0.25 < mm < hi - 0.25], key=lambda mm: sgn * mm)
        rr = Ruta([r.lat.iloc[0], *lats, r.lat.iloc[-1]], [r.lon.iloc[0], *mer, r.lon.iloc[-1]])
        pts = np.array([rr.punto(f * rr.longitud)[:2] for f in np.linspace(0, 1, 60)])
        ax.plot(pts[:, 1], pts[:, 0], color=S3, lw=1.5, ls=(0, (5, 2)), zorder=5)
        ax.plot(mer, lats, ls="none", marker="D", ms=4, color=S3, zorder=6)
        manejadores.insert(2, Line2D([], [], color=S3, lw=1.5, ls=(0, (5, 2)), marker="D", ms=4,
                                     label=t("Optimum through 0.5° waypoints at 10° meridians",
                                             "Óptimo por puntos de 0,5° en meridianos de 10°")))
ax.legend(handles=manejadores, loc="lower left", fontsize=6.5, frameon=True, facecolor=SUPERFICIE, edgecolor="none")
guardar(fig, "fig6_ejemplo_tramo")
