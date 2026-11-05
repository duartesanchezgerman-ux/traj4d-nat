"""Fase 1 — análisis físico base del B789 (crucero). Ejecutar:  PYTHONPATH=src python scripts/f1_analisis.py

Hipótesis declaradas (se sustituirán por datos reales en las fases 2–3):
  * Distancia de crucero 5.770 km (ortodrómica MAD–JFK; el tramo de crucero real es menor).
  * Masa al inicio del crucero 225 t (supuesta: 787-9 con unas 250 plazas y combustible de ruta + reservas).
  * Viento uniforme a lo largo de la ruta: 0 / -30 m/s (oeste, en contra) / +30 m/s (este, a favor).
  * ΔISA = 0. Coste por tiempo nulo (mínimo combustible) salvo indicación.
"""
import warnings

import numpy as np

from traj4d import atmosfera as atm
from traj4d.crucero import (Envolvente, integrar_crucero, optimo_local, perfil_constante,
                            perfil_escalones, perfil_optimo, rango_breguet)
from traj4d.rendimiento import OpenAPModel, PollSchumannModel

warnings.filterwarnings("ignore")
ps, oa, env = PollSchumannModel(), OpenAPModel(), Envolvente()
R, M0 = 5_770e3, 225e3

print("=" * 78)
print("A. Envolvente: FL máximo válido a M0,85 según la masa (Poll–Schumann)")
fls = np.arange(290, 432, 2)
for m in (250e3, 230e3, 210e3, 190e3, 170e3):
    ok = env.valido(m, 0.85, fls)
    print(f"  {m/1e3:5.0f} t -> FL máx {int(fls[ok].max()) if ok.any() else '—'}")

print("=" * 78)
print("B. Óptimo local (FL, M) que maximiza el SR sobre el suelo, según masa y viento")
print("   masa   |   w=-30 m/s        |   w=0              |   w=+30 m/s")
for m in (230e3, 210e3, 190e3, 170e3):
    row = []
    for w in (-30.0, 0.0, 30.0):
        fl, M, J = optimo_local(ps, env, m, viento=w)
        sr = 1.0 / J
        row.append(f"FL{fl:.0f} M{M:.3f} {sr:5.1f}m/kg")
    print(f"  {m/1e3:4.0f} t | " + " | ".join(row))
print("  -> con viento en contra el Mach óptimo SUBE y con viento a favor BAJA: es la "
      "variación de velocidad de la que sale el 4,2 % de Wells et al.")

print("=" * 78)
print(f"C. Crucero de {R/1e3:.0f} km desde {M0/1e3:.0f} t — perfiles verticales y de velocidad")
for w, etiqueta in ((0.0, "sin viento"), (-30.0, "en contra 30 m/s (oeste)"), (30.0, "a favor 30 m/s (este)")):
    casos = {
        "FL350 const, M0,85": perfil_constante(350, 0.85),
        "FL370 const, M0,85": perfil_constante(370, 0.85),
        "Escalones 2000 ft desde FL350, M0,85": perfil_escalones(ps, env, 350, 0.85),
        "Ascenso continuo óptimo, M0,85 fijo": perfil_optimo(ps, env, mach_fijo=0.85),
        "Ascenso continuo + Mach óptimos": perfil_optimo(ps, env),
    }
    print(f"  [{etiqueta}]")
    base = None
    for nombre, perfil in casos.items():
        r = integrar_crucero(ps, M0, R, perfil, viento=w, ds=20e3)
        if base is None and nombre.startswith("Escalones"):
            pass
        res = (r.combustible, r.tiempo_h, r.fl[0], r.fl[-1], r.mach.mean())
        casos[nombre] = res
    base = casos["Escalones 2000 ft desde FL350, M0,85"][0]
    for nombre, (f, t, fli, flf, mm) in casos.items():
        print(f"    {nombre:40s} comb {f/1e3:6.2f} t  ({100*(f/base-1):+5.2f} % vs escalones)  "
              f"{t:4.2f} h  FL{fli:.0f}->FL{flf:.0f}  M̄ {mm:.3f}")

print("=" * 78)
print("D. Error de Breguet (L/D y TSFC del punto inicial) frente a la integración, FL370 M0,85, sin viento")
r = integrar_crucero(ps, M0, R, perfil_constante(370, 0.85), ds=10e3)
e0 = ps.estado(M0, 0.85, 370)
V = float(e0["TAS"])
R_b = rango_breguet(M0, r.m[-1], V, float(e0["LD"]), float(e0["TSFC"]))
em = ps.estado(0.5 * (M0 + r.m[-1]), 0.85, 370)
R_bm = rango_breguet(M0, r.m[-1], V, float(em["LD"]), float(em["TSFC"]))
print(f"  Combustible integrado {r.combustible/1e3:.2f} t para {R/1e3:.0f} km")
print(f"  Breguet con el punto inicial: {R_b/1e3:.0f} km ({100*(R_b/R-1):+.1f} %)  |  "
      f"con la masa media: {R_bm/1e3:.0f} km ({100*(R_bm/R-1):+.1f} %)")

print("=" * 78)
print("E. Sesgo entre modelos: mismo crucero (FL370 const, M0,85, sin viento)")
for mod in (ps, oa):
    r = integrar_crucero(mod, M0, R, perfil_constante(370, 0.85), ds=20e3)
    print(f"  {mod.nombre:22s} combustible {r.combustible/1e3:6.2f} t")
