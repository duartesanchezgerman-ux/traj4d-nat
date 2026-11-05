"""Fase 1b: coste del tiempo (cost index) y física lateral con un chorro sintético.
Ejecutar:  PYTHONPATH=src python scripts/f1b_analisis.py

Hipótesis declaradas: masa al inicio del crucero 225 t; ΔISA = 0; chorro sintético (no es un dato
meteorológico: sirve para aislar la física antes de usar ERA5). CI en kg/min (convención Airbus).
"""
import warnings

import numpy as np

from traj4d.crucero import (Envolvente, integrar_crucero, integrar_ruta, optimo_local,
                            perfil_escalones, perfil_mach_optimo_fl_escalones, perfil_optimo)
from traj4d.rendimiento import PollSchumannModel
from traj4d.ruta import AEROPUERTOS, Ruta
from traj4d.viento import ChorroSintetico

warnings.filterwarnings("ignore")
ps, env = PollSchumannModel(), Envolvente()
R, M0 = 5_770e3, 225e3

print("=" * 86)
print("A. Mach ECON (FL libre, sin viento) según el cost index [kg/min]")
cis = [0, 10, 20, 30, 40, 60, 80, 100, 150]
print("   masa  | " + " ".join(f"CI{c:>4d}" for c in cis))
for m in (230e3, 200e3, 170e3):
    print(f"  {m/1e3:4.0f} t | " + " ".join(f"{optimo_local(ps, env, m, 0.0, c)[1]:6.3f}" for c in cis))
ci_impl = next(c for c in range(0, 400, 5) if optimo_local(ps, env, 200e3, 0.0, c)[1] >= 0.85)
print(f"  -> Con este modelo, volar a M0,85 a 200 t equivale a un CI ≈ {ci_impl} kg/min.")

print("=" * 86)
print(f"B. Crucero de {R/1e3:.0f} km desde {M0/1e3:.0f} t. Comparación A COSTE TOTAL "
      "(combustible + CI·min), referencia = escalones + M0,85")
for ci in (0, ci_impl // 2, ci_impl):
    print(f"  --- CI = {ci} kg/min ---")
    for w, et in ((-30.0, "contra 30"), (0.0, "sin viento"), (30.0, "favor 30")):
        casos = {
            "Escalones + M0,85 (ref.)": perfil_escalones(ps, env, 350, 0.85, ci=ci),
            "Escalones + Mach ECON": perfil_mach_optimo_fl_escalones(ps, env, 350, ci=ci),
            "Continuo FL+Mach ECON": perfil_optimo(ps, env, ci=ci),
        }
        res = {k: integrar_crucero(ps, M0, R, p, viento=w, ds=20e3) for k, p in casos.items()}
        ref = res["Escalones + M0,85 (ref.)"]
        linea = []
        for k, r in res.items():
            if k.startswith("Escalones + M0,85"):
                continue
            dc = 100 * (r.combustible / ref.combustible - 1)
            dk = 100 * (r.coste(ci) / ref.coste(ci) - 1)
            linea.append(f"{k}: comb {dc:+5.2f} %  coste {dk:+5.2f} %  t {60*(r.tiempo_h-ref.tiempo_h):+5.1f} min")
        print(f"   [{et:10s}] " + " | ".join(linea))

print("=" * 86)
print("C. Física lateral: MAD<->JFK atravesando un chorro sintético (75 m/s, 45°N, FL340), CI =", ci_impl)
jet = ChorroSintetico()
gc = Ruta.ortodromica("LEMD", "KJFK", n=40)
print(f"  Ortodrómica: {gc.longitud/1e3:.0f} km, latitud máxima {gc.lat.max():.1f}°N")
(la_m, lo_m), (la_j, lo_j) = AEROPUERTOS["LEMD"], AEROPUERTOS["KJFK"]
alternativas = {
    "Ortodrómica": None,
    "Por 37N 40W (sur del chorro)": (37.0, -40.0),
    "Por 45N 40W (núcleo del chorro)": (45.0, -40.0),
    "Por 53N 40W (norte del chorro)": (53.0, -40.0),
}
for sentido, (o, d) in (("Oeste MAD->JFK", ("LEMD", "KJFK")), ("Este JFK->MAD", ("KJFK", "LEMD"))):
    print(f"  [{sentido}]")
    base = None
    for nombre, wp in alternativas.items():
        if wp is None:
            ruta = Ruta.ortodromica(o, d, n=40)
        else:
            (a1, b1), (a2, b2) = AEROPUERTOS[o], AEROPUERTOS[d]
            ruta = Ruta([a1, wp[0], a2], [b1, wp[1], b2])
        r = integrar_ruta(ps, M0, ruta, perfil_mach_optimo_fl_escalones(ps, env, 350, ci=ci_impl), jet, ds=20e3)
        base = base or r
        print(f"    {nombre:34s} dist {ruta.longitud/1e3:6.0f} km  comb {r.combustible/1e3:6.2f} t "
              f"({100*(r.combustible/base.combustible-1):+5.2f} %)  t {r.tiempo_h:4.2f} h  "
              f"coste {100*(r.coste(ci_impl)/base.coste(ci_impl)-1):+5.2f} %")
