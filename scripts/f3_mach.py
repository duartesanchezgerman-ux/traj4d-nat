"""Fase 3: Mach real (ADS-C, Tag 15) frente al Mach óptimo del modelo Poll–Schumann.
Ejecutar:  PYTHONPATH=src python scripts/f3_mach.py | tee resultados/resultados_f3_mach.txt

El Mach que vuela un avión depende del cost index (CI) de la aerolínea, de la masa y del viento.
Ninguno de los tres se observa directamente. Por eso no se pide al modelo que "acierte" un Mach:
se calcula, para cada nivel de vuelo observado, el Mach ECON del modelo con CI = 0, 25, 50, 75 y
100 kg/min y una masa de referencia de 0,75·MTOW (con sensibilidad a 0,65 y 0,85), sin viento.
Después se busca el CI que reproduce la mediana observada ("CI implícito").
Interpretación: si el CI implícito es plausible (decenas de kg/min), el modelo es coherente con la
operación real. Si hiciera falta un CI enorme, el Mach de diseño del modelo estaría sesgado a la baja.
"""
import warnings

import numpy as np
import pandas as pd
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params

from traj4d.crucero import Envolvente, optimo_local
from traj4d.rendimiento import PollSchumannModel
from traj4d.vuelos import asignar_tipo, cargar_informes

warnings.filterwarnings("ignore")
P = load_aircraft_engine_params()
df = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                  pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
df = df[df.mach.notna() & df.alt_ft.between(30000, 43000)]
df["fl"] = (df.alt_ft / 1000).round() * 10          # redondeo a 1000 ft
CIS = [0, 25, 50, 75, 100]

for tipo in ("B789", "A359"):
    s = df[df.typecode == tipo]
    mod = PollSchumannModel(tipo)
    env = Envolvente(tipo)
    mtow = P[tipo].amass_mtow
    print("=" * 96)
    print(f"{tipo}: {len(s)} informes con Mach | mediana {s.mach.median():.3f} | IQR {s.mach.quantile(.25):.3f}–{s.mach.quantile(.75):.3f}"
          f" | Mach de diseño Poll–Schumann {P[tipo].m_des:.3f}")
    print("   FL  |   n  | Mach obs (p25/p50/p75) | Mach ECON del modelo a 0,75·MTOW con CI = " + "/".join(map(str, CIS)) + " | CI implícito")
    for fl, g in s.groupby("fl"):
        if len(g) < 100:
            continue
        econ = [optimo_local(mod, None, 0.75 * mtow, 0.0, ci, fl_fijo=fl)[1] for ci in CIS]
        obs = g.mach.median()
        ci_imp = next((c for c in range(0, 301, 5) if optimo_local(mod, None, 0.75 * mtow, 0.0, c, fl_fijo=fl)[1] >= obs), np.inf)
        print(f"  {fl:4.0f} | {len(g):5d} | {g.mach.quantile(.25):.3f} / {obs:.3f} / {g.mach.quantile(.75):.3f}    | "
              + " ".join(f"{m:.3f}" for m in econ) + f"   | {ci_imp:>4} kg/min")
    obs = s.mach.median()
    fl_med = float(s.fl.median())
    print(f"  Sensibilidad a la masa (FL{fl_med:.0f}, mediana observada {obs:.3f}): CI implícito con", end=" ")
    for frac in (0.65, 0.75, 0.85):
        ci_imp = next((c for c in range(0, 301, 5) if optimo_local(mod, None, frac * mtow, 0.0, c, fl_fijo=fl_med)[1] >= obs), np.inf)
        print(f"{frac:.2f}·MTOW -> {ci_imp} kg/min;", end=" ")
    print()
print("\nNota: sin viento. Con viento en contra el Mach ECON sube y con viento a favor baja (Fase 1b), así que el CI")
print("implícito tiene una incertidumbre de unos ±10–20 kg/min según el régimen de vientos del día.")
