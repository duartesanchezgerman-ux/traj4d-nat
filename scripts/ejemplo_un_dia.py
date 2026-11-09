"""Quick reproducibility check on one study day (called by scripts/ejemplo_un_dia.sh).

Compares the one-day results with the values of the full 70-day run (docs/ejemplo_un_dia_referencia.json)
and draws a one-day version of Figures 2 and 5 of the article.

    PYTHONPATH=src python scripts/ejemplo_un_dia.py 2023-07-09 resultados/ejemplo_2023-07-09
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

dia, out = sys.argv[1], Path(sys.argv[2])
ref = json.load(open("docs/ejemplo_un_dia_referencia.json", encoding="utf-8"))["dias"].get(dia)
if ref is None:
    sys.exit(f"{dia} is not one of the 70 study days (docs/dias_estudio_era5.txt).")

v = pd.read_csv(out / "validacion_era5_adsc.csv")
v = v[pd.to_datetime(v.t, unit="s").dt.strftime("%Y-%m-%d") == dia]
t = pd.read_csv(out / "comparacion_tramos.csv")
e = np.hypot(v.u_era - v.u_obs, v.v_era - v.v_obs)
obt = {
    "informes_viento": len(v), "error_vectorial_rms_ms": float(np.sqrt((e**2).mean())),
    "tramos": len(t), "exceso_medio_pct": float(t.exceso_vs_opt_pct.mean()),
    "exceso_mediana_pct": float(t.exceso_vs_opt_pct.median()),
    "ganancia_opt_vs_gc_media_pct": float(t.ganancia_opt_vs_gc_pct.mean()),
}
etiquetas = {
    "informes_viento": "valid wind reports", "error_vectorial_rms_ms": "ERA5 vector RMS error [m/s]",
    "tramos": "oceanic segments", "exceso_medio_pct": "mean fuel excess, flown vs optimal [%]",
    "exceso_mediana_pct": "median fuel excess [%]", "ganancia_opt_vs_gc_media_pct": "mean gain, optimal vs great circle [%]",
}
ok = True
print(f"\nOne-day check, {dia}: this run vs the full 70-day run")
for k, r in ref.items():
    igual = obt[k] == r if isinstance(r, int) else abs(obt[k] - r) <= 0.0015
    ok &= igual
    print(f"  {etiquetas[k]:<42} {obt[k]:>9.3f}  reference {r:>9.3f}  {'OK' if igual else 'DIFFERENT'}")

fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
vo, ve = np.hypot(v.u_obs, v.v_obs), np.hypot(v.u_era, v.v_era)
ax[0].scatter(vo, ve, s=4, alpha=.5)
m = float(max(vo.max(), ve.max())) * 1.05
ax[0].plot([0, m], [0, m], lw=.8, color="k")
ax[0].set(xlabel="Measured wind speed [m/s]", ylabel="ERA5 wind speed [m/s]",
          title=f"ERA5 vs on-board wind ({len(v)} reports)")
ax[1].hist(t.exceso_vs_opt_pct, bins=np.arange(-0.5, 2.05, 0.1))
ax[1].axvline(t.exceso_vs_opt_pct.mean(), color="k", lw=.8)
ax[1].set(xlabel="Fuel excess, flown vs wind-optimal [%]", ylabel="Segments",
          title=f"Flown vs optimal ({len(t)} segments)")
fig.suptitle(f"traj4d one-day example, {dia}", fontsize=10)
fig.tight_layout()
fig.savefig(out / f"ejemplo_{dia}.png", dpi=150)
print(f"\nFigure: {out / f'ejemplo_{dia}.png'}")
print("RESULT:", "all values match the full run." if ok else "some values differ; see above.")
sys.exit(0 if ok else 1)
