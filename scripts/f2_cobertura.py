"""Fase 2: calidad del ADS-C y cobertura del Atlántico Norte.
Ejecutar:  PYTHONPATH=src python scripts/f2_cobertura.py | tee resultados/resultados_f2_cobertura.txt

Entradas: datos/procesado/adsc_atlantico.csv.gz y datos/procesado/aeronaves_adsc.csv
Salida adicional: datos/procesado/segmentos_nat.csv (un registro por tramo observado)
"""
import warnings

import pandas as pd
from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params

from traj4d.vuelos import asignar_tipo, cargar_informes, cierre_triangulo, resumen_vuelos, segmentar_vuelos

warnings.filterwarnings("ignore")
PS = set(load_aircraft_engine_params().keys())

df = asignar_tipo(cargar_informes("datos/procesado/adsc_atlantico.csv.gz"),
                  pd.read_csv("datos/procesado/aeronaves_adsc.csv", dtype=str))
print("=" * 80)
print(f"Informes en la caja NAT: {len(df)} | con icao24: {df.icao24.notna().sum()} "
      f"(recuperados por matrícula: {int(df.icao24_por_matricula.sum())}) | con tipo: {df.typecode.notna().sum()}")
print(f"Con viento medido (Tag 16): {df.viento_kt.notna().sum()} | con Mach (Tag 15): {df.mach.notna().sum()}")

print("=" * 80)
c = cierre_triangulo(df)
q = c.residuo_kt.quantile([.5, .9, .95, .99])
print(f"Cierre del triángulo de viento (informes con 14+15+16): n = {len(c)}")
print(f"  residuo vectorial [kt]: p50 {q[.5]:.2f} | p90 {q[.9]:.2f} | p95 {q[.95]:.2f} | p99 {q[.99]:.2f} | "
      f">10 kt: {100 * (c.residuo_kt > 10).mean():.1f} %")
for tc in ("B789", "A359"):
    s = df[(df.typecode == tc) & df.mach.notna()]
    print(f"  Mach real {tc}: mediana {s.mach.median():.3f} | p10 {s.mach.quantile(.1):.3f} | p90 {s.mach.quantile(.9):.3f} | n = {len(s)}")

print("=" * 80)
df = segmentar_vuelos(df[df.icao24.notna()])
v = resumen_vuelos(df)
v["modelable"] = v.typecode.isin(PS)
seg = v[(v.span_lon >= 15) & (v.n >= 4) & v.modelable]
print(f"Segmentos (hueco > 3 h = vuelo nuevo): {len(v)}")
for s_ in (15, 20, 25, 30):
    print(f"  tramos con extensión >= {s_}° de longitud y >= 4 informes (tipo modelable): "
          f"{((v.span_lon >= s_) & (v.n >= 4) & v.modelable).sum()}")
print(f"Unidad de análisis (>= 15°): {len(seg)} tramos | B789 {(seg.typecode == 'B789').sum()} | "
      f"A359 {(seg.typecode == 'A359').sum()} | intervalo mediano {seg.intervalo_mediano_min.median():.1f} min | "
      f"longitud mediana de los tramos: {seg.lon_min.median():.1f}° a {seg.lon_max.median():.1f}°")
print(f"Longitud de los informes (p10/p50/p90): {df.lon.quantile(.1):.1f} / {df.lon.median():.1f} / {df.lon.quantile(.9):.1f}")

print("=" * 80)
met = df[df.viento_kt.notna()].groupby(df.fecha.dt.date).size()
tab = pd.DataFrame({"tramos": seg.groupby("fecha").size(),
                    "b789_a359": seg[seg.typecode.isin(["B789", "A359"])].groupby("fecha").size(),
                    "vientos_medidos": met}).fillna(0).astype(int)
tab["mes"] = pd.to_datetime(tab.index.astype(str)).month
print("Días de estudio propuestos (más tramos):")
print("  INVIERNO"); print(tab[tab.mes.isin([11, 12, 1, 2])].sort_values("tramos", ascending=False).head(6).to_string())
print("  VERANO"); print(tab[tab.mes.isin([6, 7, 8])].sort_values("tramos", ascending=False).head(4).to_string())
seg.to_csv("datos/procesado/segmentos_nat.csv")
