"""Fase 4: reconstrucción de cruces completos del Atlántico Norte uniendo ADS-B (adsb.lol) y ADS-C (OpenSky).

Para cada tramo ADS-C (un avión, informes oceánicos consecutivos) se buscan en las trazas ADS-B del mismo avión:
  * el último punto ADS-B antes del primer informe ADS-C (lado de entrada al océano) y el despegue de ese vuelo;
  * el primer punto ADS-B después del último informe ADS-C (lado de salida) y el aterrizaje.
Entre el último ADS-B de entrada y el primer ADS-C (y entre el último ADS-C y el primer ADS-B de salida) quedan
huecos sin posiciones. Se miden su duración y su distancia; la reconstrucción del recorrido dentro de los huecos
se valida aparte (cierre en tiempo con el viento de ERA5 y el Mach medido).

Formato de adsb.lol (readsb, trace_full): JSON comprimido con gzip, con "timestamp" (inicio del día, s) y "trace",
una lista de puntos [dt, lat, lon, alt, gs, track, flags, vrate, aircraft|None, fuente, alt_geom, geom_rate, ias, roll].
alt es un entero en pies o la cadena "ground". Si "aircraft" no es None trae, cuando el avión lo emite por Mode S
EHS, mach, ws (viento, kt), wd (dirección de procedencia), oat y tat.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .ruta import distancia

HUECO_MAX_H = 14.0          # un vuelo transatlántico no dura más; separa vuelos distintos del mismo avión
HUECO_TIERRA_S = 90 * 60    # un hueco mayor en la traza separa vuelos distintos (para buscar el despegue/aterrizaje)
ALT_TIERRA_FT = 3000.0      # primer/último punto de la traza por debajo de esta altitud: despegue/aterrizaje visto


def leer_traza(ruta) -> pd.DataFrame:
    """Lee un trace_full de adsb.lol. Devuelve t (s UNIX), lat, lon, alt_ft (NaN en tierra), en_tierra, gs_kt,
    track_deg y, si existen, mach, ws_kt, wd_deg y oat_c."""
    b = Path(ruta).read_bytes()
    d = json.loads(gzip.decompress(b) if b[:2] == b"\x1f\x8b" else b)
    t0 = float(d["timestamp"])
    filas = []
    for p in d["trace"]:
        ac = p[8] if len(p) > 8 and isinstance(p[8], dict) else {}
        alt = p[3]
        filas.append((t0 + p[0], p[1], p[2], np.nan if alt == "ground" or alt is None else float(alt), alt == "ground",
                      p[4], p[5], ac.get("mach"), ac.get("ws"), ac.get("wd"), ac.get("oat")))
    df = pd.DataFrame(filas, columns=["t", "lat", "lon", "alt_ft", "en_tierra", "gs_kt", "track_deg",
                                      "mach", "ws_kt", "wd_deg", "oat_c"])
    for c in ("gs_kt", "track_deg", "mach", "ws_kt", "wd_deg", "oat_c"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df.attrs.update(icao24=d.get("icao"), matricula=d.get("r"), tipo=d.get("t"))
    return df


def trazas_de(icao24: str, dias: list[str], base="datos/adsblol") -> pd.DataFrame:
    """Une las trazas de un avión en varios archivos-día (sin duplicados, ordenadas por tiempo)."""
    partes = []
    for d in dias:
        f = Path(base) / d / f"{icao24.lower()}.json.gz"
        if f.exists():
            partes.append(leer_traza(f))
    if not partes:
        return pd.DataFrame(columns=["t", "lat", "lon", "alt_ft", "en_tierra"])
    return pd.concat(partes, ignore_index=True).drop_duplicates("t").sort_values("t").reset_index(drop=True)


def _extremo_en_tierra(adsb: pd.DataFrame, t_ref: float, antes: bool):
    """Despegue (antes=True) o aterrizaje del vuelo, si la traza ADS-B llega hasta cerca del suelo.

    Se toma el tramo CONTINUO de la traza a ese lado del océano (se corta en huecos de más de HUECO_TIERRA_S, que
    separan vuelos distintos del mismo avión) y, dentro de él, el último punto en tierra o por debajo de
    ALT_TIERRA_FT antes del océano (despegue) o el primero después (aterrizaje). Antes (hasta el 3-11-2026) se exigía un punto marcado "ground", que solo
    emiten algunos aeropuertos con cobertura de superficie: se veían el despegue y el aterrizaje en un 6 % de los
    tramos; con este criterio, en torno a la mitad."""
    if antes:
        sub = adsb[(adsb.t < t_ref) & (adsb.t > t_ref - HUECO_MAX_H * 3600)].sort_values("t")
    else:
        sub = adsb[(adsb.t > t_ref) & (adsb.t < t_ref + HUECO_MAX_H * 3600)].sort_values("t")
    if sub.empty:
        return None
    cortes = np.where(np.diff(sub.t.to_numpy()) > HUECO_TIERRA_S)[0]
    if len(cortes):
        sub = sub.iloc[cortes[-1] + 1:] if antes else sub.iloc[:cortes[0] + 1]
    es_bajo = (sub.en_tierra | (sub.alt_ft <= ALT_TIERRA_FT)).to_numpy()
    if not es_bajo.any():
        return None
    # La última racha de puntos bajos antes del océano es la salida de ESTE vuelo; la primera después, su llegada
    # (así no se toma el vuelo siguiente si el avión vuelve a salir sin que se corte la traza).
    k = np.flatnonzero(es_bajo)
    if antes:
        ini = k[-1]
        while ini > 0 and es_bajo[ini - 1]:
            ini -= 1
        racha = sub.iloc[ini:k[-1] + 1]
    else:
        fin = k[0]
        while fin < len(es_bajo) - 1 and es_bajo[fin + 1]:
            fin += 1
        racha = sub.iloc[k[0]:fin + 1]
    tierra = racha[racha.en_tierra]
    if len(tierra):                                   # despegue: último punto en tierra; aterrizaje: primero
        return tierra.iloc[-1] if antes else tierra.iloc[0]
    alt = racha.alt_ft.to_numpy()                     # sin puntos en tierra: el punto más bajo de la racha
    if np.all(np.isnan(alt)):
        return racha.iloc[-1] if antes else racha.iloc[0]
    i = len(alt) - 1 - int(np.nanargmin(alt[::-1])) if antes else int(np.nanargmin(alt))
    return racha.iloc[i]


def unir_cruce(adsc: pd.DataFrame, adsb: pd.DataFrame) -> dict:
    """Une un tramo ADS-C (columnas t, lat, lon, alt_ft) con la traza ADS-B del mismo avión.

    Devuelve un diccionario con los puntos de enlace, la duración (h) y la distancia (km) de los dos huecos, y el
    despegue y el aterrizaje si se ven. Los valores ausentes son NaN/None."""
    a = adsc.sort_values("t")
    t_ini, t_fin = float(a.t.iloc[0]), float(a.t.iloc[-1])
    r = dict(t_adsc_ini=t_ini, t_adsc_fin=t_fin, lat_adsc_ini=a.lat.iloc[0], lon_adsc_ini=a.lon.iloc[0],
             lat_adsc_fin=a.lat.iloc[-1], lon_adsc_fin=a.lon.iloc[-1])
    vuelo = adsb[~adsb.en_tierra] if len(adsb) else adsb
    antes = vuelo[(vuelo.t < t_ini) & (vuelo.t > t_ini - HUECO_MAX_H * 3600)] if len(vuelo) else vuelo
    despues = vuelo[(vuelo.t > t_fin) & (vuelo.t < t_fin + HUECO_MAX_H * 3600)] if len(vuelo) else vuelo
    for lado, sub, fila_fn, ref in (("entrada", antes, lambda s: s.iloc[-1], (a.lat.iloc[0], a.lon.iloc[0], t_ini)),
                                    ("salida", despues, lambda s: s.iloc[0], (a.lat.iloc[-1], a.lon.iloc[-1], t_fin))):
        if len(sub):
            p = fila_fn(sub)
            r.update({f"t_adsb_{lado}": float(p.t), f"lat_adsb_{lado}": float(p.lat), f"lon_adsb_{lado}": float(p.lon),
                      f"alt_adsb_{lado}": float(p.alt_ft),
                      f"hueco_{lado}_h": abs(ref[2] - float(p.t)) / 3600,
                      f"hueco_{lado}_km": float(distancia(ref[0], ref[1], float(p.lat), float(p.lon))) / 1e3})
        else:
            r.update({f"t_adsb_{lado}": np.nan, f"hueco_{lado}_h": np.nan, f"hueco_{lado}_km": np.nan})
    dep = _extremo_en_tierra(adsb, t_ini, True) if len(adsb) else None
    arr = _extremo_en_tierra(adsb, t_fin, False) if len(adsb) else None
    r.update(alt_dep=float(0.0 if dep is not None and dep.en_tierra else dep.alt_ft) if dep is not None else np.nan,
             alt_arr=float(0.0 if arr is not None and arr.en_tierra else arr.alt_ft) if arr is not None else np.nan,
             t_despegue=float(dep.t) if dep is not None else np.nan,
             lat_dep=float(dep.lat) if dep is not None else np.nan, lon_dep=float(dep.lon) if dep is not None else np.nan,
             t_aterrizaje=float(arr.t) if arr is not None else np.nan,
             lat_arr=float(arr.lat) if arr is not None else np.nan, lon_arr=float(arr.lon) if arr is not None else np.nan)
    return r
