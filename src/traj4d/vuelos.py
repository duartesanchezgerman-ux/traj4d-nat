"""Tablas de informes ADS-C → vuelos: tipo de aeronave, segmentación, calidad del triángulo de viento."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import atmosfera as atm
from .ruta import distancia

KT = atm.KT


def cargar_informes(ruta_csv_gz, ruta_todas=None) -> pd.DataFrame:
    """Lee una tabla de informes ADS-C y corrige la vuelta de hora (ver corregir_vuelta_de_hora).

    La corrección se decide con TODOS los informes de posición (adsc_posiciones.csv.gz, en la misma carpeta si
    existe, o `ruta_todas`), para que los vecinos de un informe no dependan del recorte geográfico de la tabla."""
    df = pd.read_csv(ruta_csv_gz, dtype={"icao24": "string", "icao24_bruto": "string", "matricula": "string",
                                         "vuelo": "string", "tipo_informe": "string", "precision": "string"})
    todas = Path(ruta_todas) if ruta_todas else Path(ruta_csv_gz).with_name("adsc_posiciones.csv.gz")
    if todas.exists() and todas.resolve() != Path(ruta_csv_gz).resolve():
        ref = corregir_vuelta_de_hora(pd.read_csv(todas, usecols=["t", "icao24", "matricula", "lat", "lon"],
                                                  dtype={"icao24": "string", "matricula": "string"}))
        ref = ref[ref.t_corregida_1h]
        malas = set(zip(ref["icao24"].astype(object).where(ref["icao24"].notna(), ref["matricula"].astype(object)),
                        ref["t"] + 3600))
        clave = df["icao24"].astype(object).where(df["icao24"].notna(), df["matricula"].astype(object))
        k = pd.Series([(c, t) in malas for c, t in zip(clave, df["t"])], index=df.index)
        df["t_corregida_1h"] = k
        df.loc[k, "t"] -= 3600
        if "wp1_eta_t" in df:
            df.loc[k, "wp1_eta_t"] -= 3600
    else:
        df = corregir_vuelta_de_hora(df)
    df["fecha"] = pd.to_datetime(df["t"], unit="s", utc=True)
    return df.sort_values(["t"]).reset_index(drop=True)


def corregir_vuelta_de_hora(df: pd.DataFrame, err_max_km=150.0, err_ok_km=40.0, min_minuto=50) -> pd.DataFrame:
    """Corrige los informes que el volcado fecha una hora tarde al cambiar la hora.

    El informe básico ADS-C solo lleva los segundos desde el comienzo de la hora; la fecha y la hora
    las pone el receptor. Si el informe se genera a hh:5x y llega después de las (hh+1):00, el volcado
    lo fecha a (hh+1):5x. Se detecta así (comprobado el 3-11-2026, 0,15 % de los informes, casi todos en
    los minutos 57–59): con los demás informes del mismo avión, la posición en t queda a más de
    `err_max_km` del informe y en t − 3600 s a menos de `err_ok_km`. La posición se interpola entre los dos
    vecinos (a menos de 2 h entre sí) o, al principio y al final del vuelo, se extrapola con la velocidad de
    los dos informes más cercanos (hasta 45 min para t − 3600 s; así se detectan también los informes de los
    extremos, añadido en la auditoría del 3-11-2026). Solo se corrige si el minuto es >= `min_minuto` (el mecanismo no puede dar otros).
    El avión se identifica por icao24 o, si viene truncado, por la matrícula.
    Añade la columna booleana `t_corregida_1h`; el ETA de la etiqueta 13 se desplaza igual."""
    out = df.copy()
    out["t_corregida_1h"] = False
    minuto = (out["t"] // 60 % 60).to_numpy()
    clave = out["icao24"].astype(object).where(out["icao24"].notna(), out["matricula"].astype(object))
    for hx, g in out[clave.notna()].sort_values("t").groupby(clave[clave.notna()], sort=False):
        tt, la, lo = g["t"].to_numpy(float), g["lat"].to_numpy(float), g["lon"].to_numpy(float)
        if len(tt) < 3:
            continue
        pos = out.index.get_indexer(g.index)
        for i in range(len(tt)):
            if minuto[pos[i]] < min_minuto:
                continue
            m = np.ones(len(tt), bool); m[i] = False
            T, L, O = tt[m], la[m], lo[m]

            def error(tc, max_ext):
                """Distancia [km] del informe a la posición del avión en tc según los demás informes:
                interpolada entre dos vecinos, o extrapolada (velocidad constante de los dos informes
                más cercanos) hasta `max_ext` s fuera del intervalo observado."""
                j = np.searchsorted(T, tc)
                if 0 < j < len(T):
                    a, b = j - 1, j
                    if T[b] - T[a] > 2 * 3600:
                        return np.nan
                elif len(T) >= 2:
                    a, b = (len(T) - 2, len(T) - 1) if j == len(T) else (0, 1)
                    lejos = tc - T[b] if j == len(T) else T[a] - tc
                    if T[b] - T[a] > 2 * 3600 or T[b] == T[a] or lejos > max_ext:
                        return np.nan
                else:
                    return np.nan
                w = (tc - T[a]) / (T[b] - T[a])
                return float(distancia(la[i], lo[i], L[a] + w * (L[b] - L[a]), O[a] + w * (O[b] - O[a]))) / 1e3
            e0, e1 = error(tt[i], 3 * 3600), error(tt[i] - 3600, 45 * 60)
            if e0 > err_max_km and e1 < err_ok_km:
                out.loc[g.index[i], "t_corregida_1h"] = True
    k = out["t_corregida_1h"]
    out.loc[k, "t"] -= 3600
    if "wp1_eta_t" in out:
        out.loc[k, "wp1_eta_t"] -= 3600
    out.attrs["n_corregidos_1h"] = int(k.sum())
    return out


def asignar_tipo(df: pd.DataFrame, aeronaves: pd.DataFrame) -> pd.DataFrame:
    """Añade typecode/operador. Se cruza por icao24 y, si falta (ICAO truncado), por matrícula
    normalizada, pero solo cuando la matrícula corresponde a un único icao24 en la base."""
    a = aeronaves.copy()
    a["icao24"] = a["icao24"].astype("string").str.lower()
    a["reg_n"] = a["registration"].astype("string").str.upper().str.replace("-", "", regex=False)
    por_icao = a.drop_duplicates("icao24").set_index("icao24")
    unicos = a.dropna(subset=["reg_n"]).groupby("reg_n")["icao24"].nunique()
    por_reg = a[a["reg_n"].isin(unicos[unicos == 1].index)].drop_duplicates("reg_n").set_index("reg_n")

    out = df.copy()
    out["reg_n"] = out["matricula"].str.upper().str.replace("-", "", regex=False)
    falta = out["icao24"].isna()
    out.loc[falta, "icao24"] = out.loc[falta, "reg_n"].map(por_reg["icao24"])
    out["icao24_por_matricula"] = falta & out["icao24"].notna()
    for col, src in (("typecode", "typecode"), ("operador", "operator"), ("operador_icao", "operatorIcao")):
        out[col] = out["icao24"].map(por_icao[src])
    return out


def segmentar_vuelos(df: pd.DataFrame, hueco_h: float = 3.0) -> pd.DataFrame:
    """Un vuelo = informes consecutivos de una misma aeronave separados por menos de `hueco_h` horas."""
    df = df.sort_values(["icao24", "t"]).copy()
    ic = df["icao24"].astype(object).to_numpy()
    t = df["t"].to_numpy(dtype=float)
    nuevo = np.ones(len(df), dtype=bool)
    nuevo[1:] = (ic[1:] != ic[:-1]) | (np.diff(t) > hueco_h * 3600)
    df["vuelo_id"] = np.cumsum(nuevo)
    return df


def resumen_vuelos(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby("vuelo_id")
    r = g.agg(icao24=("icao24", "first"), typecode=("typecode", "first"), operador=("operador", "first"),
              t0=("t", "min"), t1=("t", "max"), n=("t", "size"),
              lon0=("lon", "first"), lon1=("lon", "last"), lon_min=("lon", "min"), lon_max=("lon", "max"),
              lat_min=("lat", "min"), lat_max=("lat", "max"), fl_mediana=("alt_ft", lambda x: np.nanmedian(x) / 100),
              n_mach=("mach", "count"), n_meteo=("viento_kt", "count"))
    r["indicativo"] = g["vuelo"].agg(lambda x: x.dropna().mode().iat[0] if x.notna().any() else pd.NA)
    r["intervalo_mediano_min"] = g["t"].agg(lambda x: np.median(np.diff(np.sort(x.values))) / 60 if len(x) > 1 else np.nan)
    r["dur_h"] = (r["t1"] - r["t0"]) / 3600
    r["span_lon"] = r["lon_max"] - r["lon_min"]
    r["sentido"] = np.where(r["lon1"] > r["lon0"], "E", "O")
    r["fecha"] = pd.to_datetime(r["t0"], unit="s", utc=True).dt.date
    return r


def cierre_triangulo(df: pd.DataFrame) -> pd.DataFrame:
    """Coherencia interna de un informe con etiquetas 14+15+16: GS⃗ ≈ TAS·ĥ + W⃗.

    TAS = M·a(T_medida). El viento se da en convención meteorológica (dirección DE la que sopla).
    Devuelve el residuo vectorial en kt. El rumbo y la derrota se dan respecto al norte verdadero.
    """
    d = df.dropna(subset=["mach", "rumbo_deg", "gs_kt", "derrota_deg", "viento_kt", "viento_dir_deg", "temp_c"]).copy()
    tas_kt = d["mach"] * atm.vel_sonido(d["temp_c"] + 273.15) / KT
    h, trk, wd = np.radians(d["rumbo_deg"]), np.radians(d["derrota_deg"]), np.radians(d["viento_dir_deg"])
    ax, ay = tas_kt * np.sin(h), tas_kt * np.cos(h)
    wx, wy = -d["viento_kt"] * np.sin(wd), -d["viento_kt"] * np.cos(wd)
    gx, gy = d["gs_kt"] * np.sin(trk), d["gs_kt"] * np.cos(trk)
    d["tas_kt"] = tas_kt
    d["residuo_kt"] = np.hypot(ax + wx - gx, ay + wy - gy)
    d["residuo_gs_kt"] = np.hypot(ax + wx, ay + wy) - d["gs_kt"]
    return d
