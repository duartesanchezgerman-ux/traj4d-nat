"""Lector del volcado de texto ADS-C de OpenSky (Zenodo, doi:10.5281/zenodo.14659997).

Solo usa la librería estándar (Python ≥ 3.10), así que funciona en cualquier máquina, también sin
numpy ni pandas.

Estructura del fichero (adsc_file_structure.md): cada mensaje empieza con una línea
"Registration: … ICAO ID: … ATSU Address: CÓDIGO - Centro", sigue con "Channel Frequency", la
trama cruda "/…" y una o varias etiquetas "Tag NN …" con sus campos indentados, y termina en
"CRC: …". Los mensajes se separan con una línea en blanco.

Cada mensaje CON posición produce una fila. La posición puede venir en las etiquetas 07 (informe
básico), 09 (emergencia), 10, 18, 19 o 20 (informes por evento), que tienen todas el mismo formato.
Las etiquetas 14 (referencia tierra), 15 (referencia aire), 16 (meteo), 12 (indicativo) y 13
(ruta prevista) no llevan hora propia y se asignan a la posición de su mismo mensaje.

Convenciones:
  * Hora: el campo Timestamp (y el ETA de la etiqueta 13) del volcado NO está en UTC sino en hora
    local de Europa central (CET = UTC+1 en invierno, CEST = UTC+2 en verano), aunque ADS-C/FANS
    transmite UTC. Se comprobó el 3-11-2026 cruzando 95 informes de 61 aviones con trazas ADS-B
    (adsb.lol, UTC) en el mismo punto: desfase de −3600 s en noviembre y −7200 s en julio (±15 s).
    _ts convierte de Europe/Zurich a UTC y guarda segundos Unix. En el volcado no hay informes en
    las horas ambiguas o inexistentes de los cambios de hora; si aparecieran, se rechazan.
  * Viento: "True wind direction" es la dirección DE la que sopla (convención meteorológica).
    Componentes:  u = −V·sin(dir),  v = −V·cos(dir).
"""
from __future__ import annotations

import gzip
import re
from datetime import datetime, timedelta, timezone
from typing import Iterator
from zoneinfo import ZoneInfo

ETIQUETAS_POSICION = ("07", "09", "10", "18", "19", "20")

COLUMNAS = [
    "t", "icao24", "icao24_bruto", "matricula", "atsu", "centro", "frecuencia_hz", "tipo_informe", "n_posiciones",
    "lat", "lon", "alt_ft", "precision", "vuelo",
    "derrota_deg", "gs_kt", "vs_tierra_fpm",
    "rumbo_deg", "mach", "vs_aire_fpm",
    "viento_kt", "viento_dir_deg", "temp_c",
    "wp1_lat", "wp1_lon", "wp1_alt_ft", "wp1_eta_t", "wp2_lat", "wp2_lon", "wp2_alt_ft",
]

# El ICAO ID suele tener 6 dígitos hex, pero en ~6,5 % de los mensajes viene truncado (p. ej. "0a48"):
# es un defecto del decodificador de origen. Se guarda en bruto y icao24 queda vacío; se recupera
# después a través de la matrícula y la base de aeronaves.
_CAB = re.compile(r"^Registration:\s*(\S*)\s+ICAO ID:\s*([0-9A-Fa-f]*)\s*ATSU Address:\s*(\S+)\s*-\s*(.*?)\s*$")
_TAG = re.compile(r"^Tag (\d{2}) ")


ZONA_VOLCADO = ZoneInfo("Europe/Zurich")      # hora local en que viene el volcado (ver arriba)


def _ts(texto: str) -> float:
    """'2023-07-07 15:07:11.000' (hora local de Europa central) -> segundos Unix (UTC).

    Rechaza las horas inexistentes (salto de primavera) y las ambiguas (hora repetida de otoño)."""
    local = datetime.strptime(texto.strip(), "%Y-%m-%d %H:%M:%S.%f")
    a, b = (local.replace(tzinfo=ZONA_VOLCADO, fold=f) for f in (0, 1))
    if a.utcoffset() != b.utcoffset():
        raise ValueError(f"hora ambigua o inexistente en Europe/Zurich: {texto!r}")
    ida_vuelta = a.astimezone(timezone.utc).astimezone(ZONA_VOLCADO).replace(tzinfo=None)
    if ida_vuelta != local:
        raise ValueError(f"hora inexistente en Europe/Zurich: {texto!r}")
    return a.timestamp()


def _num(texto: str) -> float:
    """'35000.0 ft' / '456.5 kt' / '-51.0 C' / '0.837' -> float."""
    return float(texto.strip().split()[0])


class Estadisticas(dict):
    def sumar(self, k, n=1):
        self[k] = self.get(k, 0) + n


def leer_mensajes(lineas, est: Estadisticas | None = None) -> Iterator[dict]:
    """Recorre las líneas del volcado y produce un dict por cada mensaje con posición."""
    est = est if est is not None else Estadisticas()
    msg = None
    tag = None
    sub = None          # subcontexto de la etiqueta 13: "wp1" o "wp2"

    def cerrar(m):
        if m is None:
            return None
        est.sumar("mensajes")
        if m.get("lat") is None:
            est.sumar("mensajes_sin_posicion")
            return None
        est.sumar("mensajes_con_posicion")
        return m

    for linea in lineas:
        if linea.startswith("Registration:"):
            fila = cerrar(msg)
            if fila is not None:
                yield fila
            c = _CAB.match(linea.rstrip("\n"))
            if c is None:
                est.sumar("cabeceras_ilegibles")
                msg, tag = None, None
                continue
            msg = dict.fromkeys(COLUMNAS)
            bruto = (c.group(2) or "").lower()
            if len(bruto) != 6:
                est.sumar("icao24_truncado")
            msg.update(matricula=c.group(1) or None, icao24=bruto if len(bruto) == 6 else None,
                       icao24_bruto=bruto or None, atsu=c.group(3), centro=c.group(4).strip(), n_posiciones=0)
            tag, sub = None, None
            continue
        if msg is None:
            continue
        s = linea.strip()
        if not s:
            continue
        if linea.startswith("Channel Frequency:"):
            msg["frecuencia_hz"] = int(float(s.split(":", 1)[1]))
            continue
        if linea.startswith("Tag "):
            m = _TAG.match(linea)
            tag = m.group(1) if m else None
            sub = None
            if tag in ETIQUETAS_POSICION:
                msg["n_posiciones"] += 1
                if msg["n_posiciones"] > 1:
                    est.sumar("posiciones_extra_ignoradas")
                    tag = "_ignorar"
                else:
                    msg["tipo_informe"] = tag
            continue
        if linea.startswith("CRC:") or linea.startswith("/") or linea.startswith("Unknown tag"):
            tag = None if not linea.startswith("/") else tag
            continue
        if ":" not in s:
            continue
        clave, valor = (x.strip() for x in s.split(":", 1))
        try:
            if tag in ETIQUETAS_POSICION:
                if clave == "Latitude": msg["lat"] = float(valor)
                elif clave == "Longitude": msg["lon"] = float(valor)
                elif clave == "Altitude": msg["alt_ft"] = _num(valor)
                elif clave == "Timestamp": msg["t"] = _ts(valor)
                elif clave == "Position accuracy": msg["precision"] = valor
            elif tag == "12" and clave == "Flight ID":
                msg["vuelo"] = valor.strip() or None
            elif tag == "14":
                if clave == "True track": msg["derrota_deg"] = _num(valor)
                elif clave == "Ground speed": msg["gs_kt"] = _num(valor)
                elif clave == "Vertical speed": msg["vs_tierra_fpm"] = _num(valor)
            elif tag == "15":
                if clave == "True heading": msg["rumbo_deg"] = _num(valor)
                elif clave == "Mach speed": msg["mach"] = _num(valor)
                elif clave == "Vertical speed": msg["vs_aire_fpm"] = _num(valor)
            elif tag == "16":
                if clave == "Wind speed": msg["viento_kt"] = _num(valor)
                elif clave == "True wind direction": msg["viento_dir_deg"] = _num(valor)
                elif clave == "Temperature": msg["temp_c"] = _num(valor)
            elif tag == "13":
                if clave == "Next waypoint" or s.startswith("Next waypoint"):
                    sub = "wp1"
                elif clave == "Next+1 waypoint" or s.startswith("Next+1 waypoint"):
                    sub = "wp2"
                elif sub is not None:
                    k = {"Lat": "lat", "Lon": "lon", "Alt": "alt_ft", "ETA": "eta_t"}.get(clave)
                    if k == "eta_t":
                        if sub == "wp1":
                            msg["wp1_eta_t"] = _ts(valor)
                    elif k is not None:
                        msg[f"{sub}_{k}"] = _num(valor) if k == "alt_ft" else float(valor)
        except ValueError:
            est.sumar(f"valor_ilegible_tag{tag}_{clave}")

    fila = cerrar(msg)
    if fila is not None:
        yield fila


def abrir(ruta):
    return gzip.open(ruta, "rt", encoding="utf-8", errors="replace") if str(ruta).endswith(".gz") \
        else open(ruta, encoding="utf-8", errors="replace")
