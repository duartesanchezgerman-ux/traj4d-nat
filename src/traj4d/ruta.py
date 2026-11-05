"""Geometría de rutas sobre una Tierra esférica (radio medio 6.371.008,8 m).

Limitación declarada: frente al elipsoide WGS-84, la esfera comete un error de distancia del
orden de ±0,3 %. Ese error es casi idéntico para las rutas que se comparan entre sí en el mismo
corredor, así que afecta poco a los ahorros relativos. Si hiciera falta, se sustituye por
geodésicas WGS-84 (pyproj.Geod).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

R_TIERRA = 6_371_008.8

AEROPUERTOS = {  # (lat, lon) del punto de referencia del aeródromo, en grados
    "LEMD": (40.4719, -3.5626),
    "EGLL": (51.4706, -0.4619),
    "KJFK": (40.6398, -73.7789),
    "KBOS": (42.3656, -71.0096),
}


def _xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def _latlon(v):
    v = v / np.linalg.norm(v)
    return float(np.degrees(np.arcsin(v[2]))), float(np.degrees(np.arctan2(v[1], v[0])))


def distancia(lat1, lon1, lat2, lon2):
    """Distancia ortodrómica [m] (fórmula de haversine)."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_TIERRA * np.arcsin(np.sqrt(a))


def rumbo_inicial(lat1, lon1, lat2, lon2):
    """Rumbo verdadero inicial [rad, desde el norte, en sentido horario]."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    return np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))


@dataclass
class Ruta:
    """Polilínea de puntos de paso unidos por arcos de círculo máximo."""
    lat: np.ndarray
    lon: np.ndarray

    def __post_init__(self):
        self.lat, self.lon = np.asarray(self.lat, float), np.asarray(self.lon, float)
        tramos = distancia(self.lat[:-1], self.lon[:-1], self.lat[1:], self.lon[1:])
        self.s_nodos = np.concatenate([[0.0], np.cumsum(tramos)])

    @classmethod
    def ortodromica(cls, origen: str, destino: str, n: int = 1):
        (a, b), (c, d) = AEROPUERTOS[origen], AEROPUERTOS[destino]
        if n == 1:
            return cls([a, c], [b, d])
        v1, v2 = _xyz(a, b), _xyz(c, d)
        om = np.arccos(np.clip(v1 @ v2, -1, 1))
        pts = [_latlon(np.sin((1 - f) * om) * v1 + np.sin(f * om) * v2) for f in np.linspace(0, 1, n + 1)]
        return cls([p[0] for p in pts], [p[1] for p in pts])

    @property
    def longitud(self):
        return float(self.s_nodos[-1])

    def punto(self, s):
        """(lat, lon, rumbo_verdadero[rad]) a la distancia s [m] desde el origen."""
        s = float(np.clip(s, 0.0, self.longitud))
        i = int(np.clip(np.searchsorted(self.s_nodos, s, side="right") - 1, 0, len(self.lat) - 2))
        L = self.s_nodos[i + 1] - self.s_nodos[i]
        f = 0.0 if L == 0 else (s - self.s_nodos[i]) / L
        v1, v2 = _xyz(self.lat[i], self.lon[i]), _xyz(self.lat[i + 1], self.lon[i + 1])
        om = np.arccos(np.clip(v1 @ v2, -1, 1))

        def slerp(g):
            if om < 1e-12:
                return float(self.lat[i]), float(self.lon[i])
            return _latlon(np.sin((1 - g) * om) / np.sin(om) * v1 + np.sin(g * om) / np.sin(om) * v2)

        lat, lon = slerp(f)
        # Rumbo local = dirección de la tangente al arco (diferencia finita de 1e-5 del tramo)
        d = 1e-5
        if f + d <= 1.0:
            rumbo = rumbo_inicial(lat, lon, *slerp(f + d))
        else:
            rumbo = rumbo_inicial(*slerp(f - d), lat, lon)
        return lat, lon, float(rumbo)


def triangulo_viento(tas, rumbo_ruta, u, v):
    """Velocidad sobre el suelo [m/s] cuando la aeronave mantiene la derrota `rumbo_ruta`.

    u: componente este del viento y v: componente norte [m/s].
    Componente longitudinal: w_l = u·sin χ + v·cos χ. Componente transversal: w_t = u·cos χ − v·sin χ.
    Vg = w_l + sqrt(TAS² − w_t²). Para mantener la derrota, el avión corrige con un ángulo de deriva.
    """
    wl = u * np.sin(rumbo_ruta) + v * np.cos(rumbo_ruta)
    wt = u * np.cos(rumbo_ruta) - v * np.sin(rumbo_ruta)
    return wl + np.sqrt(np.maximum(tas**2 - wt**2, 0.0))


def puntos_distintos(lat, lon, min_m: float = 5e3):
    """Índices de los informes que se usan como vértices de la ruta volada.

    Se descartan los informes a menos de `min_m` del vértice anterior (informes repetidos en el mismo sitio),
    pero el PRIMERO y el ÚLTIMO se conservan siempre, para que la ruta volada empiece y acabe en los mismos
    puntos que la ortodrómica y el óptimo con los que se compara. Si el último queda a menos de `min_m` del
    vértice anterior, se quita ese vértice anterior (salvo que sea el primero).
    (Corrige un defecto hallado en la auditoría del 3-11-2026: antes se podía perder el último informe.)"""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    idx = [0]
    for i in range(1, len(lat)):
        if float(distancia(lat[idx[-1]], lon[idx[-1]], lat[i], lon[i])) > min_m:
            idx.append(i)
    if idx[-1] != len(lat) - 1:
        if len(idx) > 1:
            idx.pop()
        idx.append(len(lat) - 1)
    return np.array(idx)
