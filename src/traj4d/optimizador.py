"""Optimizador lateral por programación dinámica sobre una rejilla alineada con la ortodrómica.

Rejilla: N+1 etapas equiespaciadas a lo largo de la ortodrómica origen→destino. En cada etapa
intermedia hay 2K+1 nodos desplazados perpendicularmente k·Δ (k = −K…K). El origen y el destino
son nodos únicos.
Arcos: de (i, j) a (i+1, j') con |j' − j| ≤ salto_max. Esto limita el ángulo de giro entre tramos.
Coste de un arco: J = combustible + CI·tiempo, evaluado en el punto medio del arco con el viento
4D en el instante estimado de paso.

La DP hacia delante guarda en cada nodo la mejor etiqueta (coste, masa, tiempo, predecesor).
Como la masa y el tiempo dependen del camino, esta DP es una aproximación estándar de "etiqueta
dominante": el error es de segundo orden porque la masa varía poco dentro de un arco. La ruta
resultante se re-evalúa SIEMPRE con el integrador RK4 completo (`integrar_ruta`), y esa es la cifra
que se reporta.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import atmosfera as atm
from .ruta import R_TIERRA, Ruta, distancia, rumbo_inicial


def destino_v(lat, lon, rumbo, d):
    """Versión vectorizada de `destino` (arrays numpy)."""
    p1, l1, dl = np.radians(lat), np.radians(lon), np.asarray(d) / R_TIERRA
    p2 = np.arcsin(np.sin(p1) * np.cos(dl) + np.cos(p1) * np.sin(dl) * np.cos(rumbo))
    l2 = l1 + np.arctan2(np.sin(rumbo) * np.sin(dl) * np.cos(p1), np.cos(dl) - np.sin(p1) * np.sin(p2))
    return np.degrees(p2), (np.degrees(l2) + 540) % 360 - 180


def destino(lat, lon, rumbo, d):
    """Punto a distancia d [m] con rumbo inicial `rumbo` [rad] (esfera)."""
    la, lo = destino_v(lat, lon, rumbo, d)
    return float(la), float(lo)


def _muestra_lote(campo, t, fl, lat, lon):
    if hasattr(campo, "muestra_lote"):
        u, v, dT = campo.muestra_lote(t, np.full(len(t), fl), lat, lon)
        if np.isnan(u).any():
            raise ValueError("One of the requested xi is out of bounds (rejilla fuera del dominio del campo)")
        return u, v, dT
    out = np.array([campo.muestra(a, fl, b, c) for a, b, c in zip(t, lat, lon)], float).reshape(-1, 3)
    return out[:, 0], out[:, 1], out[:, 2]


@dataclass
class ResultadoDP:
    ruta: Ruta
    desplazamientos_km: np.ndarray
    coste_dp: float
    lat_nodos: np.ndarray
    lon_nodos: np.ndarray


def optimizar_lateral(modelo, campo, lat0, lon0, lat1, lon1, fl, mach, m0, t0=0.0, ci=0.0,
                      n_etapas=24, K=20, delta_km=25.0, salto_max=2) -> ResultadoDP:
    gc = Ruta([lat0, lat1], [lon0, lon1])
    L = gc.longitud
    lat = np.full((n_etapas + 1, 2 * K + 1), np.nan)
    lon = np.full_like(lat, np.nan)
    for i in range(n_etapas + 1):
        la, lo, chi = gc.punto(L * i / n_etapas)
        if i in (0, n_etapas):
            lat[i, K], lon[i, K] = la, lo
            continue
        for j, k in enumerate(range(-K, K + 1)):
            lat[i, j], lon[i, j] = destino(la, lo, chi + np.pi / 2, k * delta_km * 1e3)  # +k = derecha (estribor)

    INF = np.inf
    coste = np.full(lat.shape, INF)
    masa = np.full(lat.shape, np.nan)
    tiempo = np.full(lat.shape, np.nan)
    pred = np.full(lat.shape, -1, dtype=int)
    coste[0, K], masa[0, K], tiempo[0, K] = 0.0, m0, t0

    for i in range(n_etapas):
        origenes = np.array([K]) if i == 0 else np.arange(2 * K + 1)
        destinos = np.array([K]) if i + 1 == n_etapas else np.arange(2 * K + 1)
        J, JJ = np.meshgrid(origenes, destinos, indexing="ij")
        J, JJ = J.ravel(), JJ.ravel()
        ok = np.isfinite(coste[i, J])
        if 0 < i and i + 1 < n_etapas:
            ok &= np.abs(JJ - J) <= salto_max
        if i == 0:
            ok &= np.abs(JJ - K) <= salto_max
        if i + 1 == n_etapas:
            ok &= np.abs(J - K) <= salto_max
        J, JJ = J[ok], JJ[ok]
        la_a, lo_a, la_b, lo_b = lat[i, J], lon[i, J], lat[i + 1, JJ], lon[i + 1, JJ]
        d = distancia(la_a, lo_a, la_b, lo_b)
        chi = rumbo_inicial(la_a, lo_a, la_b, lo_b)
        lam, lom = destino_v(la_a, lo_a, chi, d / 2)
        u, v, dT = _muestra_lote(campo, tiempo[i, J], fl, lam, lom)
        e = modelo.estado(masa[i, J], np.full(len(J), mach), np.full(len(J), fl), dT)
        tas = np.asarray(e["TAS"], float)
        wl = u * np.sin(chi) + v * np.cos(chi)
        wt = u * np.cos(chi) - v * np.sin(chi)
        vg = wl + np.sqrt(np.maximum(tas**2 - wt**2, 1.0))
        dt = d / vg
        comb = np.asarray(e["FF"], float) * dt
        c = coste[i, J] + comb + ci / 60.0 * dt
        # para cada destino, el mejor origen (en empate, el primero: igual que la versión escalar)
        orden = np.lexsort((J, c, JJ))
        primero = np.r_[True, JJ[orden][1:] != JJ[orden][:-1]]
        sel = orden[primero]
        mejora = c[sel] < coste[i + 1, JJ[sel]]
        sel = sel[mejora]
        coste[i + 1, JJ[sel]] = c[sel]
        masa[i + 1, JJ[sel]] = masa[i, J[sel]] - comb[sel]
        tiempo[i + 1, JJ[sel]] = tiempo[i, J[sel]] + dt[sel]
        pred[i + 1, JJ[sel]] = J[sel]

    j = K
    camino = [j]
    for i in range(n_etapas, 0, -1):
        j = pred[i, j]
        camino.append(j)
    camino = camino[::-1]
    la = np.array([lat[i, j] for i, j in enumerate(camino)])
    lo = np.array([lon[i, j] for i, j in enumerate(camino)])
    return ResultadoDP(Ruta(la, lo), (np.array(camino) - K) * delta_km, float(coste[n_etapas, K]), lat, lon)
