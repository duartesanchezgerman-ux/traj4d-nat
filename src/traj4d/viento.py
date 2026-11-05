"""Campos meteorológicos: viento (u, v) y temperatura en función de (t, FL, lat, lon).

Interfaz común:  campo.muestra(t_unix, fl, lat, lon) -> (u [m/s], v [m/s], dT [K])
   u: componente este; v: componente norte; dT: desviación de temperatura respecto a ISA a ese FL.

Implementaciones:
  - CampoUniforme   : constante (para verificar).
  - ChorroSintetico : corriente en chorro gaussiana zonal (para pruebas controladas sin datos).
  - CampoERA5       : NetCDF de ERA5 en niveles de presión, descargado con scripts/datos/descargar_era5.py.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import atmosfera as atm


class Campo:
    def muestra(self, t, fl, lat, lon):
        raise NotImplementedError


@dataclass
class CampoUniforme(Campo):
    u: float = 0.0
    v: float = 0.0
    dT: float = 0.0

    def muestra(self, t, fl, lat, lon):
        return self.u, self.v, self.dT


@dataclass
class ChorroSintetico(Campo):
    """Chorro del oeste: u = U_max · exp(−((lat − lat0)/σ_lat)²) · exp(−((FL − FL0)/σ_FL)²).

    Valores por defecto: un chorro invernal típico del Atlántico Norte, con un máximo de unos
    75 m/s (≈145 kt) a FL340 y 45°N. Sirve para pruebas controladas; no es un dato meteorológico.
    """
    u_max: float = 75.0
    lat0: float = 45.0
    sigma_lat: float = 6.0
    fl0: float = 340.0
    sigma_fl: float = 80.0

    def muestra(self, t, fl, lat, lon):
        u = self.u_max * np.exp(-((lat - self.lat0) / self.sigma_lat) ** 2) * np.exp(-((fl - self.fl0) / self.sigma_fl) ** 2)
        return float(u), 0.0, 0.0


class CampoERA5(Campo):
    """Interpolación lineal en (tiempo, log p, lat, lon) sobre ficheros NetCDF de ERA5.

    Acepta los nombres de coordenadas del CDS actual (valid_time, pressure_level) y los del antiguo
    (time, level). Variables esperadas: u, v, t.
    """

    def __init__(self, ficheros):
        import xarray as xr
        from scipy.interpolate import RegularGridInterpolator

        ds = xr.open_mfdataset(ficheros, combine="by_coords") if isinstance(ficheros, (list, tuple)) \
            else xr.open_dataset(ficheros)
        ren = {k: v for k, v in {"valid_time": "time", "pressure_level": "level"}.items() if k in ds.dims or k in ds.coords}
        ds = ds.rename(ren).sortby(["time", "level", "latitude", "longitude"])
        tiempo = ds["time"].values.astype("datetime64[s]").astype(np.int64).astype(float)
        logp = np.log(ds["level"].values.astype(float) * 100.0)        # hPa -> Pa
        ejes = (tiempo, logp, ds["latitude"].values.astype(float), ds["longitude"].values.astype(float))
        kw = dict(bounds_error=True, method="linear")
        self._u = RegularGridInterpolator(ejes, ds["u"].values.astype(float), **kw)
        self._v = RegularGridInterpolator(ejes, ds["v"].values.astype(float), **kw)
        self._t = RegularGridInterpolator(ejes, ds["t"].values.astype(float), **kw)
        self.limites = {k: (float(a.min()), float(a.max())) for k, a in zip(("t", "logp", "lat", "lon"), ejes)}

    def muestra_lote(self, t, fl, lat, lon):
        """Versión vectorizada. Los puntos fuera del dominio (tiempo, nivel o área) devuelven NaN."""
        t, fl, lat, lon = (np.asarray(a, float) for a in (t, fl, lat, lon))
        h = atm.fl_a_m(fl)
        x = np.column_stack([t, np.log(atm.presion(h)), lat, lon])
        dentro = np.ones(len(x), bool)
        for k, (a, b) in enumerate(self.limites.values()):
            dentro &= (x[:, k] >= a) & (x[:, k] <= b)
        u = np.full(len(x), np.nan); v = u.copy(); T = u.copy()
        if dentro.any():
            u[dentro], v[dentro], T[dentro] = self._u(x[dentro]), self._v(x[dentro]), self._t(x[dentro])
        return u, v, T - atm.temperatura(h)

    def muestra(self, t, fl, lat, lon):
        h = atm.fl_a_m(fl)
        x = np.array([[float(t), float(np.log(atm.presion(h))), float(lat), float(lon)]])
        T = float(self._t(x)[0])
        return float(self._u(x)[0]), float(self._v(x)[0]), T - float(atm.temperatura(h))
