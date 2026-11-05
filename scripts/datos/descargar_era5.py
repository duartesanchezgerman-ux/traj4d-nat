"""Descarga ERA5 en niveles de presión (viento U/V, temperatura, geopotencial) sobre el Atlántico Norte.

Requisitos previos (una sola vez):
  1. pip install cdsapi
  2. En https://cds.climate.copernicus.eu entra en tu perfil y copia tu "Personal Access Token".
  3. Crea el fichero ~/.cdsapirc con estas dos líneas:
         url: https://cds.climate.copernicus.eu/api
         key: TU_TOKEN
  4. Acepta la licencia del conjunto de datos: abre "ERA5 hourly data on pressure levels from 1940
     to present", pestaña "Download", y acepta "Terms of use" al final de la página.
     Sin esto la descarga falla con un error 403.

Uso (prueba de un solo día, unos pocos MB):
    python scripts/datos/descargar_era5.py --dia 2024-01-15
Rango de fechas (un fichero por día):
    python scripts/datos/descargar_era5.py --desde 2024-01-01 --hasta 2024-01-31
Lista de días y área reducida (N O S E), p. ej. los días de estudio de la Fase 2:
    python scripts/datos/descargar_era5.py --dias 2023-11-01 2023-11-03 --area 70 -60 30 0
Modo agrupado (RECOMENDADO para muchos días): una sola petición por mes con todos sus días. La cola de
Copernicus cobra por petición, no por día. El archivo mensual se divide después en un archivo por
día (mismo formato que antes) y se borra:
    python scripts/datos/descargar_era5.py --agrupar-por-mes --dias 2023-07-11 2023-07-12 ... --area 70 -60 30 0
"""
import argparse
from datetime import date, timedelta
from pathlib import Path

# Niveles de crucero: 150 hPa ≈ FL450, 350 hPa ≈ FL265. Cubre toda la envolvente del 787-9.
NIVELES = ["150", "175", "200", "225", "250", "300", "350"]
VARIABLES = ["u_component_of_wind", "v_component_of_wind", "temperature", "geopotential"]
AREA = [70, -90, 20, 10]          # N, O, S, E: Atlántico Norte con margen (MAD/LHR -> costa este de EE. UU.)
HORAS = [f"{h:02d}:00" for h in range(24)]


def descargar(dia: date, salida: Path, area=None):
    import cdsapi
    destino = salida / f"era5_pl_{dia.isoformat()}.nc"
    if destino.exists():
        print(f"ya existe {destino}")
        return
    peticion = {
        "product_type": ["reanalysis"],
        "variable": VARIABLES,
        "year": [f"{dia.year}"], "month": [f"{dia.month:02d}"], "day": [f"{dia.day:02d}"],
        "time": HORAS,
        "pressure_level": NIVELES,
        "area": list(area) if area else AREA,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }
    cdsapi.Client().retrieve("reanalysis-era5-pressure-levels", peticion).download(str(destino))
    print(f"OK {destino}")


def reintentar(accion, espera_s=300, max_intentos=48):
    """Copernicus rechaza peticiones nuevas si ya tienes el máximo en cola ("Maximum number queued
    requests"). Las peticiones de un script cortado con Ctrl+C SIGUEN en la cola del servidor.
    En ese caso se espera y se reintenta (por defecto cada 5 min, hasta 4 h). Otros errores se propagan."""
    import time
    for intento in range(1, max_intentos + 1):
        try:
            return accion()
        except Exception as e:  # noqa: BLE001
            if "queued requests" not in str(e) and "temporarily" not in str(e).lower():
                raise
            print(f"  cola de Copernicus llena (intento {intento}/{max_intentos}); reintento en {espera_s // 60} min. "
                  f"Puedes cancelar peticiones antiguas en https://cds.climate.copernicus.eu/requests")
            time.sleep(espera_s)
    raise RuntimeError("La cola de Copernicus sigue llena tras varias horas; revisa tus peticiones en la web.")


def descargar_mes(dias, salida: Path, area=None):
    """Una petición con todos los días pendientes del mismo mes. Después divide por día."""
    import cdsapi
    import xarray as xr
    pendientes = [d for d in dias if not (salida / f"era5_pl_{d.isoformat()}.nc").exists()]
    if not pendientes:
        print(f"ya existen todos los días de {dias[0]:%Y-%m}")
        return
    y, m = pendientes[0].year, pendientes[0].month
    assert all((d.year, d.month) == (y, m) for d in pendientes)
    tmp = salida / f"_mes_{y}-{m:02d}.nc"
    peticion = {
        "product_type": ["reanalysis"], "variable": VARIABLES,
        "year": [f"{y}"], "month": [f"{m:02d}"], "day": sorted({f"{d.day:02d}" for d in pendientes}),
        "time": HORAS, "pressure_level": NIVELES, "area": list(area) if area else AREA,
        "data_format": "netcdf", "download_format": "unarchived",
    }
    print(f"petición {y}-{m:02d}: {len(pendientes)} días")
    reintentar(lambda: cdsapi.Client().retrieve("reanalysis-era5-pressure-levels", peticion).download(str(tmp)))
    dividir_por_dia(tmp, salida, pendientes)
    tmp.unlink()


def dividir_por_dia(archivo: Path, salida: Path, dias):
    """Divide un NetCDF multidía en un archivo por día. Comprueba que cada día tiene sus 24 horas."""
    import numpy as np
    import xarray as xr
    with xr.open_dataset(archivo) as ds:
        dim = "valid_time" if "valid_time" in ds.dims else "time"
        for d in dias:
            sub = ds.sel({dim: slice(np.datetime64(f"{d}T00:00"), np.datetime64(f"{d}T23:59"))})
            n = sub.sizes[dim]
            if n != 24:
                raise RuntimeError(f"{d}: {n} horas en lugar de 24; no se escribe el archivo")
            destino = salida / f"era5_pl_{d.isoformat()}.nc"
            sub.load().to_netcdf(destino)
            print(f"OK {destino}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dia")
    ap.add_argument("--dias", nargs="+")
    ap.add_argument("--area", nargs=4, type=float, metavar=("N", "O", "S", "E"))
    ap.add_argument("--agrupar-por-mes", action="store_true")
    ap.add_argument("--desde")
    ap.add_argument("--hasta")
    ap.add_argument("--salida", default="datos/era5")
    a = ap.parse_args()
    out = Path(a.salida)
    out.mkdir(parents=True, exist_ok=True)
    if a.dias and a.agrupar_por_mes:
        from itertools import groupby
        dias = sorted(date.fromisoformat(d) for d in a.dias)
        for _, grupo in groupby(dias, key=lambda d: (d.year, d.month)):
            descargar_mes(list(grupo), out, a.area)
    elif a.dias:
        for d in a.dias:
            descargar(date.fromisoformat(d), out, a.area)
    elif a.dia:
        descargar(date.fromisoformat(a.dia), out, a.area)
    else:
        d, fin = date.fromisoformat(a.desde), date.fromisoformat(a.hasta)
        while d <= fin:
            descargar(d, out, a.area)
            d += timedelta(days=1)
