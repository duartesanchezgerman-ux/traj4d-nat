"""Verificación de la Fase 1b: geometría de ruta, triángulo de viento, ERA5 y cost index.
Ejecutar:  PYTHONPATH=src pytest -q tests
"""
import numpy as np
import pytest

from traj4d import atmosfera as atm
from traj4d.crucero import (integrar_crucero, integrar_ruta, masa_final_breguet, optimo_local,
                            perfil_constante)
from traj4d.rendimiento import CongeladoModel, PollSchumannModel
from traj4d.ruta import Ruta, distancia, rumbo_inicial
from traj4d.viento import CampoERA5, CampoUniforme


# ------------------------------------------------------------------ Ruta
def test_ortodromica_mad_jfk():
    r1 = Ruta.ortodromica("LEMD", "KJFK")
    r20 = Ruta.ortodromica("LEMD", "KJFK", n=20)
    assert r1.longitud == pytest.approx(5_770e3, rel=5e-3)          # ≈ 3.115 NM
    assert r20.longitud == pytest.approx(r1.longitud, rel=1e-9)     # subdividir no cambia el arco
    lat, lon, chi = r20.punto(0.0)
    assert np.degrees(chi) == pytest.approx(np.degrees(rumbo_inicial(40.4719, -3.5626, 40.6398, -73.7789)), abs=1e-3)
    lat, lon, _ = r20.punto(r20.longitud)
    assert (lat, lon) == pytest.approx((40.6398, -73.7789), abs=1e-6)


def test_punto_medio_esta_a_mitad_de_distancia():
    r = Ruta.ortodromica("EGLL", "KJFK", n=7)
    lat, lon, _ = r.punto(r.longitud / 2)
    d1 = distancia(51.4706, -0.4619, lat, lon)
    assert d1 == pytest.approx(r.longitud / 2, rel=1e-9)


# -------------------------------------------- Ruta + campo = caso 1D
def test_ruta_sin_viento_igual_que_crucero_1d():
    mod = PollSchumannModel()
    r = Ruta.ortodromica("LEMD", "KJFK", n=10)
    a = integrar_ruta(mod, 225e3, r, perfil_constante(370, 0.85), CampoUniforme(), ds=20e3)
    b = integrar_crucero(mod, 225e3, r.longitud, perfil_constante(370, 0.85), ds=20e3)
    assert a.m[-1] == pytest.approx(b.m[-1], rel=1e-12)


def test_viento_longitudinal_en_ecuador():
    """Ruta hacia el oeste por el ecuador (derrota 270°) con u = +20 m/s: viento en contra puro."""
    mod = PollSchumannModel()
    r = Ruta([0.0, 0.0], [-10.0, -50.0])
    a = integrar_ruta(mod, 225e3, r, perfil_constante(370, 0.85), CampoUniforme(u=20.0), ds=20e3)
    b = integrar_crucero(mod, 225e3, r.longitud, perfil_constante(370, 0.85), viento=-20.0, ds=20e3)
    assert a.m[-1] == pytest.approx(b.m[-1], rel=1e-9)


def test_viento_cruzado_breguet_exacto():
    """Viento transversal puro: Vg = sqrt(TAS² − w_t²). Con el modelo congelado, Breguet es exacto."""
    mod = CongeladoModel(LD=19.0, TSFC=1.5e-5)
    r = Ruta([0.0, 0.0], [-10.0, -50.0])
    res = integrar_ruta(mod, 225e3, r, perfil_constante(370, 0.85), CampoUniforme(v=30.0), ds=20e3)
    tas = float(atm.mach_a_tas(0.85, atm.fl_a_m(370)))
    m1 = masa_final_breguet(225e3, r.longitud, np.sqrt(tas**2 - 30.0**2), 19.0, 1.5e-5)
    assert res.m[-1] == pytest.approx(m1, rel=1e-9)


# ------------------------------------------------------------- ERA5
def test_lector_era5_sintetico(tmp_path):
    """Fichero con la estructura del CDS actual (valid_time, pressure_level descendente, latitud
    descendente). Comprueba el renombrado, la ordenación y que ΔISA se recupera en un nivel de la rejilla."""
    import xarray as xr
    fls = np.array([300.0, 340.0, 380.0])
    p_hpa = atm.presion(atm.fl_a_m(fls)) / 100.0
    lat = np.array([60.0, 50.0, 40.0, 30.0])
    lon = np.array([-80.0, -40.0, 0.0])
    t = np.array(["2024-01-15T00", "2024-01-15T06"], dtype="datetime64[ns]")
    shp = (2, 3, 4, 3)
    T = np.broadcast_to((atm.temperatura(atm.fl_a_m(fls)) + 5.0)[None, :, None, None], shp).copy()
    ds = xr.Dataset({"u": (("valid_time", "pressure_level", "latitude", "longitude"), np.full(shp, 10.0)),
                     "v": (("valid_time", "pressure_level", "latitude", "longitude"), np.full(shp, -5.0)),
                     "t": (("valid_time", "pressure_level", "latitude", "longitude"), T)},
                    coords={"valid_time": t, "pressure_level": p_hpa, "latitude": lat, "longitude": lon})
    f = tmp_path / "era5.nc"
    ds.to_netcdf(f)
    c = CampoERA5(str(f))
    t_unix = float(np.datetime64("2024-01-15T03", "s").astype(np.int64))
    u, v, dT = c.muestra(t_unix, 340.0, 45.0, -20.0)
    assert (u, v) == pytest.approx((10.0, -5.0))
    assert dT == pytest.approx(5.0, abs=1e-9)


# ------------------------------------------------------- Cost index
def test_mach_econ_crece_con_ci():
    mod = PollSchumannModel()
    machs = [optimo_local(mod, None, 200e3, 0.0, ci, fl_fijo=370)[1] for ci in (0, 30, 60, 120)]
    assert all(b >= a for a, b in zip(machs, machs[1:])) and machs[-1] > machs[0]


def test_lector_era5_lote_igual_que_puntual_y_nan_fuera(tmp_path):
    import xarray as xr
    fls = np.array([300.0, 340.0, 380.0])
    p_hpa = atm.presion(atm.fl_a_m(fls)) / 100.0
    lat = np.array([60.0, 50.0, 40.0]); lon = np.array([-40.0, -20.0, 0.0])
    t = np.array(["2024-01-15T00", "2024-01-15T06"], dtype="datetime64[ns]")
    rng = np.random.default_rng(0)
    shp = (2, 3, 3, 3)
    ds = xr.Dataset({k: (("valid_time", "pressure_level", "latitude", "longitude"), rng.normal(0, 20, shp) + (220 if k == "t" else 0))
                     for k in ("u", "v", "t")},
                    coords={"valid_time": t, "pressure_level": p_hpa, "latitude": lat, "longitude": lon})
    f = tmp_path / "e.nc"; ds.to_netcdf(f)
    c = CampoERA5(str(f))
    t0 = float(np.datetime64("2024-01-15T00", "s").astype(np.int64))
    T = np.array([t0 + 3600, t0 + 7200, t0 + 99 * 3600])
    FL = np.array([320.0, 350.0, 350.0]); LA = np.array([45.0, 55.0, 45.0]); LO = np.array([-30.0, -10.0, -10.0])
    u, v, dT = c.muestra_lote(T, FL, LA, LO)
    for k in range(2):
        assert (u[k], v[k], dT[k]) == pytest.approx(c.muestra(T[k], FL[k], LA[k], LO[k]), rel=1e-12)
    assert np.isnan(u[2]) and np.isnan(dT[2])     # fuera del intervalo de tiempo


def test_dividir_era5_por_dia_preserva_datos(tmp_path):
    """El archivo mensual se divide en diarios idénticos a la porción original, cada uno con 24 h."""
    import importlib.util
    import xarray as xr
    from datetime import date
    spec = importlib.util.spec_from_file_location("dl", "scripts/datos/descargar_era5.py")
    dl = importlib.util.module_from_spec(spec); spec.loader.exec_module(dl)
    t = np.array([np.datetime64("2023-07-11T00") + np.timedelta64(h, "h") for h in range(24)] +
                 [np.datetime64("2023-07-13T00") + np.timedelta64(h, "h") for h in range(24)])
    rng = np.random.default_rng(1)
    ds = xr.Dataset({"u": (("valid_time", "pressure_level", "latitude", "longitude"), rng.normal(size=(48, 2, 3, 3)))},
                    coords={"valid_time": t, "pressure_level": [250.0, 200.0], "latitude": [50.0, 45.0, 40.0], "longitude": [-30.0, -20.0, -10.0]})
    f = tmp_path / "_mes.nc"; ds.to_netcdf(f)
    dl.dividir_por_dia(f, tmp_path, [date(2023, 7, 11), date(2023, 7, 13)])
    for k, d in enumerate(("2023-07-11", "2023-07-13")):
        sub = xr.open_dataset(tmp_path / f"era5_pl_{d}.nc")
        assert sub.sizes["valid_time"] == 24
        np.testing.assert_array_equal(sub.u.values, ds.u.values[24 * k:24 * (k + 1)])
    with pytest.raises(RuntimeError):
        dl.dividir_por_dia(f, tmp_path, [date(2023, 7, 12)])     # día sin datos -> error, no archivo vacío
