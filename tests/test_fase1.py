"""Verificación de la Fase 1. Ejecutar:  PYTHONPATH=src pytest -q tests"""
import numpy as np
import pytest

from traj4d import atmosfera as atm
from traj4d.crucero import (Envolvente, integrar_crucero, masa_final_breguet, perfil_constante)
from traj4d.rendimiento import CongeladoModel, PollSchumannModel


# ----------------------------------------------------------------- Atmósfera
def test_isa_valores_tabulados():
    # ISA: FL350 (10668 m) -> 238.4 hPa, 218.8 K; tropopausa 11 km -> 226.32 hPa, 216.65 K
    assert atm.presion(atm.fl_a_m(350)) == pytest.approx(23842, rel=2e-4)
    assert atm.temperatura(atm.fl_a_m(350)) == pytest.approx(218.81, abs=0.01)
    assert atm.presion(11000.0) == pytest.approx(22632, rel=1e-4)
    assert atm.presion(0.0) == pytest.approx(101325.0)


def test_isa_continuidad_tropopausa():
    assert atm.presion(11000 - 1e-6) == pytest.approx(atm.presion(11000 + 1e-6), rel=1e-9)


# ------------------------------------------------- Integrador frente a Breguet
@pytest.mark.parametrize("viento", [0.0, -40.0, 40.0])
def test_integrador_reproduce_breguet(viento):
    """Con L/D, TSFC y V constantes, el integrador debe reproducir Breguet exactamente,
    también con viento uniforme (V se sustituye por la velocidad sobre el suelo)."""
    mod = CongeladoModel(LD=19.0, TSFC=1.5e-5)
    m0, R, fl, M = 225e3, 5_770e3, 370, 0.85
    res = integrar_crucero(mod, m0, R, perfil_constante(fl, M), viento=viento, ds=10e3)
    V = float(atm.mach_a_tas(M, atm.fl_a_m(fl))) + viento
    m1 = masa_final_breguet(m0, R, V, 19.0, 1.5e-5)
    assert res.m[-1] == pytest.approx(m1, rel=1e-10)


def test_convergencia_rk4_orden4():
    """El error del modelo real debe caer ~16x al dividir el paso a la mitad."""
    mod = PollSchumannModel()
    p = perfil_constante(370, 0.85)
    ref = integrar_crucero(mod, 225e3, 5_770e3, p, ds=1e3).m[-1]
    e1 = abs(integrar_crucero(mod, 225e3, 5_770e3, p, ds=200e3).m[-1] - ref)
    e2 = abs(integrar_crucero(mod, 225e3, 5_770e3, p, ds=100e3).m[-1] - ref)
    assert e1 / e2 > 10.0
    assert e2 < 1.0   # kg


# ------------------------------------- Poll–Schumann frente a la referencia
def test_ps_coincide_con_referencia_pycontrails():
    from pycontrails.models.ps_model import ps_nominal_grid
    from scipy.optimize import brentq
    lv = np.array([300.0, 250.0, 200.0])
    ref = ps_nominal_grid("B789", level=lv)
    mod = PollSchumannModel()
    for i, hpa in enumerate(lv):
        fl = brentq(lambda f: atm.presion(atm.fl_a_m(f)) - hpa * 100, 0, 600)
        e = mod.estado(float(ref.aircraft_mass[i]), ref.attrs["mach_number"], fl)
        assert float(e["FF"]) == pytest.approx(float(ref.fuel_flow[i]), rel=1e-3)


# ----------------------------------------------------------------- Envolvente
def test_envolvente_sentido_fisico():
    env = Envolvente()
    assert env.valido(170e3, 0.85, 390)          # peso medio-bajo a FL390: válido
    assert not env.valido(250e3, 0.85, 430)      # muy pesado a FL430: fuera de envolvente
    assert not env.valido(170e3, 0.85, 450)      # por encima del techo certificado
