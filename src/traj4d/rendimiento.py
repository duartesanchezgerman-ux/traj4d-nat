"""Modelos de prestaciones en crucero (vuelo nivelado estacionario: empuje = resistencia).

Interfaz común:
    m.estado(masa, mach, fl, dT) -> dict con FF [kg/s], D [N], LD, TSFC [kg/(N·s)], TAS [m/s], CL

Implementaciones:
  - OpenAPModel        : OpenAP (TU Delft), polar parabólica + modelo de consumo derivado de Acropole.
  - PollSchumannModel  : Poll & Schumann (Aeronaut. J. 2021, 2024), implementado en pycontrails.
  - CongeladoModel     : L/D y TSFC constantes -> solución exacta de Breguet (solo para verificación).

Tener DOS modelos independientes es deliberado: la diferencia entre ellos es una cota
inferior de la incertidumbre de modelado, y hay que declararla en el artículo.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import atmosfera as atm


class ModeloRendimiento:
    nombre = "base"

    def estado(self, masa, mach, fl, dT=0.0) -> dict:
        raise NotImplementedError

    def ff(self, masa, mach, fl, dT=0.0):
        return self.estado(masa, mach, fl, dT)["FF"]


# --------------------------------------------------------------------------- OpenAP
class OpenAPModel(ModeloRendimiento):
    def __init__(self, tipo: str = "B789"):
        import openap
        self.nombre = f"OpenAP-{tipo}"
        self._drag = openap.Drag(tipo)
        self._fuel = openap.FuelFlow(tipo)
        self._thrust = openap.Thrust(tipo)
        self.S = openap.prop.aircraft(tipo)["wing"]["area"]

    def estado(self, masa, mach, fl, dT=0.0):
        masa, mach, fl = np.broadcast_arrays(*(np.asarray(x, float) for x in (masa, mach, fl)))
        h = atm.fl_a_m(fl)
        tas = atm.mach_a_tas(mach, h, dT)
        tas_kt, alt_ft = tas / atm.KT, fl * 100.0
        D = np.asarray(self._drag.clean(mass=masa, tas=tas_kt, alt=alt_ft, dT=dT), float)
        FF = np.asarray(self._fuel.enroute(mass=masa, tas=tas_kt, alt=alt_ft, dT=dT), float)
        p = atm.presion(h)
        CL = masa * atm.G0 / (0.5 * atm.GAMMA * p * mach**2 * self.S)
        return dict(FF=FF, D=D, LD=masa * atm.G0 / D, TSFC=FF / D, TAS=tas, CL=CL)

    def empuje_max_crucero(self, mach, fl, dT=0.0):
        h = atm.fl_a_m(fl)
        tas_kt = atm.mach_a_tas(mach, h, dT) / atm.KT
        return np.asarray(self._thrust.cruise(tas=tas_kt, alt=np.asarray(fl) * 100.0, dT=dT), float)


# ----------------------------------------------------------------- Poll–Schumann
class PollSchumannModel(ModeloRendimiento):
    """Crucero nivelado, dV/dt = 0, ángulo de trayectoria 0.

    Usa las primitivas de pycontrails.models.ps_model.ps_model en el mismo orden que PSFlight,
    con deterioro de motor por defecto del 2,5 % y Q_fuel = 43,13 MJ/kg.
    """

    def __init__(self, tipo: str = "B789", deterioro: float = 0.025, q_fuel: float = 43.13e6):
        from pycontrails.models.ps_model import ps_model as pm
        from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params
        self.pm = pm
        self.p = load_aircraft_engine_params()[tipo]
        self.nombre = f"PollSchumann-{tipo}"
        self.det, self.q = deterioro, q_fuel
        self.S = self.p.wing_surface_area

    def estado(self, masa, mach, fl, dT=0.0):
        pm, P = self.pm, self.p
        masa, mach, fl = np.broadcast_arrays(*(np.asarray(x, float) for x in (masa, mach, fl)))
        h = atm.fl_a_m(fl)
        p = atm.presion(h)
        T = atm.temperatura(h, dT)
        zero = np.zeros_like(masa)
        cl = pm.lift_coefficient(P.wing_surface_area, masa, p, mach, zero)
        rn = pm.reynolds_number(P.wing_surface_area, mach, T, p)
        cd0 = pm.zero_lift_drag_coefficient(pm.skin_friction_coefficient(rn), P.psi_0)
        e_ls = pm.oswald_efficiency_factor(cd0, P)
        cdw = pm.wave_drag_coefficient(mach, cl, P)
        cd = pm.airframe_drag_coefficient(cd0, cdw, cl, e_ls, P.wing_aspect_ratio)
        F = pm.thrust_force(masa, cl, cd, zero, zero)
        ct = pm.engine_thrust_coefficient(F, mach, p, P.wing_surface_area)
        ct_eta_b = pm.thrust_coefficient_at_max_efficiency(mach, P.m_des, P.c_t_des)
        eta = pm.overall_propulsion_efficiency(mach, ct, ct_eta_b, P, self.det)
        FF = pm.fuel_mass_flow_rate(p, T, mach, ct, eta, P.wing_surface_area, self.q)
        return dict(FF=np.asarray(FF), D=np.asarray(F), LD=cl / cd, TSFC=np.asarray(FF) / np.asarray(F),
                    TAS=atm.mach_a_tas(mach, h, dT), CL=cl, eta=eta)


# -------------------------------------------------------------------- Congelado
@dataclass
class CongeladoModel(ModeloRendimiento):
    LD: float = 18.0
    TSFC: float = 1.45e-5   # kg/(N·s)
    nombre: str = "Congelado"

    def estado(self, masa, mach, fl, dT=0.0):
        masa = np.asarray(masa, float)
        tas = atm.mach_a_tas(mach, atm.fl_a_m(fl), dT)
        D = masa * atm.G0 / self.LD
        return dict(FF=self.TSFC * D, D=D, LD=np.full_like(D, self.LD), TSFC=np.full_like(D, self.TSFC),
                    TAS=tas * np.ones_like(D), CL=np.nan * D)
