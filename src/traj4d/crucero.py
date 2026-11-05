"""Integración del crucero (modelo de masa puntual), envolvente operativa y coste del tiempo.

Ecuación (vuelo nivelado estacionario, con la derrota impuesta por la ruta):

    dm/ds = -FF(m, M, FL, ΔISA) / Vg,     Vg = w_l + sqrt(TAS² − w_t²)

Se integra en la distancia sobre el suelo s. El control (FL, M) se retiene dentro de cada paso,
igual que un piloto que cambia de nivel o de Mach en instantes discretos. Los ascensos en
escalón se tratan como instantáneos; esto subestima el coste de cada escalón en unas decenas de kg.

Coste del tiempo (cost index, CI, en kg/min, convención Airbus). Coste total = combustible + CI·minutos.
El óptimo local ECON minimiza el coste por unidad de distancia sobre el suelo:  J = (FF + CI/60) / Vg.

Breguet (jet, con L/D, TSFC y V constantes):  R = V·(L/D) / (g·TSFC) · ln(m0/m1)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from . import atmosfera as atm
from .rendimiento import ModeloRendimiento


# ------------------------------------------------------------------- Breguet
def masa_final_breguet(m0, rango, V, LD, TSFC):
    return m0 * np.exp(-rango * atm.G0 * TSFC / (V * LD))


def rango_breguet(m0, m1, V, LD, TSFC):
    return V * LD / (atm.G0 * TSFC) * np.log(m0 / m1)


# ---------------------------------------------------------------- Envolvente
@dataclass
class Envolvente:
    """Envolvente de crucero con los límites del modelo Poll–Schumann:
    - masa <= máxima admisible por el CL máximo utilizable (proxy del margen de buffet)
    - empuje disponible >= empuje necesario para 300 ft/min de ascenso residual
    - M <= M máximo a esa altitud (MMO/VMO)
    - FL <= FL máximo certificado (431 para el B789 en Poll–Schumann)
    """
    tipo: str = "B789"
    roc_residual_fpm: float = 300.0
    p: object = field(init=False)

    def __post_init__(self):
        from pycontrails.models.ps_model.ps_aircraft_params import load_aircraft_engine_params
        self.p = load_aircraft_engine_params()[self.tipo]

    def valido(self, masa, mach, fl, dT=0.0):
        from pycontrails.models.ps_model import ps_operational_limits as L
        P = self.p
        masa, mach, fl, dT = np.broadcast_arrays(*(np.asarray(x, float) for x in (masa, mach, fl, dT)))
        h = atm.fl_a_m(fl)
        p = atm.presion(h)
        T = atm.temperatura(h, dT)
        tas = atm.mach_a_tas(mach, h, dT)
        theta = np.degrees(np.arcsin(self.roc_residual_fpm * atm.FT / 60.0 / tas))
        m_max = L.max_allowable_aircraft_mass(p, mach, P.m_des, P.c_l_do, P.wing_surface_area, P.amass_mtow)
        exceso = L.get_excess_thrust_available(mach, T, p, masa.copy(), theta, P)
        m_lim = L.max_mach_number_by_altitude(fl * 100.0, p, P.max_mach_num, P.p_i_max, P.p_inf_co)
        return (masa <= m_max) & (np.nan_to_num(exceso, nan=-1.0) >= 0.0) & (mach <= m_lim) & (fl <= P.fl_max)


# ------------------------------------------------------- Viento local
# Un "viento local" es una función FL -> (w_l, w_t, dT): componentes longitudinal y transversal
# respecto a la derrota y ΔISA, todo en la posición y el instante actuales.
VientoLocal = Callable[[np.ndarray], tuple]


def _viento_constante(wl=0.0, wt=0.0, dT=0.0) -> VientoLocal:
    return lambda fl: (np.full(np.shape(fl), wl, float), np.full(np.shape(fl), wt, float),
                       np.full(np.shape(fl), dT, float))


def _vg(tas, wl, wt):
    return wl + np.sqrt(np.maximum(tas**2 - wt**2, 0.0))


# ------------------------------------------------------- Optimizador local
FLS = np.arange(290, 432, 2)
MACHS = np.round(np.arange(0.78, 0.8801, 0.005), 3)


def optimo_local(modelo: ModeloRendimiento, env: Envolvente | None, masa, viento=0.0, ci=0.0,
                 fls=FLS, machs=MACHS, mach_fijo: float | None = None, fl_fijo: float | None = None):
    """(FL, M) que minimizan J = (FF + CI/60) / Vg.

    viento: escalar (componente longitudinal, en m/s) o VientoLocal. ci: kg/min.
    Devuelve (fl, mach, J [kg/m]).
    """
    vl = viento if callable(viento) else _viento_constante(float(viento))
    if mach_fijo is not None:
        machs = np.array([mach_fijo])
    if fl_fijo is not None:
        fls = np.array([fl_fijo])
    wl, wt, dT = vl(np.asarray(fls, float))
    FL, M = np.meshgrid(fls, machs, indexing="ij")
    WL, WT, DT = (np.broadcast_to(np.asarray(a)[:, None], FL.shape) for a in (wl, wt, dT))
    e = modelo.estado(np.full(FL.shape, float(masa)), M, FL, DT)
    J = (e["FF"] + ci / 60.0) / _vg(e["TAS"], WL, WT)
    if env is not None:
        J = np.where(env.valido(masa, M, FL, DT), J, np.inf)
    i = np.unravel_index(np.argmin(J), J.shape)
    return float(FL[i]), float(M[i]), float(J[i])


# ---------------------------------------------------------------- Perfiles
# Perfil: (s, m, memoria, viento_local) -> (FL, M)
Perfil = Callable[[float, float, dict, VientoLocal], tuple[float, float]]


def perfil_constante(fl, mach) -> Perfil:
    return lambda s, m, mem, vl: (fl, mach)


def perfil_escalones(modelo, env, fl_ini, mach, paso=20, ci=0.0, fl_max=430) -> Perfil:
    """Ascenso en escalones de `paso` centenas de pies cuando el nivel superior tiene menor coste
    J y es válido. Es la operación real en RVSM (normalmente 2000 ft)."""
    def f(s, m, mem, vl):
        fl = mem.setdefault("fl", fl_ini)
        sup = fl + paso
        if sup <= fl_max:
            cand = np.array([fl, sup], float)
            wl, wt, dT = vl(cand)
            e = modelo.estado(np.full(2, m), mach, cand, dT)
            J = (e["FF"] + ci / 60.0) / _vg(e["TAS"], wl, wt)
            ok = env is None or bool(env.valido(m, mach, sup, dT[1]))
            if ok and J[1] < J[0]:
                mem["fl"] = fl = sup
        return fl, mach
    return f


def perfil_optimo(modelo, env, ci=0.0, mach_fijo=None) -> Perfil:
    """Crucero ascendente continuo con (FL, M) óptimos en cada paso. Es la cota ideal."""
    return lambda s, m, mem, vl: optimo_local(modelo, env, m, vl, ci, mach_fijo=mach_fijo)[:2]


def perfil_mach_optimo_fl_escalones(modelo, env, fl_ini, ci=0.0, paso=20, fl_max=430) -> Perfil:
    """Escalones de nivel como en la operación real + Mach ECON en cada paso (lo que hace un FMS)."""
    def f(s, m, mem, vl):
        fl = mem.setdefault("fl", fl_ini)
        sup = fl + paso
        if sup <= fl_max:
            _, _, J0 = optimo_local(modelo, env, m, vl, ci, fl_fijo=fl)
            _, _, J1 = optimo_local(modelo, env, m, vl, ci, fl_fijo=sup)
            if J1 < J0:
                mem["fl"] = fl = sup
        return fl, optimo_local(modelo, env, m, vl, ci, fl_fijo=fl)[1]
    return f


# --------------------------------------------------------------- Resultado
@dataclass
class ResultadoCrucero:
    s: np.ndarray        # m
    m: np.ndarray        # kg
    t: np.ndarray        # s
    fl: np.ndarray
    mach: np.ndarray
    lat: np.ndarray | None = None
    lon: np.ndarray | None = None

    @property
    def combustible(self):
        return self.m[0] - self.m[-1]

    @property
    def tiempo_h(self):
        return self.t[-1] / 3600.0

    def coste(self, ci):
        """Coste total en kg equivalentes de combustible: combustible + CI [kg/min] · minutos."""
        return self.combustible + ci * self.t[-1] / 60.0


# --------------------------------------------------------------- Integradores
def _rk4(modelo, m0, n, h, perfil, viento_en, t0=0.0, geo=None):
    """Núcleo común de RK4. viento_en(s, t) -> VientoLocal. geo(s) -> (lat, lon) o None."""
    s, m, t = 0.0, float(m0), float(t0)
    mem: dict = {}
    S, Mm, T, F, MA, LA, LO = [s], [m], [t], [], [], [], []
    if geo:
        la, lo = geo(s); LA.append(la); LO.append(lo)

    for _ in range(n):
        fl, mach = perfil(s, m, mem, viento_en(s, t))

        def deriv(si, mi, ti):
            wl, wt, dT = (float(np.asarray(x).ravel()[0]) for x in viento_en(si, ti)(np.array([fl])))
            e = modelo.estado(mi, mach, fl, dT)
            vg = float(_vg(float(e["TAS"]), wl, wt))
            return -float(e["FF"]) / vg, 1.0 / vg

        k1 = deriv(s, m, t)
        k2 = deriv(s + h / 2, m + h / 2 * k1[0], t + h / 2 * k1[1])
        k3 = deriv(s + h / 2, m + h / 2 * k2[0], t + h / 2 * k2[1])
        k4 = deriv(s + h, m + h * k3[0], t + h * k3[1])
        m += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        t += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        s += h
        S.append(s); Mm.append(m); T.append(t); F.append(fl); MA.append(mach)
        if geo:
            la, lo = geo(s); LA.append(la); LO.append(lo)

    F.append(F[-1]); MA.append(MA[-1])
    return ResultadoCrucero(np.array(S), np.array(Mm), np.array(T) - t0, np.array(F), np.array(MA),
                            np.array(LA) if geo else None, np.array(LO) if geo else None)


def integrar_crucero(modelo: ModeloRendimiento, m0: float, distancia: float, perfil: Perfil,
                     viento: float = 0.0, dT: float = 0.0, ds: float = 10_000.0) -> ResultadoCrucero:
    """Crucero sin geometría: viento longitudinal uniforme y ΔISA uniforme."""
    n = int(np.ceil(distancia / ds))
    vl = _viento_constante(float(viento), 0.0, dT)
    return _rk4(modelo, m0, n, distancia / n, perfil, lambda s, t: vl)


def integrar_ruta(modelo: ModeloRendimiento, m0: float, ruta, perfil: Perfil, campo, t0: float = 0.0,
                  ds: float = 10_000.0) -> ResultadoCrucero:
    """Crucero a lo largo de una Ruta geográfica, atravesando un Campo meteorológico 4D.

    El viento se muestrea en la posición y el instante de cada evaluación de RK4: es un 4D real.
    """
    n = int(np.ceil(ruta.longitud / ds))

    def viento_en(s, t):
        lat, lon, chi = ruta.punto(s)
        sx, cx = np.sin(chi), np.cos(chi)

        def vl(fls):
            fls = np.atleast_1d(np.asarray(fls, float))
            out = np.array([campo.muestra(t, f, lat, lon) for f in fls])
            u, v, dT = out[:, 0], out[:, 1], out[:, 2]
            return u * sx + v * cx, u * cx - v * sx, dT
        return vl

    return _rk4(modelo, m0, n, ruta.longitud / n, perfil, viento_en, t0=t0,
                geo=lambda s: ruta.punto(s)[:2])
