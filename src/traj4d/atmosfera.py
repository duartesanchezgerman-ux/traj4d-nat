"""Atmósfera estándar internacional (ISA) con desviación de temperatura ΔISA.

Convención: el nivel de vuelo (FL) define la PRESIÓN (altitud presión). La desviación
ΔISA modifica solo la temperatura, y con ella la velocidad del sonido y la TAS para un
Mach dado. Así es como vuela un avión real: altímetro en 1013,25 hPa por encima de la
altitud de transición.
Válida hasta 20 km (troposfera + tropopausa isoterma).
"""
from __future__ import annotations

import numpy as np

G0 = 9.80665          # m/s²
R_AIR = 287.05287     # J/(kg·K)
GAMMA = 1.4
T0 = 288.15           # K
P0 = 101325.0         # Pa
LAPSE = -0.0065       # K/m
H_TROPO = 11000.0     # m
T_TROPO = T0 + LAPSE * H_TROPO                                  # 216.65 K
P_TROPO = P0 * (T_TROPO / T0) ** (-G0 / (LAPSE * R_AIR))        # ≈ 22632 Pa

FT = 0.3048
KT = 0.514444


def fl_a_m(fl):
    """Nivel de vuelo → altitud presión en metros."""
    return np.asarray(fl, dtype=float) * 100.0 * FT


def presion(h):
    """Presión ISA [Pa] a la altitud presión h [m]."""
    h = np.asarray(h, dtype=float)
    p_trop = P0 * (1 + LAPSE * h / T0) ** (-G0 / (LAPSE * R_AIR))
    p_strat = P_TROPO * np.exp(-G0 * (h - H_TROPO) / (R_AIR * T_TROPO))
    return np.where(h <= H_TROPO, p_trop, p_strat)


def temperatura(h, dT=0.0):
    """Temperatura [K] = T_ISA(h) + ΔISA."""
    h = np.asarray(h, dtype=float)
    return np.where(h <= H_TROPO, T0 + LAPSE * h, T_TROPO) + dT


def vel_sonido(T):
    return np.sqrt(GAMMA * R_AIR * np.asarray(T, dtype=float))


def mach_a_tas(mach, h, dT=0.0):
    """TAS [m/s]."""
    return np.asarray(mach) * vel_sonido(temperatura(h, dT))
