"""Verificación del optimizador lateral (programación dinámica)."""
import numpy as np
import pytest

from traj4d.crucero import integrar_ruta, perfil_constante
from traj4d.optimizador import destino, optimizar_lateral
from traj4d.rendimiento import PollSchumannModel
from traj4d.ruta import AEROPUERTOS, Ruta, distancia, rumbo_inicial
from traj4d.viento import CampoUniforme, ChorroSintetico

PS = PollSchumannModel()


def test_destino_inverso_de_distancia_y_rumbo():
    la, lo = destino(45.0, -30.0, np.radians(70.0), 500e3)
    assert distancia(45.0, -30.0, la, lo) == pytest.approx(500e3, rel=1e-9)
    assert np.degrees(rumbo_inicial(45.0, -30.0, la, lo)) == pytest.approx(70.0, abs=1e-9)


def test_sin_viento_da_la_ortodromica():
    r = optimizar_lateral(PS, CampoUniforme(), 40.47, -3.56, 40.64, -73.78, 350, 0.85, 225e3,
                          n_etapas=12, K=6, delta_km=50.0)
    assert np.all(r.desplazamientos_km == 0.0)


def test_viento_zonal_uniforme_en_ecuador_da_la_ortodromica():
    """Ruta sobre el ecuador con viento zonal uniforme: por simetría N–S, el óptimo es el ecuador."""
    r = optimizar_lateral(PS, CampoUniforme(u=-40.0), 0.0, -10.0, 0.0, -50.0, 350, 0.85, 225e3,
                          n_etapas=10, K=6, delta_km=50.0)
    assert np.all(r.desplazamientos_km == 0.0)


def test_chorro_sintetico_mejora_las_alternativas_manuales():
    """MAD→JFK hacia el oeste contra un chorro de 75 m/s a 45°N: la DP debe igualar o mejorar la mejor
    alternativa manual de la Fase 1b (desvíos por 37N, 45N o 53N a 40W). Todo se evalúa con RK4 a FL350, M0,85."""
    jet = ChorroSintetico()
    (a, b), (c, d) = AEROPUERTOS["LEMD"], AEROPUERTOS["KJFK"]
    prof = perfil_constante(350, 0.85)
    alternativas = [Ruta.ortodromica("LEMD", "KJFK", n=40)] + [Ruta([a, la, c], [b, -40.0, d]) for la in (37.0, 45.0, 53.0)]
    mejor_manual = min(integrar_ruta(PS, 225e3, rr, prof, jet, ds=20e3).combustible for rr in alternativas)
    gc = integrar_ruta(PS, 225e3, alternativas[0], prof, jet, ds=20e3).combustible
    dp = optimizar_lateral(PS, jet, a, b, c, d, 350, 0.85, 225e3, n_etapas=24, K=25, delta_km=40.0)
    f_dp = integrar_ruta(PS, 225e3, dp.ruta, prof, jet, ds=20e3).combustible
    assert f_dp <= mejor_manual * 1.001
    assert f_dp < gc * 0.98            # rodear el chorro ahorra claramente frente a la ortodrómica


def test_ruta_volada_conserva_extremos():
    """El último informe se conserva aunque esté a menos de 5 km del anterior (auditoría 3-11-2026)."""
    import numpy as np
    from traj4d.ruta import puntos_distintos
    lat = np.array([50.0, 50.0, 50.0, 50.0, 50.0])
    lon = np.array([-30.0, -29.0, -28.0, -27.99, -27.98])     # los dos últimos, a < 1 km
    idx = puntos_distintos(lat, lon)
    assert idx[0] == 0 and idx[-1] == 4 and list(idx) == [0, 1, 4]
    assert list(puntos_distintos([50.0, 50.0], [-30.0, -29.9999])) == [0, 1]
