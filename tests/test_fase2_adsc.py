"""Verificación del lector ADS-C con mensajes reales copiados del volcado de OpenSky."""
import io

import pytest

from traj4d.adsc import Estadisticas, leer_mensajes

MUESTRA = """Registration: N790AN  ICAO ID: aab812    ATSU Address: YQXE2YA - Gander Oceanic
Channel Frequency: 1526377000
/YQXE2YA.ADS.N790AN07248C0740BE894706D11D0C041331E798200D238E3F1C71C94707C0222226E38E49470010553D3E2848F8
Tag 07 Basic report: 248C0740BE894706D11D
    Latitude: 51.39404296875
    Longitude: -33.61919403076172
    Altitude: 38000.0 ft
    Timestamp: 2023-07-07 15:07:16.250
    Position accuracy: <0.25 nm
    NAV redundancy: OK
    TCAS: OK
Tag 12 Flight ID: 041331E79820
    Flight ID: AAL199
Tag 13 Predicted Route: 238E3F1C71C94707C0222226E38E494700
  Next waypoint:
    Lat: 49.99998092651367
    Lon: -40.00001907348633
    Alt: 38000.0 ft
    ETA: 2023-07-07 15:40:20.250
  Next+1 waypoint:
    Lat: 47.99995422363281
    Lon: -49.99998092651367
    Alt: 38000.0 ft
CRC: 4AB2

Registration: N27015  ICAO ID: a2a915    ATSU Address: YQXE2YA - Gander Oceanic
Channel Frequency: 1526070000
/YQXE2YA.ADS.N2701514288A9F554B488B89639D0D271F571C63C88B89E325B32EE37F088B805CE6
Tag 03 Acknowledgement: 01
    Contract Number: 1
CRC: 6DD7

Registration: EC-NVX  ICAO ID: 3451c6    ATSU Address: PIKCPYA - Shanwick
Channel Frequency: 1526377000
/PIKCPYA.ADS.EC-NVX0000
Tag 20 Waypoint Change Event: 2AA9FFAA9A884D07CE1D
    Latitude: 52.0
    Longitude: -20.0
    Altitude: 37004.0 ft
    Timestamp: 2024-12-01 10:00:00.000
    Position accuracy: <0.05 nm
    NAV redundancy: OK
    TCAS: OK
Tag 14 Earth Reference data: 5CA0E4400C
    True track: 260.5078125 deg
    Ground speed: 456.5 kt
    Vertical speed: 48.0 ft/min
Tag 15 Air Reference Data: 0CD9A28004
    True heading: 262.123046875 deg
    Mach speed: 0.837
    Vertical speed: 16.0 ft/min
Tag 16 Meteorological data: 4AB93E68
    Wind speed: 74.5 kt
    True wind direction: 281.328125 deg
    Temperature: -51.0 C
CRC: 6A4F
"""


def test_icao_truncado_se_conserva_y_se_marca():
    txt = MUESTRA.replace("ICAO ID: aab812", "ICAO ID: 0a48")
    est = Estadisticas()
    a = list(leer_mensajes(io.StringIO(txt), est))[0]
    assert a["icao24"] is None and a["icao24_bruto"] == "0a48" and a["matricula"] == "N790AN"
    assert est["icao24_truncado"] == 1 and "cabeceras_ilegibles" not in est


def test_lee_mensajes_y_campos():
    est = Estadisticas()
    filas = list(leer_mensajes(io.StringIO(MUESTRA), est))
    assert est["mensajes"] == 3 and est["mensajes_con_posicion"] == 2 and est["mensajes_sin_posicion"] == 1
    a, b = filas
    assert a["icao24"] == "aab812" and a["matricula"] == "N790AN" and a["centro"] == "Gander Oceanic"
    assert a["tipo_informe"] == "07" and a["vuelo"] == "AAL199" and a["frecuencia_hz"] == 1526377000
    assert (a["lat"], a["lon"], a["alt_ft"]) == pytest.approx((51.39404296875, -33.61919403076172, 38000.0))
    assert a["t"] == pytest.approx(1688742436.25 - 7200)   # 15:07:16.250 CEST = 13:07:16.250 UTC
    assert (a["wp1_lat"], a["wp1_lon"], a["wp2_lon"]) == pytest.approx((49.99998092651367, -40.00001907348633, -49.99998092651367))
    assert a["wp1_eta_t"] - a["t"] == pytest.approx(33 * 60 + 4)
    assert a["mach"] is None and a["viento_kt"] is None
    assert b["tipo_informe"] == "20" and b["centro"] == "Shanwick"
    assert (b["gs_kt"], b["derrota_deg"], b["mach"], b["rumbo_deg"]) == pytest.approx((456.5, 260.5078125, 0.837, 262.123046875))
    assert (b["viento_kt"], b["viento_dir_deg"], b["temp_c"]) == pytest.approx((74.5, 281.328125, -51.0))
    assert b["vuelo"] is None


def test_hora_del_volcado_es_local_de_europa_central():
    """El Timestamp del volcado viene en CET/CEST, no en UTC (comprobado con ADS-B el 3-11-2026)."""
    import pytest
    from datetime import datetime, timezone
    from traj4d.adsc import _ts
    utc = lambda x: datetime.fromtimestamp(_ts(x), timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    assert utc("2023-07-07 15:07:11.000") == "2023-07-07 13:07:11"     # verano: UTC+2
    assert utc("2023-11-01 06:52:46.000") == "2023-11-01 05:52:46"     # invierno: UTC+1
    assert utc("2024-01-01 00:30:00.000") == "2023-12-31 23:30:00"     # cambia de día
    for x in ("2024-03-31 02:30:00.000", "2023-10-29 02:30:00.000"):   # inexistente / ambigua
        with pytest.raises(ValueError):
            _ts(x)


def test_corrige_informe_fechado_una_hora_tarde():
    """Informe generado a las 10:59 y fechado 11:59 por la vuelta de hora del receptor."""
    import numpy as np
    import pandas as pd
    from traj4d.vuelos import corregir_vuelta_de_hora
    t0 = 1_700_000_000 - 1_700_000_000 % 3600          # comienzo de una hora
    t = t0 + np.arange(0, 4 * 3600, 600.0)               # cada 10 min, 4 h, a 0,1° de longitud por minuto
    lon = -40 + (t - t0) / 60 * 0.1
    df = pd.DataFrame(dict(t=t, lat=50.0, lon=lon, icao24="abc123", matricula="X", wp1_eta_t=t + 900))
    malo = pd.DataFrame(dict(t=[t0 + 3599 + 3600], lat=[50.0], lon=[-40 + 59.98 * 0.1], icao24=["abc123"],
                             matricula=["X"], wp1_eta_t=[t0 + 3599 + 3600 + 900]))
    out = corregir_vuelta_de_hora(pd.concat([df, malo], ignore_index=True))
    assert out.t_corregida_1h.sum() == 1
    assert out.loc[out.t_corregida_1h, "t"].iat[0] == t0 + 3599
    assert out.loc[out.t_corregida_1h, "wp1_eta_t"].iat[0] == t0 + 3599 + 900


def test_corrige_vuelta_de_hora_en_el_ultimo_informe():
    """El último informe de un vuelo, fechado una hora tarde, se detecta por extrapolación."""
    import numpy as np
    import pandas as pd
    from traj4d.vuelos import corregir_vuelta_de_hora
    t0 = 1_700_000_000 - 1_700_000_000 % 3600
    t = t0 + np.arange(0, 3 * 3600, 600.0)                 # hasta t0 + 2:50
    lon = -40 + (t - t0) / 60 * 0.1
    ultimo_real = t0 + 2 * 3600 + 3540                      # generado a las 2:59 ...
    df = pd.DataFrame(dict(t=np.r_[t, ultimo_real + 3600], lat=50.0,   # ... y fechado a las 3:59
                           lon=np.r_[lon, -40 + (ultimo_real - t0) / 60 * 0.1], icao24="abc123", matricula="X"))
    out = corregir_vuelta_de_hora(df)
    assert out.t_corregida_1h.sum() == 1 and out.t.max() == ultimo_real


def test_union_cruce_despegue_y_aterrizaje():
    """Unión ADS-C + ADS-B sintética: huecos y extremos del vuelo."""
    import numpy as np
    import pandas as pd
    from traj4d.cruces import unir_cruce
    t = np.arange(0, 3600 * 2, 60.0)                       # 2 h de ADS-B antes del océano: despega y sube
    antes = pd.DataFrame(dict(t=t, lat=50.0, lon=-5 - t / 3600 * 4, alt_ft=np.minimum(t / 60 * 2000, 37000.0)))
    antes.loc[0, "alt_ft"] = np.nan
    despues = pd.DataFrame(dict(t=6 * 3600 + t, lat=45.0, lon=-50 - t / 3600 * 4,
                                alt_ft=np.maximum(37000.0 - t / 60 * 400, 500.0)))
    adsb = pd.concat([antes, despues], ignore_index=True)
    adsb["en_tierra"] = False
    adsb.loc[0, "en_tierra"] = True
    adsc = pd.DataFrame(dict(t=[3 * 3600.0, 5 * 3600.0], lat=[50.0, 46.0], lon=[-20.0, -40.0], alt_ft=37000.0))
    r = unir_cruce(adsc, adsb)
    assert r["t_despegue"] == 0.0 and abs(r["hueco_entrada_h"] - (3 * 3600 - t[-1]) / 3600) < 1e-9
    assert r["hueco_salida_h"] == 1.0 and np.isfinite(r["t_aterrizaje"]) and r["alt_arr"] <= 3000
