"""Comprueba si tu cuenta de OpenSky tiene acceso a la base histórica (Trino).

Requisitos: pip install pyopensky
Credenciales en el fichero de configuración de pyopensky. Este script te dice la ruta exacta.
En macOS suele ser: ~/Library/Application Support/pyopensky/settings.conf

    [default]
    username = TU_USUARIO
    password = TU_CONTRASEÑA

Uso:  python scripts/datos/comprobar_opensky.py
"""
import sys

try:
    from pyopensky.config import opensky_config_dir
    from pyopensky.trino import Trino
except ImportError:
    sys.exit("Falta pyopensky:  pip install pyopensky")

print(f"Fichero de configuración esperado: {opensky_config_dir / 'settings.conf'}")

try:
    trino = Trino()
    # Consulta mínima: vuelos Madrid -> Nueva York JFK de un día concreto.
    df = trino.flightlist("2024-01-15", "2024-01-16",
                          departure_airport="LEMD", arrival_airport="KJFK", limit=10)
except Exception as e:  # noqa: BLE001
    print("\nSIN ACCESO a Trino (o credenciales mal puestas).")
    print(f"Detalle: {type(e).__name__}: {e}")
    print("-> Si las credenciales son correctas, falta que OpenSky apruebe el acceso de investigación.")
    sys.exit(1)

if df is None or df.empty:
    print("\nConexión OK, pero la consulta no devolvió vuelos (revisa fechas o aeropuertos).")
else:
    print(f"\nACCESO OK. {len(df)} vuelos LEMD->KJFK el 2024-01-15:")
    print(df.head(10).to_string())
