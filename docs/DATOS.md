# Registro de datos (procedencia e integridad)

Los datos viven en `datos/` y **no entran en el repositorio**, por tamaño y por licencia.
Para reproducir el estudio hay que obtener los mismos archivos y comprobar su sha256:
`sha256sum <archivo>` en Linux, `shasum -a 256 <archivo>` en macOS.

| Archivo | Fuente | Licencia | Tamaño | sha256 |
|---|---|---|---|---|
| `datos/era5/era5_pl_2024-01-15.nc` | Copernicus CDS, ERA5 hourly data on pressure levels. u, v, t, z; 150–350 hPa; 20–70°N, 90°W–10°E; 0,25°; 24 h. Descargado con `scripts/datos/descargar_era5.py` el 1 nov 2026 | Licence to use Copernicus Products | 96.918.590 B | `a5cc4005b37db8fdd6ae8ab503784448a0da21d02d7cb28d01f674f39f03557f` |
| `datos/opensky_adsc/adsc_decoded.txt` | OpenSky Network, "ADS-C Air Traffic Data Collected by the OpenSky Network", Zenodo, doi:10.5281/zenodo.14659997 (versión de enero de 2025) | Según el registro de Zenodo | 505.078.923 B | `8275c2bd2925881c0ef64597f4369d5a3d545481e49a2da16de18b2906a8b199` |
| `datos/opensky_adsc/adsc_file_structure.md` | Ídem (descripción del formato) | Ídem | 11.064 B | `20b564ac62209b7322789027daf7333a59aa3106f25bd5eaa9ae0d961afdf2e8` |
| `datos/opensky_meta/aircraft-database-complete-2025-08.csv` | OpenSky Network, metadatos de aeronaves (s3.opensky-network.org/data-samples/metadata) | Términos de uso de OpenSky | 107.956.517 B | `41309cab2a092afc90304631532b3b9a7ce1e501ba003d9590c2ab84208daa3c` |

## Caracterización inicial (1 nov 2026)
- ADS-C: 472.769 informes de posición en 386 días (jul 2023 – ene 2025). Sin datos en ene 2024; casi nada en dic 2023 y sep 2024.
- Atlántico Norte (Shanwick, Gander Oceanic, Santa Maria Oceanic): 131.666 informes de posición de 2.979 aviones distintos.
- Cruce con la base de aeronaves: se identifican el 99,6 % (2.966) de los aviones. La matrícula del ADS-C coincide con la de la base en 2.873 casos (96,9 %); las diferencias están pendientes de revisar.
- B789: 295 aviones y 20.058 informes en el NAT. A359: 256 aviones y 25.824 informes.

## Citas obligatorias
- OpenSky: Schäfer et al., "Bringing Up OpenSky: A Large-scale ADS-B Sensor Network for Research", IPSN 2014; y el registro de Zenodo del ADS-C.
- ERA5: Hersbach et al. (2020), QJRMS 146, 1999–2049, doi:10.1002/qj.3803.

## Datos derivados (Fase 2, 1 nov 2026)
Se regeneran con los scripts. Los sha256 sirven para comprobar que la regeneración es idéntica.

| Archivo | Script | Contenido | sha256 |
|---|---|---|---|
| `datos/procesado/adsc_posiciones.csv.gz` | `scripts/datos/procesar_adsc.py` | 621.183 mensajes con posición (todas las regiones) | `fe8badc8afc4dcae46b0148998e9c4c5a3f7eec418d5b753e32a54e36faacb4f` |
| `datos/procesado/adsc_atlantico.csv.gz` | ídem | 232.959 posiciones en la caja 25–75°N, 75°W–5°W | `f7e89234de279776ccc7bf092b17df38fbfff06304a95f3968ef151a9d4b1083` |
| `datos/procesado/aeronaves_adsc.csv` | `scripts/datos/extraer_aeronaves.py` | 4.515 aeronaves de la base de OpenSky presentes en el ADS-C | — |
| `datos/procesado/segmentos_nat.csv` | `scripts/f2_cobertura.py` | 7.421 tramos oceánicos de ≥15° de longitud con ≥4 informes (tipos modelables con Poll–Schumann) | — |

Validación del lector: 720.409 mensajes y 621.183 con posición, frente a las 621.183 posiciones que publica el propio dataset (adsc_file_structure.md). 44.110 cabeceras traen el ICAO ID truncado (defecto de origen); se conservan y se recuperan por matrícula cuando la matrícula es única en la base.


## ERA5 de los días de estudio (Fase 2, descargados el 1–2 nov 2026)
Comando: `python scripts/datos/descargar_era5.py --dias 2023-11-01 2023-11-03 2023-11-04 2023-11-05 2023-11-09 2024-12-25 2023-07-08 2023-07-09 2023-07-15 2023-07-16 --area 70 -60 30 0`
Mismas variables y niveles que el archivo de prueba; área 30–70°N, 60°W–0°.

| Archivo | Tamaño | sha256 |
|---|---|---|
| `datos/era5/era5_pl_2023-07-08.nc` | 45.084.882 B | `67b79c7dc02289be1ae12dc287c7b20f344829ee5abf301dcc8f1d15342fc50a` |
| `datos/era5/era5_pl_2023-07-09.nc` | 44.144.850 B | `c1ef78e8a8e0c545d174b05807e41d949028b5ac3e01284d9434f9892ecaf444` |
| `datos/era5/era5_pl_2023-07-15.nc` | 43.873.755 B | `d884c0056422cd2b03c29ace2ee6198704f8487a71985323de80abfc9fcee64d` |
| `datos/era5/era5_pl_2023-07-16.nc` | 43.190.850 B | `09715f5b1edbc223c7400f6a5c1444a6202311eb66af6953fca191599a066257` |
| `datos/era5/era5_pl_2023-11-01.nc` | 42.486.498 B | `750409a394d21446f0fa2a53ba1ca79350a667865fa4b9abc640b61be828da32` |
| `datos/era5/era5_pl_2023-11-03.nc` | 43.048.196 B | `e03e475d3bf18189bb2f75cf608c19ce17a8939463ac65c5054c059239b63263` |
| `datos/era5/era5_pl_2023-11-04.nc` | 42.524.320 B | `42ce112ce9bafd469d057969f3fc9eb12632ac50e4f6f4ef9750a533fd3517af` |
| `datos/era5/era5_pl_2023-11-05.nc` | 42.942.128 B | `0125725954c0037e20da8a3a546c63499ff3f1ce7225e6a73e21873a9facb146` |
| `datos/era5/era5_pl_2023-11-09.nc` | 44.153.964 B | `c20163eb5f506df25c7dd2845d948ac50c8661987600df973c2f7009e491bd61` |
| `datos/era5/era5_pl_2024-12-25.nc` | 43.392.451 B | `a8433a88314a14f03a04275156a35ff02fdf74e11fa1bb0f0dccf091489b02ee` |

## Defecto de origen en el viento medido del ADS-C (Tag 16)
El 1,38 % de los informes con meteo de los días de estudio (78 de 5.633) tiene la dirección del viento girada
entre 110° y 130° con el módulo correcto. 75 de los 78 son de la familia Boeing 787 (el 4,0 % de sus informes).
Se clasifican y se excluyen en `scripts/f2_validacion_era5.py`, informando con y sin ellos.


## Referencias para validar el modelo de prestaciones (Fase 3, 2 nov 2026)
Germán las guardó en `datos/` y en el entorno de Claude están en `datos/referencias/`. Mismo contenido; se comprueba por sha256.

| Archivo | Fuente | Uso | sha256 |
|---|---|---|---|
| `1.A.3.a Aviation  -Annex 1 - Master emissions calculator - 2023 - Protected - v1.5_18_09_2024.xlsx` | EMEP/EEA Air Pollutant Emission Inventory Guidebook 2023, ©EUROCONTROL 2023. Hoja `database`: combustible CCD y LTO por tipo y distancia (NM), calculado con el método AEM de EUROCONTROL (basado en BADA) | Contraste indirecto con BADA | `6671f6baca0f3242d4c7e2138ac6f46146c5a97490d69ee10ee94eda6b517960` |

Nota: una versión anterior del análisis comparaba también con la ICAO CORSIA CERT; se retiró el 5 nov 2026 y este repositorio no contiene nada de ella.

## Datos derivados de la Fase 3 (2 nov 2026)
| Archivo | Script | Contenido |
|---|---|---|
| `datos/procesado/estructura_tramos.csv` | `scripts/f3_estructura.py` | Por tramo: combustible de la ruta volada, del óptimo libre, del óptimo restringido a puntos de 0,5° en los meridianos múltiplos de 10° y de la ruta planificada (Tag 13); latitudes planificadas y óptimas |
| `datos/procesado/validacion_modelo_ps.csv` | `scripts/f3_validacion_modelo.py` | Combustible de viaje Poll–Schumann frente a EEA/EUROCONTROL (desde el 5 nov 2026 sin columnas de la CERT) |
