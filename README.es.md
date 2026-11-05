# traj4d — ¿Cuánto se acercan al óptimo de viento las rutas oceánicas del Atlántico Norte oriental?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23986492.svg)](https://doi.org/10.5281/zenodo.23986492)

Código de investigación del artículo *How Close to Wind-Optimal Are Eastern North Atlantic Oceanic Tracks?*
(G. Duarte Sánchez et al., en preparación para el Journal of Open Aviation Science). Con los scripts de este
repositorio se regenera cada cifra, tabla y figura del artículo. El manuscrito se añadirá cuando se publique.

La documentación principal está en inglés: [README.md](README.md). Los comentarios y nombres del código están en español.

## Resumen
- Decodifica el conjunto de datos ADS-C abierto de la OpenSky Network, pasa su hora (hora local de Europa
  central) a UTC y corrige los informes fechados una hora tarde.
- Valida ERA5 frente al viento y la temperatura medidos a bordo y caracteriza un defecto de la dirección del
  viento en los informes del Boeing 787.
- Valida el modelo de prestaciones Poll–Schumann (pycontrails) frente a la calculadora EMEP/EEA (EUROCONTROL).
- Calcula, para cada tramo oceánico observado, la ruta lateral óptima con viento y el exceso de combustible de
  la ruta volada, y lo separa en estructura de puntos de paso y elección de puntos.

## Uso rápido
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
PYTHONPATH=src python -m pytest -q tests      # 28 pruebas
bash scripts/2_revision_completa.sh           # todo el análisis, con los datos ya descargados (ver README.md)
PYTHONPATH=src python scripts/auditoria_a1.py # auditoría independiente de las cifras
```

## Licencia
Código: Apache-2.0 (`LICENSE`, `NOTICE`). Los datos no se redistribuyen; cada fuente tiene sus condiciones (ver README.md). Cita: doi:10.5281/zenodo.23986492 (ver `CITATION.cff`).

## Uso de inteligencia artificial
El autor usó Claude (Anthropic), un modelo de lenguaje, como asistente para desarrollar el código y procesar los
datos. El autor diseñó el estudio, tomó las decisiones metodológicas y comprobó el código y los resultados.
