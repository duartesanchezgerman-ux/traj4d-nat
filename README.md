# traj4d — How close to wind-optimal are eastern North Atlantic oceanic tracks?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23986492.svg)](https://doi.org/10.5281/zenodo.23986492)

Research code for the article *How Close to Wind-Optimal Are Eastern North Atlantic Oceanic Tracks?*
(G. Duarte Sánchez et al., in preparation for the Journal of Open Aviation Science). Every number, table and
figure of the article can be regenerated with the scripts in this repository. The manuscript itself will be added
when it is published.

What the code does:
- decodes the open ADS-C dataset of the OpenSky Network (position, intent, Mach, wind and temperature reports),
  converts its time stamps from Central European local time to UTC (checked against ADS-B from adsb.lol) and
  corrects reports dated one hour late at the change of hour;
- validates ERA5 against the wind and temperature measured on board;
- validates the open Poll–Schumann performance model (via pycontrails) against the EMEP/EEA (EUROCONTROL)
  fuel calculator;
- characterises a wind-direction defect in Boeing 787 ADS-C reports (direction rotated by about 125° in
  southerly winds);
- computes, for each observed oceanic segment, the wind-optimal lateral path by dynamic programming and the
  fuel excess of the flown path, and splits that excess into waypoint structure and waypoint choice.

## Installation
Python ≥ 3.11.
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
PYTHONPATH=src python -m pytest -q tests          # 28 tests
```

## Data (not included)
No data are redistributed here. Obtain the files below and check them against the SHA-256 checksums in
`docs/DATOS.md` and `docs/ERA5_SHA256.txt`.

| Data | Source | Licence / terms |
|---|---|---|
| ADS-C reports | OpenSky Network, Zenodo, doi:10.5281/zenodo.14659997 | See the Zenodo record (CC BY 4.0 for the first version) |
| Aircraft database (2025-08) | OpenSky Network | OpenSky terms of use (non-commercial research) |
| ERA5 hourly data on pressure levels | Copernicus Climate Data Store, doi:10.24381/cds.bd0915c6 | CC BY 4.0 |
| EMEP/EEA master emissions calculator 2023 | European Environment Agency / EUROCONTROL | EEA re-use policy |
| ADS-B traces (time-stamp check only) | adsb.lol, `globe_history` releases on GitHub | ODbL 1.0 |

```bash
python scripts/datos/procesar_adsc.py datos/opensky_adsc/adsc_decoded.txt --salida datos/procesado
python scripts/datos/extraer_aeronaves.py datos/opensky_meta/aircraft-database-complete-2025-08.csv datos/procesado/adsc_posiciones.csv.gz datos/procesado/aeronaves_adsc.csv
bash scripts/1_descargar_era5.sh                  # ERA5 for the study days in docs/dias_estudio_era5.txt (needs ~/.cdsapirc)
```

## Reproducing the article
```bash
bash scripts/2_revision_completa.sh               # checks data, runs tests, all analyses and figures (1–3 h)
```

| Script | Article content |
|---|---|
| `scripts/f2_cobertura.py` | Coverage, segments, wind-triangle closure (Section 4.1) |
| `scripts/f2_hora_adsc.py` | Time stamps of the ADS-C export against ADS-B (Section 3.1); needs adsb.lol traces from `scripts/datos/adsblol_extraer.py` |
| `scripts/f2_validacion_era5.py` | ERA5 against in-situ wind and temperature (Table 1, Figure 2) |
| `scripts/f3_direccion_b787.py` | Boeing 787 wind-direction defect by wind direction, type, operator and month (Section 4.1) |
| `scripts/f3_mach.py` | Observed Mach number and implied cost index (Figure 3) |
| `scripts/f3_validacion_modelo.py` | Model validation against the EMEP/EEA (EUROCONTROL) calculator (Table 2, Figure 4) |
| `scripts/f2_comparacion.py` | Flown vs wind-optimal paths (Table 3, Figures 5–6) |
| `scripts/f3_estadistica.py` | Day-block bootstrap confidence intervals |
| `scripts/f3_estructura.py` | Waypoint structure vs choice (Table 4) |
| `scripts/f3_figuras.py --idioma en` | All figures (written to `articulo/en/figures/`) |
| `scripts/auditoria_a1.py` | Independent audit: re-computes each step by other means and checks the numbers quoted in the article |

Comments and identifiers in the code are in Spanish; the article is in English.

## Licence and citation
Code: Apache-2.0 (see `LICENSE` and `NOTICE`). Please cite the article and this software (doi:10.5281/zenodo.23986492; see `CITATION.cff`).
This work contains data from adsb.lol, made available under the Open Database License (ODbL).
ERA5 results contain modified Copernicus Climate Change Service information; neither the European Commission
nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains.

## Use of AI
The author used Claude (Anthropic), a large language model, as an assistant to develop the code and process
the data. The author designed the study, took the methodological decisions, and checked the code and the results.
