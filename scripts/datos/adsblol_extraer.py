"""Extrae de los volcados diarios de adsb.lol (licencia ODbL) solo las trazas de los aviones que interesan.

Fase 4 (artículo A2): para cada día de estudio D se necesitan las trazas ADS-B de los aviones con informes ADS-C
en el Atlántico Norte ese día, en los archivos del día D y del día D−1 (los vuelos hacia el este despegan de
Norteamérica la tarde anterior, en hora UTC).

El volcado de cada día es un .tar de 1–3,5 GB (a veces partido en .tar.aa, .tar.ab). Se lee en flujo, sin guardarlo
en disco: solo se escriben los archivos traces/xx/trace_full_<icao24>.json de los aviones pedidos (JSON comprimido
con gzip, tal cual vienen de readsb).

Uso:
    python scripts/datos/adsblol_extraer.py --dia 2023-11-01 --icao24 lista.txt --salida datos/adsblol
La lista tiene un icao24 (hexadecimal, minúsculas) por línea.

Atribución obligatoria (ODbL): "Contains data from adsb.lol, made available under the ODbL".
"""
import argparse
import subprocess
import sys
import tarfile
from pathlib import Path

REPOS = {"2023": "globe_history_2023", "2024": "globe_history_2024", "2025": "globe_history_2025"}


def urls_del_dia(dia: str, listas_preferidas: dict) -> list[list[str]]:
    """Candidatos de URL para un día, en orden de preferencia. Cada candidato es una lista de partes."""
    y, m, d = dia.split("-")
    tag = f"{y}.{m}.{d}"
    cands = []
    for linea in listas_preferidas.get(y, []):
        if f"v{tag}-" in linea:
            cands.append(linea.strip().split(","))
    # Los días de enero aparecen a veces en el repositorio del año anterior; y 2023 no publica la lista preferida.
    for repo in (REPOS.get(y), REPOS.get(str(int(y) - 1))):
        if repo:
            for pod in ("planes-readsb-prod-0", "planes-readsb-staging-0"):
                base = f"https://github.com/adsblol/{repo}/releases/download/v{tag}-{pod}/v{tag}-{pod}.tar"
                cands += [[base], [base + ".aa", base + ".ab"]]
    vistos, unicos = set(), []
    for c in cands:
        if tuple(c) not in vistos:
            vistos.add(tuple(c)); unicos.append(c)
    return unicos


def existe(url: str) -> bool:
    r = subprocess.run(["curl", "-sL", "-r", "0-0", "-o", "/dev/null", "-w", "%{http_code}", "-m", "30", url],
                       capture_output=True, text=True)
    return r.stdout.strip() in ("200", "206")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dia", required=True)
    ap.add_argument("--icao24", required=True, help="archivo con un icao24 por línea")
    ap.add_argument("--salida", default="datos/adsblol")
    ap.add_argument("--listas", default=str(Path.home() / "adsblol"), help="clones de los repositorios globe_history_*")
    a = ap.parse_args()

    pedidos = {l.strip().lower() for l in open(a.icao24) if l.strip()}
    out = Path(a.salida) / a.dia
    out.mkdir(parents=True, exist_ok=True)
    listas = {}
    for y, repo in REPOS.items():
        f = Path(a.listas) / repo / "PREFERRED_RELEASES.txt"
        if f.exists():
            listas[y] = f.read_text().splitlines()
    partes = next((c for c in urls_del_dia(a.dia, listas) if all(existe(u) for u in c)), None)
    if partes is None:
        print(f"{a.dia}: no se encuentra el volcado en adsb.lol", file=sys.stderr)
        sys.exit(2)
    print(f"{a.dia}: {len(partes)} parte(s): {partes[0].rsplit('/', 1)[-1]}", flush=True)

    # curl de las partes, una detrás de otra, en un único flujo
    cmd = " && ".join(f"curl -sL --retry 5 --retry-delay 10 --speed-limit 20000 --speed-time 180 '{u}'" for u in partes)
    proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, bufsize=1 << 20)
    n_escritos, n_total = 0, 0
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tar:
        for m in tar:
            n_total += 1
            nombre = m.name.rsplit("/", 1)[-1]
            if not (m.isfile() and nombre.startswith("trace_full_")):
                continue
            hexid = nombre[len("trace_full_"):].split(".")[0].lower()
            if hexid in pedidos:
                datos = tar.extractfile(m).read()
                (out / f"{hexid}.json.gz").write_bytes(datos)
                n_escritos += 1
    proc.wait()
    if proc.returncode != 0:
        print(f"{a.dia}: curl terminó con código {proc.returncode}", file=sys.stderr)
        sys.exit(3)
    faltan = len(pedidos) - n_escritos
    print(f"{a.dia}: {n_escritos} trazas guardadas de {len(pedidos)} pedidas (faltan {faltan}); {n_total} entradas leídas",
          flush=True)
    (out / "_COMPLETO").write_text(f"{n_escritos}/{len(pedidos)}\n")


if __name__ == "__main__":
    main()
