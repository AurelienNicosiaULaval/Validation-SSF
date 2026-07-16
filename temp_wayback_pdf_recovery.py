#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures as cf
import csv
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import quote

from pypdf import PdfReader

OUT = Path("notes_master_math_france_wayback")
OUT.mkdir(parents=True, exist_ok=True)

WAYBACK_HOST = "web.archive.org"
WAYBACK_IPS = ["207.241.237.3", "207.241.237.7", "207.241.237.10"]
USER_AGENT = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36"

RECORDS = [
    {
        "num": 1,
        "title": "Cours de probabilités",
        "path": "01 - Proba-stat/01.1 Probabilites M1/01 - Cours de probabilites - Yves Coudene (2015).pdf",
        "urls": [
            "https://perso.lpsm.paris/~coudene/probabilites.pdf",
            "http://perso.lpsm.paris/~coudene/probabilites.pdf",
            "https://www.lpsm.paris/~coudene/probabilites.pdf",
            "http://www.proba.jussieu.fr/pageperso/coudene/probabilites.pdf",
            "http://www.proba.jussieu.fr/pageperso/coudene/poly-proba.pdf",
        ],
        "patterns": [
            "perso.lpsm.paris/*coudene*probabilites.pdf",
            "www.proba.jussieu.fr/*coudene*probabilites.pdf",
        ],
    },
    {
        "num": 5,
        "title": "Probabilités approfondies : martingales et chaînes de Markov",
        "path": "01 - Proba-stat/01.2 Processus, Markov et martingales/05 - Probabilites approfondies martingales et chaines de Markov - Thomas Duquesne (2012).pdf",
        "urls": [
            "https://perso.lpsm.paris/~broutinn/teaching/4M011_poly_duquesne.pdf",
            "http://perso.lpsm.paris/~broutinn/teaching/4M011_poly_duquesne.pdf",
            "https://www.lpsm.paris/~broutinn/teaching/4M011_poly_duquesne.pdf",
            "http://www.proba.jussieu.fr/pageperso/broutin/teaching/4M011_poly_duquesne.pdf",
            "http://www.proba.jussieu.fr/pageperso/duquesne/4M011_poly_duquesne.pdf",
        ],
        "patterns": [
            "perso.lpsm.paris/*4M011_poly_duquesne.pdf",
            "www.proba.jussieu.fr/*4M011_poly_duquesne.pdf",
        ],
    },
    {
        "num": 11,
        "title": "Statistique, Partie 2 : approche bayésienne",
        "path": "01 - Proba-stat/01.5 Bayesien et MCMC/11 - Statistique, Partie 2 approche bayesienne - Anna Ben-Hamou; Arnaud Guyader.pdf",
        "urls": [
            "https://perso.lpsm.paris/~aguyader/files/teaching/M1/PolycopiePartie2.pdf",
            "http://perso.lpsm.paris/~aguyader/files/teaching/M1/PolycopiePartie2.pdf",
            "https://www.lpsm.paris/~aguyader/files/teaching/M1/PolycopiePartie2.pdf",
            "http://www.proba.jussieu.fr/pageperso/guyader/files/teaching/M1/PolycopiePartie2.pdf",
            "http://www.proba.jussieu.fr/pageperso/aguyader/files/teaching/M1/PolycopiePartie2.pdf",
        ],
        "patterns": [
            "perso.lpsm.paris/*PolycopiePartie2.pdf",
            "www.proba.jussieu.fr/*PolycopiePartie2.pdf",
        ],
    },
    {
        "num": 15,
        "title": "Calcul stochastique et processus de diffusion",
        "path": "01 - Proba-stat/01.3 Calcul stochastique et diffusions/15 - Calcul stochastique et processus de diffusion - Nicolas Fournier.pdf",
        "urls": [
            "https://perso.lpsm.paris/~nfournier/PolyCS.pdf",
            "http://perso.lpsm.paris/~nfournier/PolyCS.pdf",
            "https://www.lpsm.paris/~nfournier/PolyCS.pdf",
            "http://www.proba.jussieu.fr/pageperso/fournier/PolyCS.pdf",
            "http://www.proba.jussieu.fr/pageperso/nfournier/PolyCS.pdf",
        ],
        "patterns": [
            "perso.lpsm.paris/*PolyCS.pdf",
            "www.proba.jussieu.fr/*PolyCS.pdf",
        ],
    },
    {
        "num": 19,
        "title": "Modélisation et statistique bayésienne computationnelle",
        "path": "01 - Proba-stat/01.5 Bayesien et MCMC/19 - Modelisation et statistique bayesienne computationnelle - Nicolas Bousquet (2026).pdf",
        "urls": [
            "https://perso.lpsm.paris/~bousquet/poly-complet-2026-V1.pdf",
            "http://perso.lpsm.paris/~bousquet/poly-complet-2026-V1.pdf",
            "https://www.lpsm.paris/~bousquet/poly-complet-2026-V1.pdf",
            "https://perso.lpsm.paris/~bousquet/poly-complet-2026.pdf",
            "https://perso.lpsm.paris/~bousquet/cours-complet-2026.pdf",
        ],
        "patterns": [
            "perso.lpsm.paris/*bousquet*2026*.pdf",
            "perso.lpsm.paris/*poly-complet-2026*.pdf",
        ],
    },
    {
        "num": 25,
        "title": "Méthodes de tenseurs pour les problèmes en grande dimension",
        "path": "03 - EDP et calcul scientifique/25 - Methodes de tenseurs pour les problemes en grande dimension (2024).pdf",
        "urls": [
            "https://www.ljll.fr/MathModel/enseignement/cours/TenseursM2_2024.pdf",
            "http://www.ljll.fr/MathModel/enseignement/cours/TenseursM2_2024.pdf",
            "https://www.ljll.math.upmc.fr/MathModel/enseignement/cours/TenseursM2_2024.pdf",
            "http://www.ljll.math.upmc.fr/MathModel/enseignement/cours/TenseursM2_2024.pdf",
        ],
        "patterns": [
            "www.ljll.fr/*TenseursM2_2024.pdf",
            "www.ljll.math.upmc.fr/*TenseursM2_2024.pdf",
        ],
    },
    {
        "num": 33,
        "title": "Contrôle optimal : théorie et applications",
        "path": "02 - Analyse, optimisation et outils/Optimisation et controle/33 - Controle optimal theorie et applications - Emmanuel Trelat.pdf",
        "urls": [
            "https://www.ljll.fr/~trelat/fichiers/livreopt.pdf",
            "http://www.ljll.fr/~trelat/fichiers/livreopt.pdf",
            "https://www.ljll.math.upmc.fr/~trelat/fichiers/livreopt.pdf",
            "http://www.ljll.math.upmc.fr/~trelat/fichiers/livreopt.pdf",
            "http://www.ljll.math.upmc.fr/trelat/fichiers/livreopt.pdf",
        ],
        "patterns": [
            "www.ljll.fr/*trelat*livreopt.pdf",
            "www.ljll.math.upmc.fr/*trelat*livreopt.pdf",
        ],
    },
    {
        "num": 34,
        "title": "Méthodes mathématiques et numériques pour les plasmas",
        "path": "03 - EDP et calcul scientifique/34 - Methodes mathematiques et numeriques pour les plasmas - Bruno Despres (2021).pdf",
        "urls": [
            "https://www.ljll.fr/despres/BD_fichiers/m2_plasma.pdf",
            "http://www.ljll.fr/despres/BD_fichiers/m2_plasma.pdf",
            "https://www.ljll.math.upmc.fr/despres/BD_fichiers/m2_plasma.pdf",
            "http://www.ljll.math.upmc.fr/despres/BD_fichiers/m2_plasma.pdf",
        ],
        "patterns": [
            "www.ljll.fr/*m2_plasma.pdf",
            "www.ljll.math.upmc.fr/*m2_plasma.pdf",
        ],
    },
    {
        "num": 35,
        "title": "Équations aux dérivées partielles elliptiques",
        "path": "03 - EDP et calcul scientifique/35 - Equations aux derivees partielles elliptiques - Herve Le Dret (2010).pdf",
        "urls": [
            "https://www.ljll.fr/ledret/M2Elliptique/chapitre4.pdf",
            "http://www.ljll.fr/ledret/M2Elliptique/chapitre4.pdf",
            "https://www.ljll.math.upmc.fr/ledret/M2Elliptique/chapitre4.pdf",
            "http://www.ljll.math.upmc.fr/ledret/M2Elliptique/chapitre4.pdf",
        ],
        "patterns": [
            "www.ljll.fr/*M2Elliptique/chapitre4.pdf",
            "www.ljll.math.upmc.fr/*M2Elliptique/chapitre4.pdf",
        ],
    },
]


def validate(data: bytes) -> tuple[bytes, int]:
    position = data.find(b"%PDF-")
    if position < 0 or position > 4096:
        raise ValueError(f"signature PDF absente: {data[:30]!r}")
    data = data[position:]
    pages = len(PdfReader(io.BytesIO(data), strict=False).pages)
    if pages < 1:
        raise ValueError("aucune page PDF")
    return data, pages


def curl_bytes(url: str, *, timeout: int = 120) -> bytes:
    errors = []
    for ip in WAYBACK_IPS:
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            temporary = Path(handle.name)
        command = [
            "curl", "-k", "-L", "--fail", "--retry", "3", "--retry-all-errors",
            "--retry-delay", "2", "--connect-timeout", "15", "--max-time", str(timeout),
            "--resolve", f"{WAYBACK_HOST}:443:{ip}",
            "-A", USER_AGENT,
            "-H", "Accept: application/json,application/pdf,*/*;q=0.8",
            "-o", str(temporary),
            url,
        ]
        try:
            process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode == 0 and temporary.exists():
                return temporary.read_bytes()
            errors.append(f"{ip}: {(process.stderr or process.stdout)[-1000:]}")
        finally:
            temporary.unlink(missing_ok=True)
    raise RuntimeError(" | ".join(errors)[-3000:])


def cdx_query(query_url: str) -> list[dict]:
    encoded = quote(query_url, safe="")
    endpoint = (
        f"https://{WAYBACK_HOST}/cdx/search/cdx?url={encoded}"
        "&output=json&filter=statuscode:200"
        "&fl=timestamp,original,mimetype,statuscode,digest"
        "&collapse=digest&limit=100&sort=reverse"
    )
    raw = curl_bytes(endpoint, timeout=90)
    try:
        data = json.loads(raw.decode("utf-8", errors="replace"))
    except Exception as exc:
        raise RuntimeError(f"CDX non JSON: {raw[:200]!r}") from exc
    if not isinstance(data, list) or len(data) < 2:
        return []
    header = data[0]
    return [dict(zip(header, row)) for row in data[1:] if len(row) == len(header)]


def fetch_snapshot(timestamp: str, original: str) -> tuple[bytes, int, str]:
    snapshot = f"https://{WAYBACK_HOST}/web/{timestamp}id_/{original}"
    raw = curl_bytes(snapshot, timeout=180)
    data, pages = validate(raw)
    return data, pages, snapshot


def recover(record: dict) -> dict:
    errors: list[str] = []
    captures: list[dict] = []
    queries = list(record["urls"]) + list(record["patterns"])

    with cf.ThreadPoolExecutor(max_workers=min(8, len(queries))) as executor:
        future_queries = {executor.submit(cdx_query, query): query for query in queries}
        for future in cf.as_completed(future_queries):
            query = future_queries[future]
            try:
                for capture in future.result():
                    capture["query"] = query
                    captures.append(capture)
            except Exception as exc:
                errors.append(f"CDX {query}: {type(exc).__name__}: {exc}")

    unique: dict[tuple[str, str], dict] = {}
    for capture in captures:
        key = (capture.get("timestamp", ""), capture.get("original", ""))
        unique[key] = capture
    captures = list(unique.values())
    captures.sort(
        key=lambda capture: (
            "pdf" in (capture.get("mimetype") or "").lower(),
            capture.get("timestamp", ""),
        ),
        reverse=True,
    )

    for capture in captures[:80]:
        timestamp = capture.get("timestamp", "")
        original = capture.get("original", "")
        if not timestamp or not original:
            continue
        try:
            data, pages, snapshot = fetch_snapshot(timestamp, original)
            destination = OUT / record["path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
            return {
                **record,
                "status": "ok",
                "source_url": snapshot,
                "source_kind": "wayback-forced-ip",
                "pages": pages,
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
                "captures_found": len(captures),
                "detail": "",
            }
        except Exception as exc:
            errors.append(
                f"snapshot {timestamp} {original}: {type(exc).__name__}: {exc}"
            )

    return {
        **record,
        "status": "error",
        "source_url": "",
        "source_kind": "wayback-forced-ip",
        "pages": 0,
        "bytes": 0,
        "sha256": "",
        "captures_found": len(captures),
        "detail": " | ".join(errors)[-20000:],
    }


results: list[dict] = []
with cf.ThreadPoolExecutor(max_workers=5) as executor:
    futures = {executor.submit(recover, record): record for record in RECORDS}
    for future in cf.as_completed(futures):
        record = futures[future]
        try:
            result = future.result()
        except BaseException:
            import traceback
            result = {
                **record,
                "status": "error",
                "source_url": "",
                "source_kind": "wayback-forced-ip",
                "pages": 0,
                "bytes": 0,
                "sha256": "",
                "captures_found": 0,
                "detail": "UNHANDLED " + traceback.format_exc(),
            }
        results.append(result)
        print(
            f"[{result['num']:02d}] {result['status']} captures={result['captures_found']} "
            f"pages={result['pages']} bytes={result['bytes']} {result['title']}",
            flush=True,
        )

results.sort(key=lambda item: item["num"])
fields = [
    "num", "title", "path", "urls", "patterns", "status", "source_url",
    "source_kind", "pages", "bytes", "sha256", "captures_found", "detail",
]
rows = []
for result in results:
    row = dict(result)
    row["urls"] = " | ".join(row["urls"])
    row["patterns"] = " | ".join(row["patterns"])
    rows.append(row)
with (OUT / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
(OUT / "manifest.json").write_text(
    json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8"
)
errors = [result for result in results if result["status"] != "ok"]
print(f"SUCCES={len(results)-len(errors)} ECHECS={len(errors)}", flush=True)
sys.exit(1 if errors else 0)
