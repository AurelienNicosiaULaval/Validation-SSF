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

OUT = Path("notes_master_math_france_savepage")
OUT.mkdir(parents=True, exist_ok=True)
WAYBACK_HOST = "web.archive.org"
WAYBACK_IPS = ["207.241.237.3", "207.241.237.7", "207.241.237.10"]
UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/126 Safari/537.36"

RECORDS = [
    {
        "num": 1,
        "title": "Cours de probabilités",
        "path": "01 - Proba-stat/01.1 Probabilites M1/01 - Cours de probabilites - Yves Coudene (2015).pdf",
        "url": "https://perso.lpsm.paris/~coudene/probabilites.pdf",
    },
    {
        "num": 19,
        "title": "Modélisation et statistique bayésienne computationnelle",
        "path": "01 - Proba-stat/01.5 Bayesien et MCMC/19 - Modelisation et statistique bayesienne computationnelle - Nicolas Bousquet (2026).pdf",
        "url": "https://perso.lpsm.paris/~bousquet/poly-complet-2026-V1.pdf",
    },
    {
        "num": 33,
        "title": "Contrôle optimal : théorie et applications",
        "path": "02 - Analyse, optimisation et outils/Optimisation et controle/33 - Controle optimal theorie et applications - Emmanuel Trelat.pdf",
        "url": "https://www.ljll.fr/trelat/fichiers/livreopt.pdf",
    },
    {
        "num": 34,
        "title": "Méthodes mathématiques et numériques pour les plasmas",
        "path": "03 - EDP et calcul scientifique/34 - Methodes mathematiques et numeriques pour les plasmas - Bruno Despres (2021).pdf",
        "url": "https://www.ljll.fr/despres/BD_fichiers/m2_plasma.pdf",
    },
    {
        "num": 35,
        "title": "Équations aux dérivées partielles elliptiques",
        "path": "03 - EDP et calcul scientifique/35 - Equations aux derivees partielles elliptiques - Herve Le Dret (2010).pdf",
        "url": "https://www.ljll.fr/ledret/M2Elliptique/chapitre4.pdf",
    },
]


def validate(data: bytes) -> tuple[bytes, int]:
    position = data.find(b"%PDF-")
    if position < 0 or position > 4096:
        raise ValueError(f"signature PDF absente: {data[:40]!r}")
    data = data[position:]
    pages = len(PdfReader(io.BytesIO(data), strict=False).pages)
    if pages < 1:
        raise ValueError("aucune page PDF")
    return data, pages


def curl_wayback(url: str, *, method: str = "GET", timeout: int = 180, headers_only: bool = False) -> tuple[bytes, str]:
    errors = []
    for ip in WAYBACK_IPS:
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            body_path = Path(handle.name)
        with tempfile.NamedTemporaryFile(delete=False, mode="w", encoding="utf-8") as handle:
            header_path = Path(handle.name)
        command = [
            "curl", "-k", "-L", "--fail-with-body", "--retry", "2", "--retry-all-errors",
            "--connect-timeout", "15", "--max-time", str(timeout),
            "--resolve", f"{WAYBACK_HOST}:443:{ip}",
            "-A", UA,
            "-D", str(header_path),
            "-o", str(body_path),
        ]
        if method == "POST":
            command += ["-X", "POST"]
        if headers_only:
            command += ["-I"]
        command.append(url)
        try:
            process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            headers = header_path.read_text(encoding="utf-8", errors="replace") if header_path.exists() else ""
            body = body_path.read_bytes() if body_path.exists() else b""
            if process.returncode == 0:
                return body, headers
            errors.append(f"{ip}: {(process.stderr or process.stdout)[-1200:]} headers={headers[-500:]}")
        finally:
            body_path.unlink(missing_ok=True)
            header_path.unlink(missing_ok=True)
    raise RuntimeError(" | ".join(errors)[-5000:])


def save_now(target_url: str) -> str:
    encoded = quote(target_url, safe=":/?=&%~")
    endpoints = [
        f"https://{WAYBACK_HOST}/save/{encoded}",
        f"https://{WAYBACK_HOST}/save/{target_url}",
    ]
    errors = []
    for endpoint in endpoints:
        for method in ("POST", "GET"):
            try:
                body, headers = curl_wayback(endpoint, method=method, timeout=240)
                combined = headers + "\n" + body.decode("utf-8", errors="replace")
                for line in combined.splitlines():
                    lower = line.lower()
                    if lower.startswith("content-location:"):
                        location = line.split(":", 1)[1].strip()
                        if location.startswith("/"):
                            return f"https://{WAYBACK_HOST}{location}"
                        return location
                # API JSON responses sometimes contain archived_snapshots or job_id.
                try:
                    payload = json.loads(body.decode("utf-8", errors="replace"))
                    for key in ("url", "archive_url", "location"):
                        value = payload.get(key) if isinstance(payload, dict) else None
                        if isinstance(value, str) and "/web/" in value:
                            return value
                except Exception:
                    pass
                # A successful submission may not expose the snapshot immediately.
                if "200" in headers or "302" in headers or "job_id" in combined:
                    return "submitted"
            except Exception as exc:
                errors.append(f"{method} {endpoint}: {type(exc).__name__}: {exc}")
    raise RuntimeError(" | ".join(errors)[-6000:])


def cdx_latest(target_url: str) -> list[dict]:
    encoded = quote(target_url, safe="")
    endpoint = (
        f"https://{WAYBACK_HOST}/cdx/search/cdx?url={encoded}"
        "&output=json&filter=statuscode:200"
        "&fl=timestamp,original,mimetype,statuscode,digest"
        "&collapse=digest&limit=20&sort=reverse"
    )
    body, _ = curl_wayback(endpoint, timeout=90)
    data = json.loads(body.decode("utf-8", errors="replace"))
    if not isinstance(data, list) or len(data) < 2:
        return []
    header = data[0]
    return [dict(zip(header, row)) for row in data[1:] if len(row) == len(header)]


def fetch_snapshot(timestamp: str, original: str) -> tuple[bytes, int, str]:
    snapshot = f"https://{WAYBACK_HOST}/web/{timestamp}id_/{original}"
    body, _ = curl_wayback(snapshot, timeout=240)
    data, pages = validate(body)
    return data, pages, snapshot


def recover(record: dict) -> dict:
    errors = []
    submission = ""
    try:
        submission = save_now(record["url"])
    except Exception as exc:
        errors.append(f"save_now: {type(exc).__name__}: {exc}")

    # Laisser à Save Page Now le temps d'achever la capture.
    for wait_seconds in (10, 25, 45, 75):
        time.sleep(wait_seconds)
        try:
            captures = cdx_latest(record["url"])
        except Exception as exc:
            errors.append(f"cdx après {wait_seconds}s: {type(exc).__name__}: {exc}")
            continue
        captures.sort(
            key=lambda capture: (
                "pdf" in (capture.get("mimetype") or "").lower(),
                capture.get("timestamp", ""),
            ),
            reverse=True,
        )
        for capture in captures[:10]:
            try:
                data, pages, snapshot = fetch_snapshot(capture["timestamp"], capture["original"])
                destination = OUT / record["path"]
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(data)
                return {
                    **record,
                    "status": "ok",
                    "submission": submission,
                    "source_url": snapshot,
                    "pages": pages,
                    "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "detail": "",
                }
            except Exception as exc:
                errors.append(
                    f"snapshot {capture.get('timestamp')} {capture.get('original')}: {type(exc).__name__}: {exc}"
                )
    return {
        **record,
        "status": "error",
        "submission": submission,
        "source_url": "",
        "pages": 0,
        "bytes": 0,
        "sha256": "",
        "detail": " | ".join(errors)[-20000:],
    }


results = []
with cf.ThreadPoolExecutor(max_workers=5) as executor:
    futures = {executor.submit(recover, record): record for record in RECORDS}
    for future in cf.as_completed(futures):
        record = futures[future]
        try:
            result = future.result()
        except BaseException:
            import traceback
            result = {
                **record, "status": "error", "submission": "", "source_url": "",
                "pages": 0, "bytes": 0, "sha256": "", "detail": "UNHANDLED " + traceback.format_exc(),
            }
        results.append(result)
        print(
            f"[{result['num']:02d}] {result['status']} pages={result['pages']} bytes={result['bytes']} "
            f"submission={result['submission']} {result['title']}",
            flush=True,
        )

results.sort(key=lambda item: item["num"])
fields = ["num", "title", "path", "url", "status", "submission", "source_url", "pages", "bytes", "sha256", "detail"]
with (OUT / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(results)
(OUT / "manifest.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
errors = [result for result in results if result["status"] != "ok"]
print(f"SUCCES={len(results)-len(errors)} ECHECS={len(errors)}", flush=True)
sys.exit(1 if errors else 0)
