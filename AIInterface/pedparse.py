#!/usr/bin/env python3
"""
Pedigree Extractor (Azure OpenAI vision-enabled)

Accepts a local pedigree figure (PNG/JPG/JPEG/WebP/PDF) and asks Azure OpenAI to
extract a structured table of individuals with parentage, obeying the user's
founder rule:

Founder detection logic (STRICT): Only assign parents to an individual if a
continuous line connects them to a horizontal mating line from the generation
above. Individuals without such a connection—regardless of row placement—must be
set as founders with PaternalID = 0 and MaternalID = 0. Parentage must NOT be
inferred based on spatial alignment, proximity, or generation row alone.

Output: JSON (to stdout) and optional CSV file (--csv).

Requirements:
- Python 3.9+
- azure-openai SDK (pip install openai>=1.30.0) or `pip install azure-ai-inference` if using that client.
- For PDFs: `pdf2image` and poppler (optional). If unavailable, use --page-image to provide a pre-rendered image.

Environment:
- AZURE_OPENAI_API_KEY
- AZURE_OPENAI_ENDPOINT (e.g., https://YOUR-RESOURCE.openai.azure.com)
- AZURE_OPENAI_API_VERSION (default: 2024-02-15-preview)
- AZURE_OPENAI_DEPLOYMENT (the model deployment name, e.g., gpt-4.1 or gpt-4o)

Usage:
  python pedigree_extractor.py path/to/pedigree.png --csv pedigree.csv
  python pedigree_extractor.py path/to/pedigree.pdf --page 1 --csv pedigree.csv

Notes:
- The model is prompted to return STRICT JSON. We validate and print helpful
  error messages if parsing fails.
- The script also echoes any model warnings/notes.
"""
from __future__ import annotations
import argparse
import base64
import io
import json
import mimetypes
import os
import sys
from dataclasses import dataclass, asdict
from typing import List, Optional, Any

# Azure OpenAI SDK
try:
    from openai import AzureOpenAI  # pip install openai
except Exception as e:  # pragma: no cover
    AzureOpenAI = None

# Optional PDF support
try:
    from pdf2image import convert_from_path  # requires poppler
except Exception:
    convert_from_path = None

SUPPORTED_IMAGE_TYPES = {"image/png", "image/jpeg", "image/jpg", "image/webp"}
SUPPORTED_PDF_TYPES = {"application/pdf"}

SYSTEM_PROMPT = (
    "You are a meticulous clinical genetics assistant that extracts structured "
    "pedigree information from images of pedigree diagrams.\n\n"
    "Founder detection logic (STRICT) — APPLY EXACTLY: Only assign parents to an "
    "individual if a continuous line connects them to a horizontal mating line "
    "from the generation above. Individuals without such a connection—regardless "
    "of row placement—MUST be treated as founders with PaternalID = 0 and "
    "MaternalID = 0. Do NOT infer parentage based on spatial alignment, "
    "proximity, or generation row alone.\n\n"
    "Return STRICT JSON only (no prose) with this schema: {\n"
    "  \"individuals\": [\n"
    "    {\n"
    "      \"IndividualID\": string,            // unique identifier you assign (e.g., I-1, II-2, etc.)\n"
    "      \"Sex\": string,                     // Male/Female/Unknown\n"
    "      \"Affected\": string,                // Affected/Unaffected/Unknown\n"
    "      \"PaternalID\": string,              // '0' if founder or unknown\n"
    "      \"MaternalID\": string,              // '0' if founder or unknown\n"
    "      \"Notes\": string                    // optional clarifications\n"
    "    }\n"
    "  ],\n"
    "  \"warnings\": [string]                     // any ambiguities or quality issues\n"
    "}\n"
)

USER_HINT = (
    "Extract the table, obeying the strict founder rule.\n"
    "If any mating line is ambiguous or broken, set both parents to '0' and add a warning.\n"
    "Prefer roman numerals by generation if obvious; else assign stable IDs.\n"
)

@dataclass
class Person:
    IndividualID: str
    Sex: str
    Affected: str
    PaternalID: str
    MaternalID: str
    Notes: str = ""

@dataclass
class Extraction:
    individuals: List[Person]
    warnings: List[str]


def _b64_data_url(img_bytes: bytes, mime: str) -> str:
    enc = base64.b64encode(img_bytes).decode("ascii")
    return f"data:{mime};base64,{enc}"


def _load_image_from_pdf(pdf_path: str, page: int) -> bytes:
    if convert_from_path is None:
        raise RuntimeError(
            "PDF support requires pdf2image and poppler. Install pdf2image or "
            "convert the desired page to PNG/JPG and pass that instead with --page-image."
        )
    # convert single page to image (PNG) in-memory
    images = convert_from_path(pdf_path, first_page=page, last_page=page, fmt="png")
    if not images:
        raise RuntimeError(f"No images rendered from {pdf_path} page {page}")
    buf = io.BytesIO()
    images[0].save(buf, format="PNG")
    return buf.getvalue()


def _read_local_bytes(path: str) -> tuple[bytes, str]:
    guessed = mimetypes.guess_type(path)[0] or "application/octet-stream"
    with open(path, "rb") as f:
        data = f.read()
    return data, guessed


def build_client(endpoint: str | None = None, api_key: str | None = None, api_version: str | None = None) -> AzureOpenAI:
    """Build an AzureOpenAI client using the same scheme as prototype.py (key-based)."""
    from openai import AzureOpenAI
    endpoint = (endpoint or os.getenv("AZURE_OPENAI_ENDPOINT") or "https://oric-rag-demo-resource.openai.azure.com").rstrip("/")
    api_key = api_key or os.getenv("AZURE_OPENAI_API_KEY")
    api_version = api_version or os.getenv("AZURE_OPENAI_API_VERSION", "2024-02-15-preview")
    if not api_key:
        raise RuntimeError("AZURE_OPENAI_API_KEY not set. Export it or pass --api-key.")
    masked = (api_key[:4] + "…" + api_key[-4:]) if isinstance(api_key, str) and len(api_key) > 8 else ("set" if api_key else "missing")
    print("Azure OpenAI config:")
    print(f"  endpoint    = {endpoint}")
    print(f"  api_version = {api_version}")
    print(f"  api_key     = {masked}")
    # Also mirror prototype.py's global OpenAI settings for compatibility
    import openai as _openai
    _openai.api_type = "azure"
    _openai.api_base = endpoint + "/"
    _openai.api_version = api_version
    _openai.api_key = api_key

    client = AzureOpenAI(azure_endpoint=endpoint, api_key=api_key, api_version=api_version)
    return client


def call_model_with_image(client: AzureOpenAI, deployment: str, image_bytes: bytes, mime: str) -> dict:
    data_url = _b64_data_url(image_bytes, mime)
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {
            "role": "user",
            "content": [
                {"type": "text", "text": USER_HINT},
                {"type": "image_url", "image_url": {"url": data_url}},
            ],
        },
    ]
    resp = client.chat.completions.create(
        model=deployment,
        messages=messages,
        temperature=0,
        response_format={"type": "json_object"}
    )
    text = resp.choices[0].message.content
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as e:
        # Try to salvage JSON if the model wrapped it in code fences
        cleaned = text.strip()
        if cleaned.startswith("```"):
            cleaned = cleaned.strip("`\n ")
            if cleaned.startswith("json"):
                cleaned = cleaned[4:].lstrip("\n")
        try:
            payload = json.loads(cleaned)
        except Exception:
            raise RuntimeError(f"Model did not return valid JSON. Raw output was:\n{text}")
    return payload


def validate_payload(payload: dict) -> Extraction:
    if not isinstance(payload, dict) or "individuals" not in payload:
        raise RuntimeError("Response missing 'individuals' key.")
    individuals_raw = payload.get("individuals", [])
    warnings = payload.get("warnings", [])
    people: List[Person] = []
    for row in individuals_raw:
        try:
            people.append(
                Person(
                    IndividualID=str(row.get("IndividualID", "")),
                    Sex=str(row.get("Sex", "Unknown")),
                    Affected=str(row.get("Affected", "Unknown")),
                    PaternalID=str(row.get("PaternalID", "0")),
                    MaternalID=str(row.get("MaternalID", "0")),
                    Notes=str(row.get("Notes", "")),
                )
            )
        except Exception as e:
            raise RuntimeError(f"Invalid person row: {row!r}\nError: {e}")
    # Enforce founder rule post-check: any non-0 parent IDs must refer to an existing IndividualID
    ids = {p.IndividualID for p in people}
    for p in people:
        if p.PaternalID != "0" and p.PaternalID not in ids:
            p.PaternalID = "0"
            warnings.append(f"PaternalID for {p.IndividualID} not found; set to '0' per founder rule.")
        if p.MaternalID != "0" and p.MaternalID not in ids:
            p.MaternalID = "0"
            warnings.append(f"MaternalID for {p.IndividualID} not found; set to '0' per founder rule.")
    return Extraction(individuals=people, warnings=warnings)


def write_csv(extraction: Extraction, csv_path: str) -> None:
    import csv
    fieldnames = ["IndividualID", "Sex", "Affected", "PaternalID", "MaternalID", "Notes"]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for p in extraction.individuals:
            w.writerow(asdict(p))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Extract pedigree structure from a local file using Azure OpenAI vision.")
    ap.add_argument("path", help="Path to pedigree image or PDF.")
    ap.add_argument("--page", type=int, default=1, help="PDF page to process (1-indexed).")
    ap.add_argument("--csv", help="Optional CSV output path.")
    ap.add_argument("--endpoint", default=os.getenv("AZURE_OPENAI_ENDPOINT", "https://oric-rag-demo-resource.openai.azure.com"), help="Azure OpenAI endpoint (default: resource used in prototype.py)")
    ap.add_argument("--api-key", default=os.getenv("AZURE_OPENAI_API_KEY"), help="Azure OpenAI API key (or set AZURE_OPENAI_API_KEY)")
    ap.add_argument("--deployment", default=os.getenv("AZURE_OPENAI_DEPLOYMENT", "gpt-4.1"), help="Azure OpenAI deployment name (e.g., gpt-4.1, gpt-4o)")
    ap.add_argument("--api-version", default=os.getenv("AZURE_OPENAI_API_VERSION", "2024-02-15-preview"), help="API version to use (default matches prototype.py)")
    ap.add_argument("--page-image", help="If PDF tools are unavailable, supply a pre-rendered page image instead (PNG/JPG).")
    args = ap.parse_args(argv)

    # Show effective deployment too (helps catch 404 DeploymentNotFound)
    print(f"  deployment  = {args.deployment}")

    # Load bytes from file
    if args.page_image:
        img_bytes, mime = _read_local_bytes(args.page_image)
        if mime not in SUPPORTED_IMAGE_TYPES:
            raise SystemExit(f"Unsupported --page-image type: {mime}")
    else:
        data, mime = _read_local_bytes(args.path)
        if mime in SUPPORTED_PDF_TYPES:
            img_bytes = _load_image_from_pdf(args.path, args.page)
            mime = "image/png"
        elif mime in SUPPORTED_IMAGE_TYPES:
            img_bytes = data
        else:
            # Fall back by extension if mimetype failed
            ext = os.path.splitext(args.path)[1].lower()
            if ext in {".png", ".jpg", ".jpeg", ".webp"}:
                img_bytes = data
                mime = {
                    ".png": "image/png",
                    ".jpg": "image/jpeg",
                    ".jpeg": "image/jpeg",
                    ".webp": "image/webp",
                }[ext]
            elif ext == ".pdf":
                img_bytes = _load_image_from_pdf(args.path, args.page)
                mime = "image/png"
            else:
                raise SystemExit(f"Unsupported file type: {mime} ({ext})")

    # Build client and call model
    client = build_client(endpoint=args.endpoint, api_key=args.api_key, api_version=args.api_version)
    payload = call_model_with_image(client, args.deployment, img_bytes, mime)
    extraction = validate_payload(payload)

    # Print JSON to stdout
    out = {
        "individuals": [asdict(p) for p in extraction.individuals],
        "warnings": extraction.warnings,
    }
    print(json.dumps(out, indent=2, ensure_ascii=False))

    # Optionally write CSV
    if args.csv:
        write_csv(extraction, args.csv)
        print(f"\nCSV written to: {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
