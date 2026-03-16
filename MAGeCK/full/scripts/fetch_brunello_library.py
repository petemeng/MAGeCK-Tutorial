from __future__ import annotations

import argparse
from pathlib import Path
import sys

import requests

ADDGENE_URL = "https://media.addgene.org/cms/filer_public/f8/f7/f8f73d69-94f5-4382-bf5f-658c78e38f2f/broadgpp-brunello-library-contents.txt"


def validate_table(path: Path) -> None:
    with path.open() as handle:
        header = handle.readline().rstrip("\n")
        if not header:
            raise RuntimeError("Downloaded file is empty")
        lowered = header.lower()
        if not any(token in lowered for token in ["gene", "sgrna", "sequence"]):
            raise RuntimeError(f"Unexpected header: {header}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch full Brunello library annotation from Addgene")
    parser.add_argument("--url", default=ADDGENE_URL, help="Override source URL")
    parser.add_argument("--output", default="MAGeCK/full/external/brunello/library.tsv", help="Destination path")
    args = parser.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(args.url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with out.open("wb") as handle:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    validate_table(out)
    print(f"Saved Brunello library annotation to {out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
