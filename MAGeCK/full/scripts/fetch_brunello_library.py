from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

OFFICIAL_ADDGENE_URL = (
    "https://media.addgene.org/cms/filer_public/"
    "f8/f7/f8f73d69-94f5-4382-bf5f-658c78e38f2f/"
    "broadgpp-brunello-library-contents.txt"
)
GITHUB_MIRROR_URL = (
    "https://raw.githubusercontent.com/dbrookeUAB/GeCKO/master/"
    "inst/extdata/temp/broadgpp-brunello-library-contents.txt"
)


def run_curl(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(args, check=False, text=False)


def probe_size(url: str, retries: int = 8) -> int:
    for attempt in range(1, retries + 1):
        header_path = Path('/tmp/brunello_probe.headers')
        if header_path.exists():
            header_path.unlink()
        cmd = [
            'curl', '--fail', '--location', '--ipv4',
            '--connect-timeout', '10', '--max-time', '30',
            '--silent', '--show-error',
            '--range', '0-0',
            '-D', str(header_path),
            '-o', '/tmp/brunello_probe.bin',
            url,
        ]
        proc = run_curl(cmd)
        if proc.returncode == 0 and header_path.exists():
            text = header_path.read_text(errors='replace')
            match = re.search(r'Content-Range:\s*bytes\s+0-0/(\d+)', text, flags=re.I)
            if match:
                return int(match.group(1))
        time.sleep(1)
    raise RuntimeError(f'Unable to determine remote file size for {url}')


def download_in_chunks(url: str, output_path: Path, chunk_size: int = 20000, retries: int = 10) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    chunks_dir = output_path.parent / 'chunks'
    chunks_dir.mkdir(parents=True, exist_ok=True)

    total_size = probe_size(url)
    starts = range(0, total_size, chunk_size)
    for index, start in enumerate(starts, 1):
        end = min(start + chunk_size - 1, total_size - 1)
        expected = end - start + 1
        part = chunks_dir / f'{start:08d}_{end:08d}.part'
        if part.exists() and part.stat().st_size == expected:
            continue
        if part.exists():
            part.unlink()
        tmp = part.with_suffix('.tmp')
        for attempt in range(1, retries + 1):
            if tmp.exists():
                tmp.unlink()
            cmd = [
                'curl', '--fail', '--location', '--ipv4',
                '--connect-timeout', '10', '--max-time', '30',
                '--silent', '--show-error',
                '--range', f'{start}-{end}',
                '-o', str(tmp),
                url,
            ]
            proc = run_curl(cmd)
            size = tmp.stat().st_size if tmp.exists() else 0
            if proc.returncode == 0 and size == expected:
                tmp.rename(part)
                print(f'[done] {index}: {start}-{end} ({size} bytes)')
                break
            if tmp.exists():
                tmp.unlink()
            if attempt == retries:
                raise RuntimeError(
                    f'Failed to download chunk {start}-{end} from {url} after {retries} attempts'
                )
            time.sleep(1)

    with output_path.open('wb') as handle:
        for start in starts:
            end = min(start + chunk_size - 1, total_size - 1)
            part = chunks_dir / f'{start:08d}_{end:08d}.part'
            handle.write(part.read_bytes())

    if output_path.stat().st_size != total_size:
        raise RuntimeError(
            f'Unexpected assembled size: {output_path.stat().st_size} != {total_size}'
        )


def build_mageck_library(raw_path: Path, output_path: Path, metadata_path: Path) -> tuple[int, int]:
    counts: Counter[str] = Counter()
    rows = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open() as raw_handle, output_path.open('w') as out_handle, metadata_path.open('w') as meta_handle:
        reader = csv.DictReader(raw_handle, delimiter='\t')
        fieldnames = ['sgRNA', 'Gene', 'Sequence'] + list(reader.fieldnames or [])
        meta_writer = csv.DictWriter(meta_handle, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        meta_writer.writeheader()
        for row in reader:
            gene = row['Target Gene Symbol'].strip()
            sequence = row['sgRNA Target Sequence'].strip().upper()
            if not gene or not sequence or any(base not in 'ACGT' for base in sequence):
                continue
            counts[gene] += 1
            sgrna = f'{gene}_sg{counts[gene]}'
            out_handle.write(f'{sgrna}\t{sequence}\t{gene}\n')
            meta_writer.writerow({'sgRNA': sgrna, 'Gene': gene, 'Sequence': sequence, **row})
            rows += 1
    return rows, len(counts)


def validate_output(path: Path) -> None:
    with path.open() as handle:
        first = handle.readline().rstrip('\n')
        if not first:
            raise RuntimeError('Generated library is empty')
        parts = first.split('\t')
        if len(parts) != 3:
            raise RuntimeError(f'Unexpected MAGeCK library format: {first}')


def main() -> None:
    parser = argparse.ArgumentParser(description='Fetch full Brunello library annotation and convert it to MAGeCK format')
    parser.add_argument('--url', default='', help='Override download URL')
    parser.add_argument('--output', default='MAGeCK/full/external/brunello/library.tsv', help='Destination MAGeCK library path')
    parser.add_argument('--raw-output', default='MAGeCK/full/external/brunello/library_raw.txt', help='Downloaded raw annotation path')
    parser.add_argument('--metadata-output', default='MAGeCK/full/external/brunello/library_metadata.tsv', help='Expanded metadata output path')
    parser.add_argument('--chunk-size', type=int, default=20000, help='Chunk size in bytes for resilient ranged downloads')
    parser.add_argument('--retries', type=int, default=10, help='Retries per chunk')
    args = parser.parse_args()

    output_path = Path(args.output)
    raw_path = Path(args.raw_output)
    metadata_path = Path(args.metadata_output)

    sources = [args.url] if args.url else [OFFICIAL_ADDGENE_URL, GITHUB_MIRROR_URL]
    last_error: Exception | None = None
    for source in sources:
        try:
            print(f'Fetching Brunello annotation from {source}')
            download_in_chunks(source, raw_path, chunk_size=args.chunk_size, retries=args.retries)
            rows, genes = build_mageck_library(raw_path, output_path, metadata_path)
            validate_output(output_path)
            print(f'Saved MAGeCK library to {output_path} ({rows} guides across {genes} genes/groups)')
            return
        except Exception as exc:
            last_error = exc
            print(f'WARNING: failed source {source}: {exc}', file=sys.stderr)
            if raw_path.exists():
                raw_path.unlink()
    raise SystemExit(f'ERROR: could not fetch Brunello annotation: {last_error}')


if __name__ == '__main__':
    main()
