from __future__ import annotations

import argparse
from pathlib import Path


def sanitize_token(token: str) -> str:
    return '_'.join(token.split())


def main() -> None:
    parser = argparse.ArgumentParser(description='Sanitize MAGeCK count table for mageck mle by replacing whitespace in sgRNA/gene identifiers.')
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    replaced = 0
    rows = 0
    with in_path.open() as src, out_path.open('w') as dst:
        for idx, line in enumerate(src, start=1):
            fields = line.rstrip('\n').split('\t')
            if idx == 1:
                dst.write(line)
                continue
            if len(fields) < 3:
                raise SystemExit(f'Unexpected field count on line {idx}: {len(fields)}')
            original = fields[:2]
            fields[0] = sanitize_token(fields[0])
            fields[1] = sanitize_token(fields[1])
            replaced += sum(o != n for o, n in zip(original, fields[:2]))
            dst.write('\t'.join(fields) + '\n')
            rows += 1
    print(f'rows={rows}')
    print(f'replaced_tokens={replaced}')


if __name__ == '__main__':
    main()
