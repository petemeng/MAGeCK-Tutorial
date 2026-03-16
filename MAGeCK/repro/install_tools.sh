#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
python3 -m pip install --user "$ROOT/tmp/liulab-mageck-e94a4f728a7b"
python3 -m pip install --user "$ROOT/tmp/liulab-mageck-vispr-72d1b1dd9277"

echo 'Installed mageck and mageck-vispr into ~/.local/bin'
