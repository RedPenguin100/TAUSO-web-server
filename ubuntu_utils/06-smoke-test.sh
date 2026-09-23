#!/usr/bin/env bash
# The repository's smoke test, run inside the running container: eight checks that the pieces
# fit together, against a short synthetic sequence. About a minute.
set -euo pipefail
cd "$(dirname "$0")/.."

docker compose exec -T tauso-web bash -lc 'cd /app && micromamba run -n base python smoke_test.py'
