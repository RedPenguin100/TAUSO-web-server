#!/usr/bin/env bash
# The repository's own smoke test, run inside the running container.
#
# Eight checks that the pieces fit together: the tauso API the pipeline imports, bowtie on PATH,
# the model present and its features computable, a design running end to end, blank experimental
# conditions reaching the model as missing rather than zero, sugar and backbone validation, the
# results page picking the score column rather than a feature, and a job surviving a write and
# read back. It designs against a short synthetic sequence, so it finishes in about a minute.
set -euo pipefail
cd "$(dirname "$0")/.."

docker compose exec -T tauso-web bash -lc 'cd /app && micromamba run -n base python smoke_test.py'
