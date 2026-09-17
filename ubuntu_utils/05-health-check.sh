#!/usr/bin/env bash
# Is this install actually working? Read-only; changes nothing.
set -uo pipefail
cd "$(dirname "$0")/.."

DATA="${TAUSO_DATA_DIR:-./.tauso_data}"
fail=0
ok()   { printf '  \033[32mok\033[0m    %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=1; }
warn() { printf '  \033[33mwarn\033[0m  %s\n' "$1"; }

echo "Container"
if docker compose ps --status running 2>/dev/null | grep -q tauso-web; then
    ok "tauso-web is running"
else
    bad "tauso-web is not running  (docker compose logs tauso-web)"
fi

echo "Web"
code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 10 http://localhost:8501/ 2>/dev/null)
[ "$code" = "200" ] && ok "http://localhost:8501 -> 200" || bad "localhost:8501 -> ${code:-no answer}"

echo "Data directory"
[ -d "$DATA" ] && ok "$DATA exists ($(du -sh "$DATA" 2>/dev/null | cut -f1))" || bad "$DATA missing"
[ -f "$DATA/.tauso_initialized_v2" ] && ok "setup chain completed" \
    || warn "setup chain has not finished yet (first boot takes hours)"

# The two files whose absence fails late, after the long setup rather than before it.
[ -f "$DATA/human_tgcn_hsapi38.csv" ] && ok "tRNA gene counts present" \
    || bad "human_tgcn_hsapi38.csv missing -- the container will restart-loop forever"
[ -f "$DATA/general_expression.parquet" ] && ok "general expression table present" \
    || warn "general_expression.parquet missing -- the general off-target step will raise;
          fix with: docker compose exec tauso-web micromamba run -n base tauso build-general-expression"

echo "Cell lines offered"
n=$(ls "$DATA/processed_expression" 2>/dev/null | wc -l)
[ "$n" -gt 0 ] && ok "$n cell lines built" || bad "no per-cell-line expression built"

echo "Resources"
if [ -f docker-compose.override.yml ]; then
    ok "sized for this machine ($(grep -oE 'mem_limit: [0-9]+m' docker-compose.override.yml || echo '?'))"
else
    warn "no docker-compose.override.yml -- using the production box's 7g / 12 cores"
fi

echo
[ "$fail" -eq 0 ] && echo "All good." || { echo "Something above needs fixing."; exit 1; }
