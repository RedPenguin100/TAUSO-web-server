#!/usr/bin/env bash
# Copy .tauso_data from a machine that already has it, instead of rebuilding it here.
#
# Building from scratch is 4-8 hours on a slow machine; copying is bounded by the link. The full
# directory is ~43 GB, of which ~24 GB is per-cell-line expression -- one 1.8 MB gene table and one
# 13 MB transcript table per line. --lean skips those, leaving ~19 GB, and the dropdown then offers
# only the lines whose files are present.
set -euo pipefail
cd "$(dirname "$0")/.."

LEAN=0
[ "${1:-}" = "--lean" ] && { LEAN=1; shift; }

if [ $# -eq 0 ]; then
    cat <<'TXT'
Usage: copy-data-from.sh [--lean] SOURCE

  SOURCE   rsync source: user@host:/path/to/.tauso_data/  or  /mnt/usb/tauso_data/
           Keep the trailing slash: without it rsync nests a directory inside the target.

  --lean   skip the per-cell-line expression tables (~24 GB); copy the ones you need afterwards

Examples:
  ./ubuntu_utils/copy-data-from.sh michael@desktop:~/career/TAUSO-web-server/.tauso_data/
  ./ubuntu_utils/copy-data-from.sh --lean /mnt/usb/tauso_data/
TXT
    exit 0
fi

SRC=$1
mkdir -p .tauso_data

# --partial matters: a dropped link resumes rather than restarting 43 GB.
ARGS=(-avh --progress --partial)
if [ "$LEAN" -eq 1 ]; then
    ARGS+=(--exclude 'processed_expression/*' --exclude 'processed_transcript_expression/*')
    echo "==> Lean copy: skipping the per-cell-line expression tables"
fi

rsync "${ARGS[@]}" "$SRC" ./.tauso_data/

echo
echo "Copied. Check what arrived:  ./ubuntu_utils/05-health-check.sh"
if [ "$LEAN" -eq 1 ]; then
    echo "Then copy the cell lines you want, e.g. HepG2 (ACH-000739):"
    echo "  rsync -avh SOURCE/processed_expression/ACH-000739_expression.csv ./.tauso_data/processed_expression/"
    echo "  rsync -avh SOURCE/processed_transcript_expression/ACH-000739_transcript_expression.csv ./.tauso_data/processed_transcript_expression/"
fi
