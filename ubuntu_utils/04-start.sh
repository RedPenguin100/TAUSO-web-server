#!/usr/bin/env bash
# Build the image and start the server. Safe to re-run: it rebuilds and recreates.
set -euo pipefail
cd "$(dirname "$0")/.."

if [ ! -f .env.local ]; then
    echo "==> No .env.local; copying the template"
    cp .env.local.example .env.local
    echo "    Defaults run on localhost only. Edit it if this machine needs a tunnel or mail."
fi

[ -f docker-compose.override.yml ] || echo "NOTE: no override file. Run 03-size-for-this-machine.sh first."

echo "==> Building"
docker compose build tauso-web

echo "==> Starting"
docker compose up -d tauso-web

cat <<'TXT'

Started. The first boot on an empty data directory runs the whole setup chain -- genome, Bowtie
index, DepMap tables -- which is 4-8 hours on a slow machine, most of it the index building
single-threaded. Later boots skip it and take seconds.

Watch it:            docker compose logs -f tauso-web
Then open:           http://localhost:8501

Run this under tmux so a dropped ssh session does not take the setup with it.
TXT
