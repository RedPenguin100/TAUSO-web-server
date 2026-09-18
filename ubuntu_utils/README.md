# Bringing TAUSO up on an Ubuntu machine

Scripts for standing this server up on a fresh Ubuntu box. The numbered ones run in order; the
other two are helpers you may or may not need.

Everything here assumes you are in the repository root or running the scripts by path — they
`cd` to the repo themselves, so `./ubuntu_utils/04-start.sh` works from anywhere in the tree.

## The short version

```bash
git clone https://github.com/RedPenguin100/TAUSO-web-server.git
cd TAUSO-web-server
./ubuntu_utils/01-install-docker.sh     # then log out and back in
./ubuntu_utils/02-fix-permissions.sh
./ubuntu_utils/03-size-for-this-machine.sh
./ubuntu_utils/04-start.sh
```

Then open <http://localhost:8501>. On an empty data directory the first boot builds everything
before the UI answers — see *How long* below.

## The scripts

| Script | What it does | When |
|---|---|---|
| `01-install-docker.sh` | Docker Engine + Compose v2 from Docker's own apt repo, and adds you to the `docker` group | Once, on a new machine |
| `02-fix-permissions.sh` | Makes docker usable without `sudo` and the data directory writable by the container | Once, and after any `sudo docker` mistake |
| `03-size-for-this-machine.sh` | Writes `docker-compose.override.yml` sized to this machine's cores and RAM | Once, and again if the hardware changes |
| `04-start.sh` | Creates `.env.local` if missing, builds, starts | Every time you deploy |
| `05-health-check.sh` | Read-only check of container, web, data, and the two files that fail late | Whenever something looks wrong |
| `06-smoke-test.sh` | Runs the repository's own eight checks inside the container | After a deploy, or when a design fails for no clear reason |
| `copy-data-from.sh` | rsync `.tauso_data` from a machine that has it, instead of rebuilding | Instead of the 4–8 hour build |
| `mount-usb.sh` | Lists block devices and mounts one — Ubuntu Server auto-mounts nothing | If the data arrives on a drive |

`01` needs `sudo`; the rest do not.

## Two routes to the data directory

`.tauso_data` is ~43 GB and the repo does not carry it. Either:

**Build it** — do nothing, `04-start.sh` triggers it on first boot. Downloads the genome, DepMap
tables and half-life data, then builds the 2.7 GB Bowtie index. 4–8 hours on a slow CPU, most of
it that index, single-threaded. You get the default cohort, so the cell-line dropdown is short.

**Copy it** — `copy-data-from.sh`, bounded by your link rather than the CPU, and you get whatever
cell lines the source machine has. `--lean` skips the ~24 GB of per-cell-line expression tables
and leaves ~19 GB; copy back only the lines you actually use.

## Configuration

`04-start.sh` copies `.env.local.example` to `.env.local` if it is absent. The defaults run on
localhost with no mail and no tunnel, which is enough to work.

Fill it in only if this machine needs more. **A Cloudflare tunnel token belongs to exactly one
machine** — two connectors sharing a token both serve that hostname and Cloudflare alternates
between them, so requests would land on whichever machine answered first, with different data.
Give a second machine its own tunnel and its own hostname.

`.env.local` and `docker-compose.override.yml` are both gitignored. Never commit either.

## How long the first boot takes

| Stage | Slow laptop |
|---|---|
| `setup-genome` | 30–60 min |
| `setup-bowtie` | 1–2 h |
| `setup-depmap` | 20–40 min |
| `build-cell-context` | 10–30 min |
| `build-cohort-transcript-expression` | 20–60 min |
| everything else | ~10 min |

Run it under `tmux` — a dropped ssh session otherwise takes the setup with it.

```bash
sudo apt install -y tmux
tmux
./ubuntu_utils/04-start.sh
# detach with Ctrl+B then D; come back with: tmux attach
```

## What a job costs

Measured on HBB, 3,932 nt, with 9 workers: ~90 s, holding under 1.5 GB. Memory barely tracks gene
size — an 18,000 nt target holds about the same as an 800 nt one — so 4 GB of RAM is comfortable.
Time does track length, at roughly 25 s fixed plus 17 ms per nucleotide.

Parallelism tops out near 3× however many cores you give it: about a third of the work is serial,
and laptop CPUs drop their clock under all-core load.

## When something is wrong

Run `05-health-check.sh` first. The failures that actually happen:

**The container restarts over and over, with nothing useful in the log.** `human_tgcn_hsapi38.csv`
is missing. It ships in `assets/` and the entrypoint copies it in, so this should be fixed — if it
recurs, check that file exists in `.tauso_data`.

**"No space left on device" while the disk is empty.** `/dev/shm` ran out, not the volume. Run
`03-size-for-this-machine.sh`, which sets `shm_size`.

**A job dies with `BrokenProcessPool`.** The kernel OOM-killed a worker. Lower `TAUSO_CORES` in
`docker-compose.override.yml` before touching anything else.

**`Permission denied` on `/home/mambauser/.tauso_data`.** The image runs as uid 57439, not as you
and not as root, so a directory only you can write is closed to it. Run `02-fix-permissions.sh`.
Do not reach for `sudo docker`: it works by running the container as root, and then everything it
writes is root-owned and your next command without `sudo` fails instead.

**Edits to `app.py` do not show up.** The app is baked into the image, not mounted, so
`docker compose restart` reuses the old one. Rebuild: `./ubuntu_utils/04-start.sh`.

**Result emails link to the wrong host.** `PUBLIC_BASE_URL` in `.env.local` is baked into those
links. Set it to whatever this machine actually answers on.
