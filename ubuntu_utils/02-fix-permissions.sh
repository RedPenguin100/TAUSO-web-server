#!/usr/bin/env bash
# Make docker usable without sudo, and make .tauso_data writable by the container (uid 57439).
#
# Never `sudo docker ...` for this stack -- user mode only. Sudo works, then leaves root-owned
# files your user cannot touch, which fails much later in the setup chain.
#
# A shell's groups are fixed at login, so `usermod -aG docker` cannot affect the shell that ran
# it. A console login on tty1 from before the script keeps saying "permission denied" while a
# fresh ssh session on the same machine works fine. Log out and back in there -- re-running this
# script will not help.
set -euo pipefail

SELF=$(readlink -f "$0")   # before the cd: the re-exec below starts from elsewhere
cd "$(dirname "$SELF")/.."

DATA=.tauso_data
CONTAINER_UID=57439

echo "==> docker group"
if id -nG | tr ' ' '\n' | grep -qx docker; then
    echo "    $USER is in the docker group, and this shell has picked it up"
else
    if id -nG "$USER" | tr ' ' '\n' | grep -qx docker; then
        echo "    $USER is already in the docker group, but this shell predates it"
    else
        sudo usermod -aG docker "$USER"
        echo "    added $USER to the docker group"
    fi
    # `sg` gives one shell the new group, so the checks below test the socket rather than fail
    # for a reason this script just created.
    if [ -z "${TAUSO_FIXPERMS_REEXEC:-}" ] && command -v sg >/dev/null 2>&1; then
        echo "    re-running under the new group"
        export TAUSO_FIXPERMS_REEXEC=1
        exec sg docker -c "$(printf '%q ' "$SELF" "$@")"
    fi
    RELOGIN=1
fi

echo "==> $DATA"
mkdir -p "$DATA"

# 777 here so both uid 57439 and you can write without sudo. Contents keep normal modes.
chmod 777 "$DATA" 2>/dev/null || sudo chmod 777 "$DATA"
echo "    $(stat -c '%A %U:%G' "$DATA")  <- container uid $CONTAINER_UID can write here"

# Residue of an earlier `sudo docker`. Only files neither of you can read are a problem.
stuck=$(find "$DATA" -user root ! -perm -o+r -print -quit 2>/dev/null || true)
if [ -n "$stuck" ]; then
    echo "    root-owned and unreadable files found (left by sudo docker); taking ownership"
    sudo chown -R "$USER:$USER" "$DATA"
    chmod 777 "$DATA"
fi

echo "==> checking"
# Daemon reachability and volume writability fail for different reasons, so report them apart.
if ! docker info >/dev/null 2>&1; then
    echo "    cannot talk to the docker daemon. In order of likelihood:"
    echo "      - this shell still has the old groups:  log out and back in, or run 'newgrp docker'"
    echo "      - the daemon is not running:            sudo systemctl enable --now docker"
    echo "      - docker came from snap, not apt:       snap list docker"
    echo "        (the snap build confines the socket differently; 01-install-docker.sh uses apt)"
    exit 1
fi
echo "    docker answers without sudo"

if ! docker image inspect tauso-web-server-tauso-web >/dev/null 2>&1; then
    echo "    the image is not built yet, which is fine -- run 04-start.sh"
elif docker run --rm -v "$PWD/$DATA:/home/mambauser/.tauso_data" --entrypoint bash \
       tauso-web-server-tauso-web -lc 'touch /home/mambauser/.tauso_data/.writetest \
       && rm -f /home/mambauser/.tauso_data/.writetest' 2>/dev/null; then
    echo "    the container can write to $DATA"
else
    echo "    the container could NOT write to $DATA -- see the mode printed above"
    exit 1
fi

if [ -n "${RELOGIN:-}" ]; then
    echo
    echo "NOTE: the checks above ran under a temporary shell that had the docker group. Your own"
    echo "      shell still does not. Log out and back in (or run 'newgrp docker') before running"
    echo "      docker yourself, or you will get 'permission denied' despite the results above."
fi

echo
echo "Never 'sudo docker compose ...' for this stack -- user mode only, or you get root-owned"
echo "files your user cannot touch and a setup failure much later."
