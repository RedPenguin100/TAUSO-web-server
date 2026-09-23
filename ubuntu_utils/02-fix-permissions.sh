#!/usr/bin/env bash
# Make docker usable without sudo, and make the data directory writable by the container.
#
# These are one problem wearing two hats. The image runs as mambauser, uid 57439 -- not root and
# not you. A bind-mounted directory that only uid 1000 can write gives "Permission denied" on
# /home/mambauser/.tauso_data the moment setup tries to create anything. Reaching for `sudo docker`
# makes that go away by running the container as root, and then every file it writes is root-owned,
# so the next command you run without sudo fails instead. The fix is to stop using sudo and open
# the directory to the container's uid.
set -euo pipefail

# Resolved before the cd, because the re-exec below runs this path again from a different cwd.
SELF=$(readlink -f "$0")
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
    # A shell's groups are fixed at login, so this one cannot reach the docker socket whatever we
    # just did -- and the check at the bottom would then fail for a reason this script created,
    # reporting it as a missing image. `sg` starts one shell that does have the group, so re-exec
    # there and the checks mean something on the first run.
    if [ -z "${TAUSO_FIXPERMS_REEXEC:-}" ] && command -v sg >/dev/null 2>&1; then
        echo "    re-running under the new group"
        export TAUSO_FIXPERMS_REEXEC=1
        exec sg docker -c "$(printf '%q ' "$SELF" "$@")"
    fi
    RELOGIN=1
fi

echo "==> $DATA"
mkdir -p "$DATA"

# The container writes as 57439, you read and rsync as yourself, and neither of you should need
# sudo for it. 777 on this one directory is what makes that true; its contents keep normal modes.
chmod 777 "$DATA" 2>/dev/null || sudo chmod 777 "$DATA"
echo "    $(stat -c '%A %U:%G' "$DATA")  <- container uid $CONTAINER_UID can write here"

# Root-owned files are the residue of an earlier `sudo docker`. Only the ones neither you nor the
# container can get at are a problem: a root-owned file left world-readable is read perfectly well
# by both, and taking ownership of it would be noise.
stuck=$(find "$DATA" -user root ! -perm -o+r -print -quit 2>/dev/null || true)
if [ -n "$stuck" ]; then
    echo "    root-owned and unreadable files found (left by sudo docker); taking ownership"
    sudo chown -R "$USER:$USER" "$DATA"
    chmod 777 "$DATA"
fi

echo "==> checking"
# Reaching the daemon and writing to the volume are separate failures with separate fixes, so they
# are reported separately. The old script rolled both into "the image may not be built yet", which
# named the one cause that was usually not it.
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
echo "Never run 'sudo docker compose ...' for this stack. It works, and it leaves behind files"
echo "your own user cannot touch, which surfaces much later as a failure in the setup chain."
