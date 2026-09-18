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
cd "$(dirname "$0")/.."

DATA=.tauso_data
CONTAINER_UID=57439

echo "==> docker group"
if id -nG | tr ' ' '\n' | grep -qx docker; then
    echo "    $USER is in the docker group, and this shell has picked it up"
else
    sudo usermod -aG docker "$USER"
    echo "    added $USER to the docker group"
    echo "    this shell still has the old groups -- run 'newgrp docker' or log out and back in"
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
if docker run --rm -v "$PWD/$DATA:/home/mambauser/.tauso_data" --entrypoint bash \
     tauso-web-server-tauso-web -lc 'touch /home/mambauser/.tauso_data/.writetest \
     && rm -f /home/mambauser/.tauso_data/.writetest' 2>/dev/null; then
    echo "    the container can write to $DATA"
else
    echo "    could not verify -- the image may not be built yet, which is fine; run 04-start.sh"
fi

echo
echo "Never run 'sudo docker compose ...' for this stack. It works, and it leaves behind files"
echo "your own user cannot touch, which surfaces much later as a failure in the setup chain."
