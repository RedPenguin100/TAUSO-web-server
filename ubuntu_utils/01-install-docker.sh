#!/usr/bin/env bash
# Docker Engine from Docker's own apt repository.
#
# Not snap: the snap runs confined and bind mounts outside $HOME fail or mount the wrong thing,
# and this stack bind-mounts .tauso_data. Not `apt install docker.io` either: it trails by a
# release or two and does not reliably bring the Compose v2 plugin, which every command here uses.
set -euo pipefail

echo "==> Removing any distro Docker packages"
sudo apt-get remove -y docker docker-engine docker.io containerd runc 2>/dev/null || true

echo "==> Adding Docker's repository"
sudo apt-get update
sudo apt-get install -y ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" \
  | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

echo "==> Installing"
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

echo "==> Letting $USER run docker without sudo"
# Without this the bind-mounted data directory ends up root-owned and the container cannot write.
sudo usermod -aG docker "$USER"

echo
docker --version
echo "Compose: $(docker compose version 2>/dev/null || echo 'MISSING - the plugin did not install')"
echo
echo "Log out and back in (or run: newgrp docker) before the group change takes effect."
