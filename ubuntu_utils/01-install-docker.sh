#!/usr/bin/env bash
# Docker Engine from Docker's own apt repo. Not snap (confined; breaks the .tauso_data bind
# mount) and not docker.io (trails, and no reliable Compose v2 plugin).
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
sudo usermod -aG docker "$USER"

echo
docker --version
echo "Compose: $(docker compose version 2>/dev/null || echo 'MISSING - the plugin did not install')"
echo
echo "Log out and back in (or run: newgrp docker) before the group change takes effect."
