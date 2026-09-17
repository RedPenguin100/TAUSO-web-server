#!/usr/bin/env bash
# Claude Code, native installer. Optional: nothing in this stack needs it.
#
# The native install auto-updates in the background, which suits a machine you rarely log into.
# The apt repository exists too but needs `apt upgrade claude-code` by hand.
set -euo pipefail

curl -fsSL https://claude.ai/install.sh | bash

echo
echo "Installed to ~/.local/bin/claude. Start a new shell, then:"
echo "  claude --version"
echo
echo "Logging in on a headless box: run 'claude', open the URL it prints on a machine that has a"
echo "browser, and paste the code back. Simpler still, ssh in from a desktop and the prompt opens"
echo "there. Claude Code needs a Pro, Max, Team or Enterprise account."
