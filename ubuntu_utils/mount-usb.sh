#!/usr/bin/env bash
# Find and mount a USB drive. Ubuntu Server has no desktop session, so nothing auto-mounts:
# the drive appears as a block device only.
set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Block devices:"
    lsblk -o NAME,SIZE,FSTYPE,LABEL,MOUNTPOINT
    echo
    echo "Usage: $0 /dev/sdb1 [mountpoint]     (default mountpoint: /mnt/usb)"
    echo "Unmount with: sudo umount /mnt/usb   -- always, before pulling the drive out"
    exit 0
fi

DEV=$1
MNT=${2:-/mnt/usb}
sudo mkdir -p "$MNT"
sudo mount "$DEV" "$MNT" || {
    echo "Mount failed. A missing driver is the usual cause:"
    echo "  exFAT: sudo apt install exfatprogs"
    echo "  NTFS:  sudo apt install ntfs-3g"
    exit 1
}
echo "Mounted $DEV at $MNT"
df -h "$MNT" | tail -1
