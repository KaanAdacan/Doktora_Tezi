#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec watch -n 15 bash ./status_local.sh
