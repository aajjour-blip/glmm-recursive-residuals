#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(dirname "$0")"
python3 finalize_glmm_1725547.py
printf '\a'
