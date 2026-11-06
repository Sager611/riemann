#!/bin/sh
# rerun.sh: the smallest certified window c = 2 (lambda = sqrt 2) from scratch, one single-threaded process, under 10 s.
# Prints one JSON line: U < Lodd, U < L2, rig [true, true], even_simple true.
set -e
cd "$(dirname "$0")"
cc -O2 -ffp-contract=off -o /tmp/ws_rscan rscan.c
rm -rf /tmp/ws_scan /tmp/ws_rerun.jsonl; mkdir -p /tmp/ws_scan
OMP_NUM_THREADS=1 WZ_RSCAN=/tmp/ws_rscan WZ_SCAN=/tmp/ws_scan uv run --no-project --with mpmath --with numpy --with python-flint python3 run_windows.py /tmp/ws_rerun.jsonl 2
