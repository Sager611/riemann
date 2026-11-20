#!/bin/bash
# npaframe queue 11: after queue 10, tiny analyses of saved laws (form factors up to 2N, types)
cd "$(dirname "$0")"
while pgrep -f run_q10.sh > /dev/null; do sleep 5; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
uv run --no-project --with numpy python3 formfactor.py logs/N8_Po.npz logs/N8_To6.npz logs/N8_To4.npz logs/N12_Po.npz logs/N12_P.npz logs/N12_T6.npz > logs/formfactor.log 2>&1
uv run --no-project --with numpy python3 analyze_law.py logs/N8_Po.npz logs/N8_To6.npz logs/N12_P.npz logs/N12_T6.npz > logs/analyze_laws.log 2>&1
