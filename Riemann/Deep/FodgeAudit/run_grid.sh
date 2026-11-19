#!/bin/bash
# usage: CASES=logs/gridA_cases.json ./run_grid.sh x out.jsonl  (runs falsify.py on every case, single process)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
uv run --no-project --with mpmath --with numpy python3 - "$1" > "$2" <<'PY'
import json, sys, time
sys.argv = [sys.argv[0]]
import importlib.util
spec = importlib.util.spec_from_file_location('f', 'falsify.py'); f = importlib.util.module_from_spec(spec); spec.loader.exec_module(f)
cases = json.load(open(sys.argv[0] if False else __import__('os').environ['CASES']))
t0 = time.time()
for cs in cases:
    try:
        r = f.case(*cs); o = r['out']
        print(json.dumps({'case': cs, 'lam': r['lam'], 'nstat': r['nstat'], **{k: [v[0], v[1]] for k, v in o.items()}}), flush=True)
    except Exception as e:
        print(json.dumps({'case': cs, 'error': repr(e)}), flush=True)
print(json.dumps({'done': len(cases), 'sec': time.time()-t0}), flush=True)
PY
