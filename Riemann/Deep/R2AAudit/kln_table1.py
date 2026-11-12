#!/usr/bin/env python3
# kln_table1.py -- r2aaudit: the KLN implementation used by Robin2Adic (Ingham4 KP) against KLN Table 1 (arXiv 2101.12263),
# at KLN's own H0 = 3.0610046e10 with delta as printed; then the rows at H0 = 3e12 (what Robin2Adic uses).
# uv run --no-project --with python-flint python3 Riemann/Deep/R2AAudit/kln_table1.py
import importlib.util, sys, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('i4', os.path.join(HERE, '..', 'Ingham4', 'cert.py')); I4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(I4)
worst = 0
for p in I4.build_kln():
    A = p.logA(p.s0).exp(); B = p.Bc(p.s0)
    ra = float((A/p.Atab).mid()); rb = float((B/p.Btab).mid()); worst = max(worst, abs(ra - 1), abs(rb - 1))
    print(p.name, 'admissible' if p.ok else 'NOT-ADMISSIBLE', 'delta %s' % p.de.str(4, radius=False), 'A %.4f / table %.4f' % (float(A.mid()), float(p.Atab.mid())), 'B %.4f / table %.4f' % (float(B.mid()), float(p.Btab.mid())))
print('max relative deviation from KLN Table 1 (H0 = 3.0610046e10):', worst)
P3 = I4.build_kln_H0('3e12')
print('rows admissible at H0 = 3e12:', len(P3), '; delta in [%.3f, %.3f]' % (min(float(p.de.mid()) for p in P3), max(float(p.de.mid()) for p in P3)),
      '; all Vinf > 1:', all(I4.Vinf(p.al, p.k, p.de) > 1 for p in P3) if hasattr(I4, 'Vinf') else 'n/a')
