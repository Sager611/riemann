# lowcert_dbg_patch.py: run in .scratch/intchebaudit with copy/lowcert.c present; writes bin/lowcert_dbg.c (prints leaves with lb < V + 3e-9).
src = open('copy/lowcert.c').read()
old = '    if (b >= Vt) { if (b < minleaf) minleaf = b; return; }'
assert src.count(old) == 1
new = '    if (b >= Vt) { if (b < minleaf) minleaf = b; if (b < Vt + 3e-9 && ntight < 400) { ntight++; fprintf(stderr, "TIGHT %.17g %.17g %.15g %d\\n", lo, hi, b, depth); } return; }'
src = src.replace(old, new).replace('static double Vt;', 'static double Vt; static int ntight = 0;')
open('bin/lowcert_dbg.c', 'w').write(src)
