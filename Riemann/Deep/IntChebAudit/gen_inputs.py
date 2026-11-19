#!/usr/bin/env python3
"""gen_inputs.py CERT.json OUTDIR X1 -- write lowcert/gcert stdin files exactly as verify.py check_lower does (imports its write_lists)."""
import sys, json, os
sys.dont_write_bytecode = True   # never write __pycache__ into Riemann/Deep/IntCheb/
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'IntCheb'))
import verify as Vf
C = json.load(open(sys.argv[1])); out = sys.argv[2]; X1 = float(sys.argv[3])
edges, masses, lib, mults, V = C["cell_edges"], C["masses"], [list(map(int, q)) for q in C["library"]], C["multipliers"], C["V"]
with open(os.path.join(out, 'g.txt'), 'w') as f:
    Vf.write_lists(f, edges, masses); f.write("%d\n" % len(lib))
    for q, c in zip(lib, mults): f.write("%.17g %d %s\n" % (c, len(q) - 1, " ".join(str(a) for a in reversed(q))))
with open(os.path.join(out, 'a.txt'), 'w') as f:
    Vf.write_lists(f, edges, masses); act = [(q, c) for q, c in zip(lib, mults) if c > 0]; f.write("%d\n" % len(act))
    for q, c in act: f.write("%.17g %d %s\n" % (c, len(q) - 1, " ".join(str(a) for a in reversed(q))))
    f.write("0 %.17g %.17g 45\n" % (X1, V))
print('wrote', out)
