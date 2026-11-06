#!/usr/bin/env python3
"""exp48 constcheck.py (self-review): every numeric token of DERIVATION.md is looked up (as text, commas removed) in
cert.py, the code files and logs/*.txt. Tokens not found must be years, identifiers, facts quoted from the TeX sources,
or roundings of printed arb balls. Run: python3 constcheck.py"""
import re, glob, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
der = open('DERIVATION.md').read()
srcs = {}
for f in ['cert.py', 'scan.c', 'conductor.py', 'lemma_check.py', 'spot.gp'] + glob.glob('logs/*.txt'): srcs[f] = open(f).read().replace(',', '')
toks = set(re.findall(r'\d+\.\d+(?:·10\^\{?-?\d+\}?)?|\d{1,3}(?:,\d{3})+|\b\d{2,}\b', der))
missing = [t for t in sorted(toks, key=lambda x: (len(x), x)) if not any(t.replace(',', '').split('·')[0] in v for v in srcs.values())]
print(len(toks), 'numeric tokens;', len(missing), 'not found in cert.py, code or logs:', ' '.join(missing))
