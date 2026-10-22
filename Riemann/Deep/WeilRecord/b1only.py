# pyright: reportMissingImports=false
"""b1only.py IN OUT: copy parity-1 matrices of a full-parity grams_c part (integrity-checked) into a single-parity file."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import ctx
import bc, cert4
src, dst = sys.argv[1], sys.argv[2]
h = cert4.read_head(src); ctx.prec = int(h['prec'][1])
G = cert4.load_keys(src, {'A0raw_1', 'B0raw_1'})
res = {'A0raw_1': G['A0raw_1'], 'B0raw_1': G['B0raw_1'], 'N': int(h['N'][1]), 'prec': int(h['prec'][1]), 'part': [int(x) for x in h['part'][1]],
       'tag': h['tag'][1], 'npieces': int(h['npieces'][1]), 'pieces': [int(x) for x in h['pieces'][1]], 'mode': h['mode'][1], 'parities': [1]}
bc.save_mats(dst, res)
cert4.load_keys(dst, {'A0raw_1', 'B0raw_1'})          # re-read with the integrity checks
print('ok', os.path.getsize(dst))
