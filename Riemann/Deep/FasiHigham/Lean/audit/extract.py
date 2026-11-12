
import re, os
base = 'Riemann/Deep/FasiHigham/Lean/'
want = {'Defs.lean': ['def qcoef', 'def cs', 'def c ', 'theorem conv '],
        'PS.lean': ['def qPS', 'theorem coeff_inv', 'theorem conv_unique'],
        'Main.lean': ['def vR', 'def uR', 'def R0', 'def KR', 'theorem c_pos'],
        'MainPS.lean': ['theorem fasi_higham_m3'],
        'Root.lean': ['theorem exists_root'],
        'Quant.lean': ['def qR ', "def qR'"],
        'Quant2.lean': ['theorem c_quant']}
out = []
for f, keys in want.items():
    lines = open(base + f).read().split('\n')
    for key in keys:
        for i, l in enumerate(lines):
            if key in l and not l.strip().startswith('--'):
                # include preceding docstring
                j = i
                while j > 0 and (lines[j-1].strip().startswith('/--') or (lines[j-1].strip() and not lines[j-1].startswith(('theorem','def','noncomputable','end','namespace','open','import','@')) and '-/' in ''.join(lines[j-1:i]) and not lines[j-1].strip().startswith(('·','rw','exact','have')))):
                    j -= 1
                    if lines[j].strip().startswith('/--'): break
                k = i
                block = []
                while k < len(lines):
                    block.append(lines[k])
                    if ':= by' in lines[k] or lines[k].rstrip().endswith(':=') or ':= PowerSeries' in lines[k] or (' := ' in lines[k] and 'def' in lines[i]):
                        break
                    k += 1
                # defs: include body lines until blank
                if 'def ' in lines[i]:
                    k += 1
                    while k < len(lines) and lines[k].strip():
                        block.append(lines[k]); k += 1
                out.append('-- ' + f + '\n' + '\n'.join(lines[j:i] + block))
                break
open('.scratch/fasilean/statements.txt', 'w').write('\n\n'.join(out) + '\n')
print(len(out), sum(x.count('\n') for x in out))
