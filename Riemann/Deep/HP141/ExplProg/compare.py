# compare.py -- ExplProg certificate output: original (B = 0.8733819666132076) vs this folder (B = exact theta141); prints changed lines only
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
old = [l.rstrip() for l in open(os.path.join(HERE, '..', '..', 'ExplProg', 'logs', 'cert_output.txt'))]
new = [l.rstrip() for l in open(os.path.join(HERE, 'logs', 'cert_output.txt')) if l.startswith('[')]
old = [l for l in old if l.startswith('[')]
def key(l): return re.sub(r'\[[-0-9.e+ ]+\+/- [0-9.e-]+\]|[-0-9.]+e[-+]?[0-9]+|\d+\.\d+', '#', l)
ko = {key(l): l for l in old}
ch = 0
for l in new:
    k = key(l)
    if k in ko and ko[k] != l:
        ch += 1; print('OLD:', ko[k][:200]); print('NEW:', l[:200])
    elif k not in ko: print('NEW-ONLY:', l[:200])
print(f'{ch} changed of {len(new)} lines; {len(old)} old lines')
