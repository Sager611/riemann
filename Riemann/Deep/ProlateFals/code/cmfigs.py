# Compare our t_n (even, lambda = sqrt2) with the blue dots of Connes–Moscovici's figures (arXiv 2112.05500 e-print:
# compareplots1-3.pdf).  Get the PDFs:  curl -s -L -o cm.tgz https://arxiv.org/e-print/2112.05500 && mkdir cm && tar -xzf cm.tgz -C cm
# Usage: uv run --no-project --with mpmath python3 cmfigs.py <dir containing cm/> <tn file: n t_n ...>
import re, zlib, math, sys, os
from collections import Counter
import mpmath as mp
D = sys.argv[1] if len(sys.argv) > 1 else '.'
TN = sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, 'tn_even_sqrt2.txt')
os.chdir(D)

import re, zlib
from collections import Counter
def parse(fn):
    data=open('cm/'+fn,'rb').read()
    ss=[zlib.decompress(m.group(1)) for m in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', data, re.S)]
    txt=ss[0].decode('latin1')
    toks=txt.replace('\n',' ').split()
    color=None; dots=[]; lines=[]; stack=[]; segs=[]; path=[]
    for t in toks:
        if t=='rg': color=tuple(round(float(x),3) for x in stack[-3:]); stack=[]
        elif t=='m': path=[('m',float(stack[-2]),float(stack[-1]))]; stack=[]
        elif t=='l': path.append(('l',float(stack[-2]),float(stack[-1]))); stack=[]
        elif t=='c': segs.append(tuple(float(x) for x in stack[-6:])); stack=[]
        elif t in ('f*','f'):
            if segs:
                xs=[p[4] for p in segs]; ys=[p[5] for p in segs]
                dots.append((color,(max(xs)+min(xs))/2,(max(ys)+min(ys))/2))
            segs=[]; stack=[]
        elif t=='S':
            if len(path)==2: lines.append((path[0][1],path[0][2],path[1][1],path[1][2]))
            segs=[]; stack=[]
        else: stack.append(t)
    bts=re.findall(r'BT(.*?)ET', txt, re.S)
    return dots,lines,bts

zs = [float(mp.zetazero(n).imag) for n in range(1, 101)]
tn = [float(l.split()[1]) for l in open(TN)]
for fn in ('compareplots1.pdf', 'compareplots2.pdf', 'compareplots3.pdf'):
    dots, lines, bts = parse(fn)
    red = sorted([d for d in dots if d[0] == (1.0, 0.0, 0.0)], key=lambda d: d[1])
    blue = sorted([d for d in dots if d[0] == (0.0, 0.0, 1.0)], key=lambda d: d[1])
    ys = [d[2] for d in red]; g = zs[:len(red)]
    n = len(ys); mg = sum(g)/n; my = sum(ys)/n
    b = sum((x-mg)*(y-my) for x, y in zip(g, ys))/sum((x-mg)**2 for x in g); a = my - b*mg
    rr = [y-(a+b*x) for x, y in zip(g, ys)]
    tb = [(d[2]-a)/b for d in blue]
    diff = [tb[i]-tn[i] for i in range(len(tb))]
    best = None
    for sh in (-2, -1, 0, 1, 2):
        ds = [tb[i]-tn[i+sh] for i in range(len(tb)) if 0 <= i+sh < len(tn)]
        rms = math.sqrt(sum(x*x for x in ds)/len(ds))
        if best is None or rms < best[1]: best = (sh, rms)
    print(f'{fn}: {len(red)} red, {len(blue)} blue; scale {abs(b):.4f} pt/unit, red-fit max resid {max(abs(x) for x in rr):.3f} pt ({max(abs(x) for x in rr)/abs(b):.3f} units)')
    print(f'   blue (CM plotted) - ours t_n: max |diff| = {max(abs(x) for x in diff):.3f} units, rms = {math.sqrt(sum(x*x for x in diff)/len(diff)):.3f}; rounding unit = {1/abs(b):.3f}; best index shift {best}')
    print('   first 8 CM-blue:', [round(x, 2) for x in tb[:8]]); print('   first 8 ours   :', [round(x, 2) for x in tn[:8]])
