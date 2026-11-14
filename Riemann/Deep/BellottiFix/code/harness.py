# Runs a copy of Bellotti's constants.py once per row of the paper's tables (arXiv 2405.12545v2, Tables 2-6),
# substituting t=T0, tupp=T1, sigma=alpha0, u,x,v,w (as the paper's appendix instructs), and records every printed constant.
# usage: python3 harness.py <constants_script.py> <tables.json> <out.json> [dps]
import sys, re, io, json, contextlib
src0 = open(sys.argv[1]).read()
tabs = json.load(open(sys.argv[2]))
dps = sys.argv[4] if len(sys.argv) > 4 else None
def texT(s):
    s = s.strip()
    if s.startswith('3\\cdot'): return '3*power(10,12)'
    m = re.match(r'\\exp\((.*)\)', s); return 'exp(%s)' % m.group(1)
def num(s):
    s = s.strip()
    m = re.match(r'([0-9.]+)\\cdot 10\^\{?(-?[0-9]+)\}?', s)
    return float(m.group(1)) * 10 ** int(m.group(2)) if m else float(s)
keys = ['Constant C', 'Constant B', 'd_1,1', 'd_1,2', 'd_2,1', 'd_2,2', 'd_3', 'd_5', 'd_4', 'c_1', 'c_2', 'c_3', 'c_4', 'c_5', 'b2', 'b1']
out = []
for i, (cb, ab, p, d, c) in enumerate(zip(tabs['CB'], tabs['ab'], tabs['uxvw'], tabs['d'], tabs['c'])):
    T0, T1, a0 = texT(cb[0]), texT(cb[1]), cb[2]
    u, x, v, w = p[2:6]
    s = src0
    for pat, rep in [(r'(?m)^t = .*$', 't = ' + T0), (r'(?m)^tupp= .*$', 'tupp= ' + T1), (r'(?m)^u=.*$', 'u=' + u),
                     (r'(?m)^x=.*$', 'x=' + x), (r'(?m)^v=.*$', 'v=' + v), (r'(?m)^w=.*$', 'w=' + w), (r'(?m)^sigma=.*$', 'sigma=' + a0)]:
        s, n = re.subn(pat, rep, s)
        assert n == 1, (pat, n)
    if dps: s, n = re.subn(r'(?m)^mpm\.mp\.dps = .*$', 'mpm.mp.dps = ' + dps, s); assert n == 1
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(s, 'row%d' % i, 'exec'), {'__name__': 'row'})
    got = {}
    for ln in buf.getvalue().splitlines():
        k, _, val = ln.partition(':')
        if k in keys: got[k] = float(val)
        elif k.strip() == 'extra': got.update(json.loads(val))
    paper = {'Constant C': num(cb[3]), 'Constant B': num(cb[4]), 'b1': num(ab[3]), 'b2': num(ab[4]),
             'd_1,1': num(d[2]), 'd_1,2': num(d[3]), 'd_2,1': num(d[4]), 'd_2,2': num(d[5]), 'd_3': num(d[6]), 'd_4': num(d[7]), 'd_5': num(d[8]),
             'c_1': num(c[2]), 'c_2': num(c[3]), 'c_3': num(c[4]), 'c_4': num(c[5]), 'c_5': num(c[6])}
    out.append({'row': i + 1, 'T0': cb[0], 'T1': cb[1], 'alpha0': a0, 'uxvw': [u, x, v, w], 'code': got, 'paper': paper})
    print(i + 1, cb[0], cb[1], 'C code %.6g paper %.6g rel %+.4f%%' % (got['Constant C'], paper['Constant C'], 100 * (got['Constant C'] / paper['Constant C'] - 1)), flush=True)
json.dump(out, open(sys.argv[3], 'w'), indent=1)

