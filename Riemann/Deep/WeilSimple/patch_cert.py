
import re, sys
src = open(sys.argv[1]).read()
def rep(old, new, cnt=1):
    global src
    assert src.count(old) == cnt, (old, src.count(old))
    src = src.replace(old, new)
rep("# cert_zeta.py: Chekmachev's periodized lattice minorant", "# cert_simple.py (weilsimple): cert_zeta.py of ../WeilZeta (copied unchanged except the parts marked WS) plus a deflated\n# even certificate: sector 'defl' proves Q(f) - eps ||f||^2 + M <f, g>^2 >= 0 for all even f in L^2(-L/2, L/2), with g the\n# approximate ground state of the even block of QW_lambda^NG (weil_matrix.py; env WS_NG, WS_M), and prints\n# U = Q(g)/||g||^2 (Arb). Then: even f orthogonal to g have Q(f) >= eps ||f||^2 (Courant-Fischer: lambda_2^even >= eps).\n# Original header of cert_zeta.py:\n# cert_zeta.py: Chekmachev's periodized lattice minorant")
rep("from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq\n", "from flint import arb, acb, arb_mat, acb_mat, ctx, fmpq\nsys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\nfrom weil_matrix import WeilMat   # WS\n\ndef ground_vec(T, n, its=8):   # WS: inverse iteration on mid(T); returns exact dyadic coefficients (midpoints)\n    rnd = random.Random(7); v = arb_mat(n, 1)\n    for i in range(n): v[i, 0] = arb(rnd.uniform(-1, 1))\n    for it in range(its):\n        y = T.solve(v, algorithm='approx')\n        nrm = sum((y[i, 0]*y[i, 0] for i in range(n)), arb(0)).sqrt()\n        for i in range(n): v[i, 0] = (y[i, 0]/nrm).mid()\n    return [v[i, 0] for i in range(n)]\n\ndef rayleigh(T, g):   # WS: rigorous enclosure of g^T T g / g^T g\n    n = len(g); num = arb(0)\n    for i in range(n):\n        s = arb(0)\n        for j in range(n): s += T[i, j]*g[j]\n        num += g[i]*s\n    return num, sum((x*x for x in g), arb(0))\n")
rep("cq = parse_c(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; ctx.prec = int(sys.argv[4])",
    "cq = parse_c(sys.argv[1]); sector = sys.argv[2]; epsstr = sys.argv[3]; ctx.prec = int(sys.argv[4])\n    defl = (sector == 'defl')   # WS\n    if defl: sector = 'even'")
rep("    tsym = time.time() - t0\n",
    "    tsym = time.time() - t0\n    gc = None; Uenc = None; NG = None   # WS: ground-state vector g and U = Q(g)/||g||^2\n    if defl:\n        NG = int(os.environ.get('WS_NG', '80')); Mw = arb(os.environ.get('WS_M', '1'))\n        WM = WeilMat(cq); Tg, _ = WM.block(NG, 'even')\n        gc = ground_vec(Tg, NG + 1)\n        num, gg = rayleigh(Tg, gc); Uenc = num/gg\n        del Tg\n        Gset = set(G)\n        assert all(l in dl for l in range(NG + 1))\n        pl = [(dl[l] if l in Gset else deff) for l in range(NG + 1)]\n    tU = time.time() - t0\n")
rep("    if sector == 'even': ucols.append(('v', None, arb(1), yv, bv)); Md.append(arb(2))\n",
    "    if sector == 'even': ucols.append(('v', None, arb(1), yv, bv)); Md.append(arb(2))\n    if defl: ucols.append(('g', None, arb(1), None, None)); Md.append(Mw)   # WS: + M g g^T in the positive part\n")
rep("        if col[0] == 'phi': continue\n", "        if col[0] in ('phi', 'g'): continue   # WS\n")
rep("        ta, ia = a[0], a[1]; tb, ib = b[0], b[1]\n",
    "        ta, ia = a[0], a[1]; tb, ib = b[0], b[1]\n        if ta == 'g' or tb == 'g':   # WS: g = sum_l gc[l] phi_l\n            if ta == 'g' and tb == 'g': return sum((x*x for x in gc), arb(0))\n            o = b if ta == 'g' else a\n            return sum((gc[l]*ipf(('phi', l), o) for l in range(len(gc))), arb(0))\n")
rep("        a = allc[p]; b = allc[qq]\n",
    "        a = allc[p]; b = allc[qq]\n        if a[0] == 'g' or b[0] == 'g':   # WS: <g, P0^{-1} b> = sum_l gc[l]/p_l <phi_l, b>, P0 phi_l = p_l phi_l\n            if a[0] == 'g' and b[0] == 'g': v = sum((gc[l]*gc[l]/pl[l] for l in range(len(gc))), arb(0))\n            else:\n                o = b if a[0] == 'g' else a\n                v = sum((gc[l]/pl[l]*ipf(('phi', l), o) for l in range(len(gc))), arb(0))\n            return v*a[2]*b[2]\n")
rep("    res = dict(char='zeta', c=str(cq), sector=sector,", "    res = dict(char='zeta', c=str(cq), sector=('defl' if defl else sector), NG=NG, M=os.environ.get('WS_M', '1') if defl else None,\n               U=(Uenc.str(12, radius=True) if Uenc is not None else None), U_up=(arb(Uenc.upper()).str(12, radius=False) if Uenc is not None else None),")
open(sys.argv[2], 'w').write(src)
print('patched ok')
