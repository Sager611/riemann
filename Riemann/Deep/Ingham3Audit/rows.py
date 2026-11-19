
# rows.py -- pbpaudit: which KLN row cert.py's Setup uses in each sigma-cell, admissibility margins, B_P > 0 and q1 > 0 (Arb, 256 bits)
import types
src = open('Riemann/Deep/Ingham3/cert.py').read().replace('ctx.prec = 128', 'ctx.prec = 256')
C = types.ModuleType('certP'); C.__dict__['__name__'] = 'certP'; C.__file__ = 'cert_copy.py'
exec(compile(src, 'cert.py@256', 'exec'), C.__dict__)
arb = C.arb
for de in (None, '1'):
    P = C.build_kln(de)
    S = C.Setup(P, {}, [], arb('0.515'), arb('0.0025'))
    groups = []
    for (s, t, p) in S.cells:
        if groups and groups[-1][2] is p: groups[-1][1] = t
        else: groups.append([s, t, p])
    print('delta =', de or 'Table 1', ':', len(S.cells), 'sigma-cells; last cell ends at B:', S.cells[-1][1] is C.B)
    for s, t, p in groups:
        marg = min((c[0] - p.smin) for c in S.cells if c[2] is p)
        print('  row %s on [%s, %s]: s_min = %s, min (s - s_min) = %s, d = %s, c1 > 0: %s, q1 > 0: %s, B_P(B) = %s' % (
            p.name, s.str(6, radius=False), t.str(8, radius=False), p.smin.str(10, radius=False), marg.str(4, radius=False), p.d.str(4, radius=False),
            bool(p.c1 > 0), bool(p.q1 > 0), p.Bc(C.B).str(5, radius=False)))
    assert all(bool(p.q1 > 0) and bool(p.Bc(C.B) > 0) for p in P)
print('all rows: q1 > 0 and B_P(B) > 0 (so B_P > 0 and decreasing on [0.515, B])')
