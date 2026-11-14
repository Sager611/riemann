
import json
L = lambda n: json.load(open('logs/%s.json' % n))
o, o30, f, f30, c3, c5 = L('orig'), L('orig30'), L('fix'), L('fix30'), L('c3only'), L('c5only')
mx = lambda a, b, k='Constant C': max(abs(y['code'][k]/x['code'][k]-1) for x, y in zip(a, b))
print('dps30 vs dps10: orig max rel diff C %.2e, fix %.2e' % (mx(o, o30), mx(f, f30)))
for k in ['c_1','c_2','c_3','c_4','c_5','b1','b2']:
    print(' ', k, 'orig dps30 vs 10 max rel %.2e' % mx(o, o30, k))
print('row: c3only/orig c5only/orig fix/orig (product)')
s = []
for a, x, y, z in zip(o, c3, c5, f):
    A = a['code']['Constant C']
    s.append('%d %.3f %.3f %.3f (%.3f)' % (a['row'], x['code']['Constant C']/A, y['code']['Constant C']/A, z['code']['Constant C']/A, x['code']['Constant C']*y['code']['Constant C']/A**2))
for i in range(0, len(s), 3): print(' | '.join(s[i:i+3]))
