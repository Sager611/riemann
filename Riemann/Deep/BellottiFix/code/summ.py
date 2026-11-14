
import json
L = lambda n: json.load(open('logs/%s.json' % n))
o, f, f16, b, c5 = L('orig'), L('fix'), L('fix16'), L('both'), L('c5only')
print('row | Cpaper | Ccode | Cfix ratio | fix16 ratio | c5only | c5: code -> fix | cond4 fix')
for a, x, y, z, w in zip(o, f, f16, b, c5):
    A = a['code']['Constant C']
    print(a['row'], '| %.6g | %.6g | %.6g %.4f | %.4f | %.4f | %.4f -> %.4f | %.4f' % (a['paper']['Constant C'], A, x['code']['Constant C'], x['code']['Constant C']/A,
          y['code']['Constant C']/A, w['code']['Constant C']/A, a['code']['c_5'], x['code']['c_5'], x['code']['cond4']))
print('min cond4 fix16:', min(r['code']['cond4'] for r in f16))
