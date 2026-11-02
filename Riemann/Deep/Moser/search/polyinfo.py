
"""Exact vertices/area of the certified quadrilateral; writes polygon.json.
Usage: polyinfo.py CERT_DIR OUT_JSON"""
import json, sys
from fractions import Fraction as F
def n(t): return ((1 - t*t) / (1 + t*t), 2*t / (1 + t*t))
def meet(a, ha, b, hb):
    d = a[0]*b[1] - a[1]*b[0]
    return ((ha*b[1] - hb*a[1]) / d, (a[0]*hb - b[0]*ha) / d)
T = [F(x) for x in ["-631/58", "-1", "7/16", "37/56"]]
H0 = [F(0), F(0), F(130783, 250000), F(446199, 1000000)]
s = F(499, 500)
H = [s*h for h in H0]
N = [n(t) for t in T]
V = [meet(N[i-1], H[i-1], N[i], H[i]) for i in range(4)]
A = sum(V[i][0]*V[(i+1)%4][1] - V[i][1]*V[(i+1)%4][0] for i in range(4)) / 2
A0 = F(7170601298323360123534337, 29109471743520000000000000)
cert = json.load(open(sys.argv[1] + '/header.json'))
assert cert['h'] == [str(h) for h in H] and cert['t'] == [str(t) for t in T], (cert['h'], [str(h) for h in H])
assert A == A0 * s * s, (A, A0*s*s)
assert A == F(1785486893883814994120173447337, 7277367935880000000000000000000)
# every vertex satisfies all constraints
for v in V:
    assert all(N[j][0]*v[0] + N[j][1]*v[1] <= H[j] for j in range(4))
out = {"description": "K = {x : n(t_i).x <= h_i, i=0..3}, n(t) = ((1-t^2)/(1+t^2), 2t/(1+t^2)); vertices V_i = line(i-1) cap line(i)",
       "t": [str(t) for t in T], "h": [str(h) for h in H], "area": str(A), "area_decimal": "%.15f" % float(A),
       "vertices": [[str(c) for c in v] for v in V], "vertices_decimal": [["%.12f" % float(c) for c in v] for v in V],
       "scale_of_keller_quadrilateral": str(s), "keller_area": str(A0)}
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print('area', A, float(A), 'ratio', float(A/A0), 'gain', float(A0 - A))
print('diag check |V0-V2|^2 =', float((V[0][0]-V[2][0])**2 + (V[0][1]-V[2][1])**2), '|V1-V3|^2 =', float(V[3][0]**2 + V[3][1]**2))
