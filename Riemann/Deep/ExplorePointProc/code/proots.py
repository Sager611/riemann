
import sys, json
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, arb, ctx, acb
ctx.prec=1500
d=json.load(open(sys.argv[1])); D=d['D']; B=[Fr(b) for b in d['B']]; s=Fr(d['s'])
H=[fmpq_poly([1]), fmpq_poly([0,2])]
for n in range(1,2*D): H.append(fmpq_poly([0,2])*H[n]-2*n*H[n-1])
P=fmpq_poly([0])
for n in range(D+1): P+= fmpq(B[n].numerator,B[n].denominator)*H[2*n]
p=fmpq_poly([P[2*k] for k in range(P.degree()//2+1)])
from flint import fmpz_poly
# real roots via complex_roots of integer-scaled poly
den=1
for k in range(p.degree()+1): den = den*p[k].q // __import__('math').gcd(den, p[k].q)
zp=fmpz_poly([int(p[k]*den) for k in range(p.degree()+1)])
rts=zp.complex_roots()
import math
for r,mult in rts:
    if abs(float(r.imag.mid()))<1e-30 or r.imag.contains(0):
        w=float(r.real.mid())
        if w>=0:
            y=math.sqrt(w); a=y*float(s)/math.sqrt(2*math.pi)
            print('real root w=%.10f  y=%.6f  alpha=%.6f'%(w,y,a))
# min of p e^{-w/2} near those roots
