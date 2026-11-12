import json, re
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq
d=json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExplorePointProc/certs/cert768.json'))
B=[Fr(b) for b in d['B']]; B[0]+=Fr(1,1000); D=d['D']
def q2(x): return fmpq(x.numerator, x.denominator)
H=[fmpq_poly([1]), fmpq_poly([0,2])]
for n in range(1,2*D): H.append(fmpq_poly([0,2])*H[n]-2*n*H[n-1])
P=fmpq_poly([0]); Q=fmpq_poly([0])
for n in range(D+1):
    P+= q2(B[n])*H[2*n]; Q+= q2(B[n]*(-1)**n)*H[2*n]
p=[P[2*k] for k in range(P.degree()//2+1)]; q=[Q[2*k] for k in range(Q.degree()//2+1)]
def parse(fn):
    out=[]
    for tok in open(fn).read().split(','):
        tok=tok.strip()
        m=re.match(r'\((-?\d+) : ℚ\)',tok)
        if m: out.append(Fr(int(m.group(1)))); continue
        a,b=tok.split('/'); out.append(Fr(int(a),int(b)))
    return out
lp=parse('/Users/adriansager/playground/riemann/.scratch/pplean/pcoef.txt'); lq=parse('/Users/adriansager/playground/riemann/.scratch/pplean/qcoef.txt')
fp=[Fr(int(c.p),int(c.q)) for c in p]; fq=[Fr(int(c.p),int(c.q)) for c in q]
print('p match', lp==fp, len(lp), len(fp)); print('q match', lq==fq, len(lq), len(fq))
A=[Fr(a) for a in d['A']]; lA=parse('/Users/adriansager/playground/riemann/.scratch/pplean/nodesA.txt'); print('A match', lA==A, len(lA))
