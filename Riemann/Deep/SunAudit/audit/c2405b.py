import sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from c2405lib import *
PR = [p for p in primes_upto(int(sys.argv[1])) if p > 3]
def coeffs(f):
    al = f(0,0,0,0); be = f(1,0,0,0)-al; ga = f(0,1,0,0)-al; de = f(0,0,1,0)-al; ep = f(0,0,0,1)-al
    ze = f(1,0,0,1)-al-be-ep   # xw coefficient
    yz = f(0,1,1,0)-al-ga-de   # yz coefficient
    # consistency check: general point
    x,y,z,w = 3,-2,5,7
    assert f(x,y,z,w) == al+be*x+ga*y+de*z+ep*w+ze*x*w+yz*y*z, 'not of expected form'
    return al,be,ga,de,ep,ze,yz
out = []
for p in PR:
    c = consts(p); l2 = c['l2']; n1 = (p-1)//2
    Bpm = lambda j,k: leg(j+k,p)+leg(j-k,p)
    if p % 4 == 1:
        a, b = c['a'], c['b']
        # 7(ii)
        f = lambda x,y,z,w: detM(p, 0, n1, lambda j,k: leg(j+k,p)-leg(j-k,p), 0, y, z, 0, use_w=False)
        al,be,ga,de,ep,ze,yz = coeffs(f)
        ok = (al,be,ga,de,ep,ze)==(0,0,0,0,0,0) and yz % (4*p**((p-5)//4)) == 0
        out.append(('7ii', p, ok, 'x_p=%s' % Fr(yz, 4*p**((p-5)//4))))
        # 15 (p=1 mod 4)
        f = lambda x,y,z,w: detM(p, 1, n1, lambda j,k: leg(j+k+1,p), x, y, z, 0, use_w=False)
        al,be,ga,de,ep,ze,yz = coeffs(f)
        q3 = l2*2**((p-3)//2)
        ok = (al == q3*((2*b-a-1)+1) and be == q3*(p*b-2*a) and ga == q3*(2*b-a-1) and de == q3*(2*b-a-1) and ep == 0 and ze == 0)
        npv = Fr(yz, 2*q3) - b + a
        ok = ok and npv.denominator == 1 and npv > 0
        out.append(('15a', p, ok, 'n_p=%s' % npv))
    else:
        hm = c['hm']
        # 8(ii)
        f = lambda x,y,z,w: (-1)**((hm+1)//2) * detM(p, 1, n1, Bpm, x, y, 0, w, use_z=False)
        al,be,ga,de,ep,ze,yz = coeffs(f)
        P = p**((p-3)//4)
        ok = al == -P and ga == P*Fr(p-1,2) and be == P*(2-l2)*hm and ep == P*(2-l2)*hm and de == 0 and yz == 0
        qp = Fr(-ze*p, 16*l2*P)
        ok = ok and qp.denominator == 1
        out.append(('8ii', p, ok, 'q_p=%s' % qp))
        # 9 (p=3 mod 4)
        f = lambda x,y,z,w: detM(p, 0, (p-3)//2, Bpm, x, y, z, w)
        al,be,ga,de,ep,ze,yz = coeffs(f)
        s = (-1)**((hm+1)//2) * Fr(p**((p+1)//4), p**2)   # p^((p-7)/4)
        mp = Fr(ze, -2*s)
        ok = al == 0 and ga == 0 and de == 0 and ep == 0 and yz == 0 and be == s*p and mp.denominator == 1
        out.append(('9b', p, ok, 'm_p=%s' % mp))
bad = [o for o in out if not o[2]]
print('existential-form checks, primes 5..%s:' % PR[-1], len(out), 'checks; failures:', bad)
print([o for o in out][:12])
# 19: QR property
res19 = []
for p in [q for q in primes_upto(int(sys.argv[1])) if q > 2]:
    n1 = (p-1)//2
    if p % 4 == 1:
        for d1 in (1,-1):
            for d2 in (1,-1):
                M = [[leg(j+k,p)+leg(j-k,p)+d1*leg(j*j+d2*k*k,p) for k in range(n1+1)] for j in range(n1+1)]
                D = 2*det_bareiss(M)
                res19.append((p, d1, d2, leg(D, p)))
    else:
        M = [[leg(j+k,p)+leg(j*j+k*k,p) for k in range(n1+1)] for j in range(n1+1)]
        D = 2*det_bareiss(M)
        res19.append((p, 0, 0, leg(D, p)))
print('19 non-QR cases (legendre != 1):', [r for r in res19 if r[3] != 1])
