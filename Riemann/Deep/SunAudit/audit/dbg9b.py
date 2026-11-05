import sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from c2405lib import *
for p in [q for q in primes_upto(110) if q % 4 == 3 and q > 3]:
    c = consts(p); hm = c['hm']
    Bpm = lambda j,k: leg(j+k,p)+leg(j-k,p)
    f = lambda x,y,z,w: detM(p, 0, (p-3)//2, Bpm, x, y, 0, w, use_z=False)
    al = f(0,0,0,0); be = f(1,0,0,0)-al; ga = f(0,1,0,0)-al; ep = f(0,0,0,1)-al
    ze = f(1,0,0,1)-al-be-ep
    s = (-1)**((hm+1)//2) * Fr(p**((p+1)//4), p**2)
    mp = Fr(ze, -2*s)
    ok = al == 0 and ga == 0 and ep == 0 and be == s*p and mp.denominator == 1 and f(3,5,0,-2) == s*(p-2*mp*(-2))*3
    print(p, 'no-z version ok:', ok, 'm_p =', mp)
