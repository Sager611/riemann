import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from c2405lib import *
random.seed(12345)
PR = [p for p in primes_upto(int(sys.argv[1])) if p > 2]
fails = {}
def rec(name, p, pt, lhs, rhs):
    if lhs != rhs:
        fails.setdefault(name, []).append((p, pt, lhs, rhs))
def L(p): return lambda t: leg(t, p)
pts = [(random.randint(-9, 9), random.randint(-9, 9), random.randint(-9, 9), random.randint(-9, 9)) for _ in range(4)] + [(0,0,0,0),(1,0,0,0),(0,1,1,0)]
for p in PR:
    c = consts(p); l2 = c['l2']
    m1 = (p % 4 == 1)
    if m1:
        a, b, a1, b1 = c['a'], c['b'], c['a1'], c['b1']
    else:
        hm = c['hm']; eh = (-1) ** ((hm - 1) // 2)
    for (x, y, z, w) in pts:
        pt = (x, y, z, w)
        Bp = lambda j, k: leg(j + k, p)
        Bm = lambda j, k: leg(j - k, p)
        Bpm = lambda j, k: leg(j + k, p) + leg(j - k, p)
        Bp1 = lambda j, k: leg(j + k - 1, p)
        Bpp1 = lambda j, k: leg(j + k + 1, p)
        Bmp1 = lambda j, k: leg(j - k + 1, p)
        n1 = (p - 1) // 2
        # #1
        if p > 3:
            d = detM(p, 0, n1, Bp, x, y, z, w)
            r = l2 * 2**n1 * (p*b*x + a*(w*x - (y+1)*(z+1))) if m1 else 2**n1 * ((y+1)*(z+1) - w*x)
            rec('1i', p, pt, d, r)
        d = detM(p, 0, n1, Bm, x, y, z, w)
        r = a1*(w*x - (y+1)*(z+1)) + l2*p*b1*x if m1 else w*x + (1+y)*(1-z)
        rec('1ii', p, pt, d, r)
        # #2
        v = w*x - (y+1)*(z+1)
        if p > 3:
            d = detM(p, 0, (p-3)//2, Bp, x, y, z, w)
            r = l2 * 2**((p-3)//2) * ((p*b - 2*a)*x + (a - 2*b)*v) if m1 else 2**((p-3)//2) * (v - 2*x)
            rec('2i', p, pt, d, r)
        d = detM(p, 0, (p-3)//2, Bm, x, y, z, w)
        r = -a1*x - l2*b1*v if m1 else x
        rec('2ii', p, pt, d, r)
        # #3
        if p >= 5:
            d = detM(p, 0, (p-5)//2, Bm, x, y, z, w)
            r = l2*(2*a1 - p*b1)*x + (a1 - 2*b1)*((1+y)*(1+z) - w*x) if m1 else w*x + (1+y)*(1-z)
            rec('3', p, pt, d, r)
        # #4
        if p >= 7 and p % 4 == 3:
            d = detM(p, 0, (p-7)//2, Bm, x, 0, 0, 0, use_w=False, use_z=False)
            rec('4', p, pt, d, ((p-2)//3)**2 * x)
        if p > 3:
            # #5
            d = detM(p, 1, n1, Bp, x, y, z, w)
            if m1:
                cp = (p+1)*b - 2
                r = l2 * 2**n1 * (a*(w - x) + b + (b - 1)*(y+z) - cp*(w*x - y*z))
            else:
                r = -2**n1 * (w + x + eh*(y + z + 2*y*z - 2*w*x))
            rec('5i', p, pt, d, r)
            d = detM(p, 1, n1, Bm, x, y, z, w)
            if m1:
                cp1 = (p+1)*b1 - 2
                r = a1*(w - x) + l2*(b1 + (b1 - 1)*(y+z) + cp1*(y*z - w*x))
            else:
                r = w + x - eh*(y+z)
            rec('5ii', p, pt, d, r)
        if p >= 5:
            # #6
            d = detM(p, 1, (p-3)//2, Bp, x, y, z, w)
            if m1:
                dp = (p+1)*b - 2*(a+1)
                r = l2 * 2**((p-3)//2) * (b - a*x + (a - 2*b)*w + (b-1)*(y+z) + dp*(y*z - w*x))
            else:
                r = 2**((p-3)//2) * (w + x + 2*(w*x - y*z) + eh*(y + z + 2*y*z - 2*w*x))
            rec('6i', p, pt, d, r)
            d = detM(p, 1, (p-3)//2, Bm, x, y, z, w)
            if m1:
                e1 = a1 - 2*b1
                r = e1 + l2*((2*a1 - p*b1)*x - b1*w) + (e1 + 1)*(y+z) + 2*(b1 - 1)*(w*x - y*z)
            else:
                r = 1 + (1 - eh*l2)*(2*(w*x - y*z) + y - z)
            rec('6ii', p, pt, d, r)
        if p > 3:
            # #7 (i)
            d = detM(p, 0, n1, Bpm, x, y, z, w)
            vv = (y+2)*z - w*x
            if m1:
                r = l2 * Fr(p) ** Fr(p+3, 4) if False else l2 * p**((p+3)//4) * x
            else:
                r = eh * p**((p-3)//4) * (p*x + (2 - l2)*hm*vv)
            rec('7i', p, pt, d, r)
            # #8
            if m1:
                M0 = [[leg(j+k, p) - leg(j-k, p) for k in range(1, n1+1)] for j in range(1, n1+1)]
                rec('8i-a', p, (), det_bareiss(M0), (-p)**((p-1)//4))
                d = detM(p, 1, n1, Bpm, x, y, z, w)
                h2 = (p-1)//2
                r = Fr((-p)**((p-1)//4), -p) * (h2*h2*w*x - (h2*y - 1)*(h2*z - 1))  # (-p)^((p-5)/4)
                rec('8i-b', p, pt, d, r)
            # #9
            d = detM(p, 0, (p-3)//2, Bpm, x, y, z, w)
            if m1:
                r = l2 * Fr(p**((p+3)//4), p**2) * (p*x - w*x + (y+2)*(z+2))  # p^((p-5)/4)
                rec('9a', p, pt, d, r)
        # #10
        if p > 3:
            d = detM(p, 1, n1, Bp1, x, y, z, w)
            r = l2 * 2**n1 * ((y*z - (w+1)*x)*a + (w*(1-x) + (y+1)*(z+1))*b) if m1 else 2**n1 * (y*z - (w+1)*x)
            rec('10i', p, pt, d, r)
            d = detM(p, 1, (p+1)//2, Bp1, x, y, z, w)
            r = l2 * 2**n1 * (p*b*((w+1)*x - y*z) + a*(w*(x-1) - (y+1)*(z+1))) if m1 else 2**n1 * (w*(1-x) + (y+1)*(z+1))
            rec('10ii', p, pt, d, r)
            # #11
            d = detM(p, 1, (p-3)//2, Bp1, x, y, z, w)
            q = Fr(2)**(Fr(p-5, 2)) if False else Fr(2**((p-3)//2), 2)  # 2^((p-5)/2)
            if m1:
                r = l2*q*((2 - l2)*a - p*b)*x + l2*q*(a + (l2 - 2)*b)*(w + y + z + 1) + l2*q*((p-1)*b + (l2 - 1)*(a + b))*(y*z - w*x)
            else:
                r = q*((l2 - 2)*x - w - y - z - 1 + (1 - l2)*(y*z - w*x))
            rec('11', p, pt, d, r)
            # #12
            d = detM(p, 0, n1, Bp1, x, y, z, w)
            q3 = Fr(2**n1, 2)  # 2^((p-3)/2)
            if m1:
                r = l2*q3*(p*b - 2*a)*(w + y + z + 1) + l2*q3*((2*b - a)*p*x + ((p-2)*a - p*b)*(y*z - w*x))
            else:
                r = q3*(2*w*(1-x) + 2*(y+1)*(z+1) + p*((w+1)*x - y*z))
            rec('12', p, pt, d, r)
            # #13
            if p % 4 == 3:
                d = detM(p, 0, (p-3)//2, Bp1, x, y, z, w)
                q5 = Fr(2**((p-3)//2), 2)
                r = q5*((Fr(p, 2)*(2 - l2) - 4)*x + (Fr(p, 2) - 2 - l2)*(w + y + z + 1)) + q5*(Fr(p, 2)*(l2 - 1) + 2 - l2)*(y*z - w*x)
                rec('13', p, pt, d, r)
            # #14
            d = detM(p, 0, n1, Bpp1, x, y, z, w)
            if m1:
                r = l2*2**n1*p*b*(x + Fr(p-2, 2)*(y*z - w*x)) + l2*2**n1*a*(w*(x + Fr(p-2, 2)) - (y+1)*(z+1))
            else:
                r = 2**n1*(w*(Fr(p-2, 2) - x) + (y+1)*(z+1))
            rec('14', p, pt, d, r)
            # #15 (p=3 mod 4 part fully specified)
            if not m1:
                d = detM(p, 1, n1, Bpp1, x, y, z, w)
                q3 = Fr(2**n1, 2)
                r = q3*(1 - eh)*(y + z + 2*(y*z - w*x)) + q3*((p-3)*(y*z - w*x + Fr(w, 2)) - 2*x + 1)
                rec('15b', p, pt, d, r)
            # #16
            d = detM(p, 0, n1, Bmp1, x, y, z, w)
            if m1:
                r = (p*b1 - a1)*(w*(1-x) + (y+1)*(z+1)) + p*(w*x - (y+1)*z + l2*(b1 - a1)*((1+w)*x - y*z))
            else:
                r = 1 - l2*p*x + w + y + (p*l2*eh - 1)*z + (p*l2*(1 + eh) - 1)*(y*z - w*x)
            rec('16', p, pt, d, r)
            # #17
            d = detM(p, 0, (p-3)//2, Bmp1, x, y, z, w)
            if m1:
                r = (p*b1 - a1)*((w+1)*x - y*z) + l2*(w*x - (y+1)*z) + l2*(b1 - a1)*(w*(1-x) + (y+1)*(z+1))
            else:
                r = x - l2*(w + y - z + 1) - eh*z + (1 + eh - l2)*(w*x - y*z)
            rec('17', p, pt, d, r)
            # #18
            d = detM(p, 1, n1, Bmp1, x, y, z, w)
            if m1:
                r = l2*Fr(p-1, 2)*((y+1)*z - w*x) + (p*b1 - a1)*((w+1)*x - y*z) + l2*(b1 - a1)*(w*(1-x) + (y+1)*(z+1))
            else:
                r = (w+1)*x - y*z + l2*((y+1)*(z-1) - w*(x+1)) + eh*Fr(p+1, 2)*(w*x - (y+1)*z)
            rec('18', p, pt, d, r)
print('primes tested:', PR)
for k in sorted(fails, key=lambda s: (int(''.join(ch for ch in s if ch.isdigit())), s)):
    v = fails[k]
    ps = sorted(set(t[0] for t in v))
    print('FAIL', k, 'primes', ps, 'first:', v[0])
names = ['1i','1ii','2i','2ii','3','4','5i','5ii','6i','6ii','7i','8i-a','8i-b','9a','10i','10ii','11','12','13','14','15b','16','17','18']
print('OK (no failure):', [n for n in names if n not in fails])
