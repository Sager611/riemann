# A: grhsub MV-claim threshold A0 = sup_{sigma in [1/2, theta]} c_eff(sigma)/2, fraction 1/A0
# c_eff(sigma) = min over density inputs valid at sigma whose T-exponent c_T(1-sigma) <= 1 (absorbed by 1/U in the dyadic sum)
from fractions import Fraction as F
theta = 0.8733458561067185768
def inputs(sig, use):
    out = []
    # Montgomery (Topics Thm 12.1): (kT)^{3(1-s)/(2-s)}; T-exp <= 1 always
    out.append(('Montgomery 3/(2-s)', 3/(2-sig), 3/(2-sig)))
    if 'L' in use and sig >= 0.8: out.append(('Montgomery large-s 2/s', 2/sig, 2/sig))
    if 'H' in use: out.append(('Huxley 12/5', 12/5, 12/5))
    if 'C' in use: out.append(('Chen-Gupta-Li 2507.08296 k^{7/3(1-s)}T^{2(1-s)}+(kT)^{30/13(1-s)}', 7/3, 30/13))
    ok = [(n, ck) for (n, ck, cT) in out if cT*(1-sig) <= 1 + 1e-15]
    return min(ok, key=lambda x: x[1])
for use in ['', 'L', 'LH', 'LHC']:
    best = (0, None, None)
    N = 200000
    for i in range(N+1):
        s = 0.5 + (theta-0.5)*i/N
        n, c = inputs(s, use)
        if c > best[0]: best = (c, s, n)
    A0 = best[0]/2
    print(f"inputs={use or 'M'}: sup c_eff = {best[0]:.12f} at sigma = {best[1]:.6f} ({best[2]}); A0 = {A0:.12f}; fraction 1/A0 = {1/A0:.12f}")
print("check 2(2-theta)/3 =", 2*(2-theta)/3, " 6/7 =", 6/7, " 5/6 =", 5/6)
# interval where CGL's 7/3 is binding: 3/(2-s) >= 7/3 <=> s >= 5/7; 2/s >= 7/3 <=> s <= 6/7
print("binding interval [5/7, 6/7] =", 5/7, 6/7, "inside [1/2, theta]:", 6/7 < theta)
# T-exponent constraint for CGL: (30/13)(1-s) <= 1 <=> s >= 17/30
print("CGL usable for s >= 17/30 =", 17/30, "; Montgomery below: max 3/(2-s) on [1/2,17/30] =", 3/(2-17/30))
