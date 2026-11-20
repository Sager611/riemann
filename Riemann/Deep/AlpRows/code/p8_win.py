#!/usr/bin/env python3
# alprows P8: finite instance of the window row Win_at (Formalism/SharpRestriction/GlobalSchurBounds/QuadraticLogDerivCongr.lean:231-246),
# which carries m >= x^(1/2-eps) into RestSize (LossesSmallShares.lean:281-289, w13). Sizes as in SizesTextAt (CORE.md s5):
# n = x^(2+8e), f = 2n (Params.lean:30), m = x^(1/2-e/16), k = x^(9/2+10e), b = exp(x^e), a = x b, X' = exp(x^(e/3)), b' = a/X'.
# Checks the conjuncts with known sizes in log space; reports the least log10 x (on a grid) where all hold, for c6 in {1, 1e-3}.
import mpmath as mp
mp.mp.dps = 30
def conj(L, e, c6):   # L = ln x
    x = mp.e ** L
    lm = (mp.mpf(1) / 2 - e / 16) * L; lk = (mp.mpf(9) / 2 + 10 * e) * L; ln_ = (2 + 8 * e) * L
    lf = mp.log(2) + ln_; lXp = x ** (e / 3)
    c = {
      'm<=x^(1/2+e/2)': lm <= (mp.mpf(1) / 2 + e / 2) * L,
      'x^(1/2-e)<=m': (mp.mpf(1) / 2 - e) * L <= lm,
      'x m<=k': L + lm <= lk,
      'k<=n^2 x^(1/2-e)': lk <= 2 * ln_ + (mp.mpf(1) / 2 - e) * L,
      'k b\'<=a (k<=X\')': lk <= lXp,
      '2f(x+1)exp(-c6 x^(e/8))<=1/x': mp.log(2) + lf + mp.log(x + 1) - c6 * x ** (e / 8) <= -L,
    }
    return c
for e in (mp.mpf('0.01'), mp.mpf(1) / 60):
    for c6 in (mp.mpf(1), mp.mpf('0.001')):
        first = None
        for l10 in range(1, 200001, 1):
            L = l10 * mp.log(10)
            c = conj(L, e, c6)
            if all(c.values()):
                first = l10; break
        bad = [k for k, v in conj(mp.mpf(10) * mp.log(10), e, c6).items() if not v]
        print('eps=%s c6=%s: all listed Win_at conjuncts hold from log10 x = %s on (grid step 1); failing at x=1e10: %s' % (mp.nstr(e, 4), mp.nstr(c6, 3), first, bad))

