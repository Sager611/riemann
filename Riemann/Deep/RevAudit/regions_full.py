# Full-precision (50 digits) region coefficients of 2602.05614 via Explicit611/sanity_authors.py's re-derivation of the paper's
# formulas (M1-M7, TeX l. 590-945), then the max over a log-grid of [t0, t1] of the ratio bound/(t^{1/6} log t): <= 0.611 ?
src = open('/Users/adriansager/playground/riemann/Riemann/Deep/Explicit611/sanity_authors.py').read()
exec(src.split('print("region')[0])
for (t0, t1, c, phi, paper) in regions:
    co = corput1(t0, t1, c, phi)
    T0, T1 = mpf(t0), mpf(t1)
    best = max((ratio(co, T0 * (T1 / T0) ** (mpf(k) / 4000), phi), k) for k in range(4001))
    print("[%.3g, %.3g] a=(%s) max ratio %s at k=%d/4000 %s" % (t0, t1, ", ".join(mp.nstr(x, 9) for x in co), mp.nstr(best[0], 10), best[1],
          "EXCEEDS 0.611" if best[0] > mpf('0.611') else "ok"))
